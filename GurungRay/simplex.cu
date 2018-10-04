#include "simplex.cuh"
#include <omp.h>
#include "iostream"
#include <math.h>

/*
 * This file contains Streaming Method of model number 1 where in a single loop the below steps are done. This showed little better performance in Tesla K40m
 * for-loop
 *   1) create streams
 *   2) cudaMemcpyAsc( host to device )
 *   3) kernel call
 *   4) cudaMemcpyAsc( device to host)
 * end-loop
 *
 * Reading has been taken using this file. But for taking the reading for FLOPs, we use the other model since taking note of time for the kernel only excluding memcpy
 * operation is easier in the other model.
 */

// Convenience function for checking CUDA runtime API results
// can be wrapped around any runtime API call. No-op in release builds.
inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
	if (result != cudaSuccess) {
		fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
		assert(result == cudaSuccess);
	}
#endif
	return result;
}


//1st Method : Most Negative Value approach
__global__ void mykernel(double *S_MAT, int S_row, int S_col, double *Result, int S_N, double *R_data, int *R_index, int offset_res) {
	//int index = threadIdx.x + (blockIdx.x * blockDim.x);
	unsigned int index = offset_res + blockIdx.x;
	if (offset_res == -1)
		return;

	if (index < (offset_res + S_N)) {
		int tid;
		int i = 0; // used for for index
		unsigned long temp_index;
		unsigned long temp_index1;
		unsigned long base = index * S_row * S_col;
		unsigned int R_base = index * blockDim.x;  // blockDim.x = 96
		__shared__ bool c;
		__shared__ int row;	//pivotRow
		//		int col = 1;
		int Last_row = S_row - 2;//Amit now this should be 2nd last row

		//__shared__ double col1[1024];	//pivotColumn
		__shared__ double col1[513];	//pivotColumn Since our benchmarks has Maximum 512 threads used
		//		col1[threadIdx.x]=0;//initializing all values
		/*************/
		if (threadIdx.x == 0) {
			c = false;
			row = -1;		//pivotRow
		}
		__syncthreads();

		//Debugging after one iteration -----------------------
		/*if (threadIdx.x==0 && c==false){
			printf("\n ----- Iteration Before update operation----- \n");
			printf("Rows = %d and Cols = %d\n", S_row, S_col);
			for (int n=0;n<S_N;n++){
			unsigned int some = n * S_row * S_col;// * sizeof(double);
			for (int r=0;r<S_row;r++){
			for (int cl=0;cl<S_col;cl++){
			printf("%f  ",S_MAT[(int)(r + cl*S_row + some)]);
			}
			printf("\n");
			}
			printf("\n");
			}
			}*/
		/*		printf("\n");
			for (int s = 0; s < N_S; s++) {
			unsigned int some=M * N * s;
			for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
			std::cout<<MAT[(int)(i+ j*M +some)]<<"  ";
			}
			printf("\n");
			}
			printf("\n");
			}*/

		//----------------------------------------------------- Correct replacement


		while (!c) {
			__syncthreads(); //has no effect
			//   ***************** Get_Pivot function begins  *****************
			// ******** First Reduction Begins **********
			__shared__ int notEntered;
			//			__shared__ double minValue;
			__shared__ int newpivotcol;
			if (threadIdx.x == 0) {
				row = -1;		//pivotRow
				//				minValue = 0;
				newpivotcol = -1;
				notEntered = 1;
				c = true;
			}
			__syncthreads();	//making sure newpivotcol is initialised to -1
			int data_size = blockDim.x;
			tid = threadIdx.x;
			if (threadIdx.x >= 2 && threadIdx.x < (S_col - 1)) { //find minimum from last row leaving last column
				temp_index = Last_row + tid * S_row + base;	//avoiding re-computation
				R_data[tid + R_base] = S_MAT[temp_index];	//	g_data[i];
				//R_data[tid + R_base] = S_MAT[(Last_row + (unsigned long)(tid * S_row) + base)];	//	g_data[i];
				R_index[tid + R_base] = tid;//tid; should be the real index of the data

			}
			else {
				R_data[tid + R_base] = INT_MAX;	//	g_data[i];
				R_index[tid + R_base] = tid;	//tid;
			}
			__syncthreads();//here will have all values in shared/Global memory from 0 to BLOCK_SIZE

			tid = threadIdx.x;
			for (i = (data_size / 2); i > 0;) {
				//__syncthreads(); //testing here no effect
				if (tid < i) {
					if ((double)R_data[tid + R_base] >(double)R_data[tid + R_base + i]) { //> give correct values
						if (R_data[tid + R_base + i] <= 0) {	//produced more closed results on no 18 ::ToDo changed < to <= Makes no difference
							R_data[tid + R_base] = R_data[tid + R_base + i];//put the smaller value to left-side
							R_index[tid + R_base] = R_index[tid + R_base + i];

							int local_notEntered;
							local_notEntered = *(volatile int*)&notEntered;
							atomicCAS(&notEntered, local_notEntered, 0);
						}
					}
				}
				__syncthreads(); //Todo try removing this
				i >>= 1;
				if ((i != 1) && (i % 2) != 0) {	//if s is odd
					i = i + 1;
				}
				__syncthreads(); //Todo try removing this
			}
			//__syncthreads(); //test this here

			/*if (threadIdx.x==0){
				printf("\n ----- After ONe Iteration ----- \n");
				int some = 0 * S_row * S_col * sizeof(double);
				//for (int r=0;r<S_row;r++){
				int r=Last_row;
				for (int cl=0;cl<S_col;cl++){
				printf("%f  ",S_MAT[(int)(r + cl*S_row + some)]);
				}
				printf("\n");
				//}
				}*/


			// if (notEntered == false && tid == 2) { // tid==0 is always true if minValue is still -1 then what?
			if (threadIdx.x == 0) { // tid==0 is always true if minValue is still -1 then what?
				if (notEntered == false) {
					//					minValue = R_data[R_base];
					newpivotcol = R_index[R_base];
				}
				//	printf("newpivotcol=%d",newpivotcol);
			}
			__syncthreads(); //waiting for all threads to have same newpivotcol value
			// ********* First Reduction Ends *************

			//  ******** Second Reduction Begins **********
			if (newpivotcol == -1) {//All Threads will follow the Same path so no issue with divergence
				c = true; //can terminate
			}
			else { //if pivot column found then Find pivot row
				// ********** Second Reduction Process ******
				__shared__ double row_min;
				__shared__ int row_num;
				__shared__ int notEntered2;
				if (threadIdx.x == 0) {
					row_min = INT_MAX;
					row_num = -1;
					notEntered2 = 1;
				}
				__syncthreads();
				int k1 = 0;
				//if (threadIdx.x >= 0 && threadIdx.x < Last_row) {
				if (threadIdx.x < Last_row) {
					k1 = threadIdx.x;	//here k1 =0 to Last_row only
					//for (int k1 = 0; k1 < Last_row; k1++) {	//Last_row = (S_row - 1)
					unsigned long temp_index2 = newpivotcol * S_row + k1 + base;
					temp_index1 = k1 + (S_col - 1) * S_row + base; //avoiding re-computation
					//					__syncthreads(); //testing has effect in the result
					//if ((S_MAT[temp_index2] > (double)0.0001) && (S_MAT[temp_index1] >= (double)0)) { //changed to >=0.0001 makes no difference to the 8 working benchmarks
					if ((S_MAT[temp_index2] >(double)0.0000000001) && (S_MAT[temp_index1] >= (double)0)) { //No Cycling both for benchmark 15 and 4 //ToDo:: changed to >=0.0001
						//if ((S_MAT[temp_index2] > 0.0001) && (S_MAT[temp_index1] >= 0)) { //No Cycling both for benchmark 15 and 4 //ToDo:: changed to >=0.0001
						R_data[k1 + R_base] = (double)S_MAT[temp_index1] / (double)S_MAT[temp_index2]; //b_i / S_MAT[pivotcol]
						R_index[k1 + R_base] = k1;
						//-------------------------------------------------
						//Since there exists some feasible value which may be in the 1st location
						int local_notEntered2;
						local_notEntered2 = *(volatile int*)&notEntered2;
						atomicCAS(&notEntered2, local_notEntered2, 0);
						//-------------------------------------------------
					}
					else {
						R_data[k1 + R_base] = INT_MAX; //to make the array size equal
						R_index[k1 + R_base] = k1; //to make the array size equal
					}
					//	__syncthreads();//testing here no effect
				}
				else { //remaining threads above Last_row(including) upto Block_Size
					k1 = threadIdx.x;
					R_data[k1 + R_base] = INT_MAX; //to make the array size equal
					R_index[k1 + R_base] = k1; //to make the array size equal
				}
				__syncthreads(); //Verified All data and index stored correctly with index as threadIdx.x
				//Now find the minValue and its index from R_data and R_index using Reduction

				int data_size2 = blockDim.x; //Now it is Block_Size
				// ************************* Second Reduction on R_data and R_index ************************
				//	if (threadIdx.x >= 0 && threadIdx.x < Last_row) {	//Now for all threads
				tid = threadIdx.x;
				//			__syncthreads(); //testing
				for (int s = (data_size2 / 2); s > 0;) {
					//__syncthreads();//testing no effect
					if (tid < s) {
						int indexValue2 = tid + R_base;
						if (R_data[indexValue2] >= R_data[indexValue2 + s]) { //changed >= to > ToDo:: Accordingly Fix other bug. Back to >= fixed no 1
							R_data[indexValue2] = R_data[indexValue2 + s];	//For arranging in ascending order we better swap the value instead of only replacing
							R_index[indexValue2] = R_index[indexValue2 + s];
						}
					}
					__syncthreads();
					s >>= 1;
					if ((s != 1) && (s % 2) != 0) {	//if s is odd
						s = s + 1;
					}
					__syncthreads(); //Todo:: try removing this. This creates unpredictable behaviour
				}
				//	__syncthreads(); //test this here no effect
				//if (notEntered2 == false && tid == 0) {
				if (tid == 0) {
					if (notEntered2 == false) {	//if at least once swaped the flag 'notEntered2' will be equal to 0
						row_min = R_data[R_base];
						row_num = R_index[R_base];
					}
					//	printf("newpivotrow=%d\n",row_num);
				}
				__syncthreads(); // Looks like this can be skipped 	//Todo:: try removing this.
				// ********** Second Reduction on R_data and R_index ******
				if (threadIdx.x == 0) {
					//					__syncthreads(); //testing has effect serious/worse effect on all benchmarks
					if (row_min == INT_MAX) {
						row = -1;
					}
					//__syncthreads();//testing has no effect
					if ((row_min != INT_MAX) && (row_num != -1)) {
						row = row_num;
					}
				}
				//__syncthreads(); // this has no effect on all benchmarks
			} //end of else of newpivotcol == -1
			__syncthreads(); // Looks like this can be skipped but here we have row synchronized
			//  ******** Second Reduction Ends **********
			//   ***************** Get_Pivot function ends  *****************

			//*******************************************************************************************************
			//			col = newpivotcol;
			//			__syncthreads(); //has no effect
			//if ((row > -1) && (col != -1)) {	//some candidate leaving variable or pivot row have been found and also some pivot column found
			if ((row > -1) && (newpivotcol != -1)) {	//some candidate leaving variable or pivot row have been found and also some pivot column found
				//*******************************************************************************************************
				tid = threadIdx.x;
				//		if (threadIdx.x >= 0 && threadIdx.x < S_row) {
				if (threadIdx.x < S_row) {
					//for (int i = 0; i < S_row; i++) {	//Data Parallel section 2
					//col1[tid] = (double)S_MAT[(unsigned long)((tid + col * S_row) + base)];//keeping the old pivotcol coeff
					col1[tid] = (double)S_MAT[(unsigned long)((tid + newpivotcol * S_row) + base)];//keeping the old pivotcol coeff
				}	//Data Parallel section 2 done
				__syncthreads(); // has no effect can be removed
				//*******************************************************************************************************
				//			if (threadIdx.x==0){
				unsigned long temp_row_base = row + base;
				//		__syncthreads();		//testing here
				//S_MAT[(unsigned long)(temp_row_base + 1 * S_row)] = S_MAT[(unsigned long)((S_row - 1) + (unsigned long)(col * S_row) + base)]; //column 1 replaced by the new Last row containing objective coefficient
				S_MAT[(unsigned long)(temp_row_base + 1 * S_row)] = S_MAT[(unsigned long)((S_row - 1) + (unsigned long)(newpivotcol * S_row) + base)]; //column 1 replaced by the new Last row containing objective coefficient
				//		__syncthreads();		//testing here
				//S_MAT[temp_row_base] = col - 1; //replacing entering variable index in leaving variable index (1-based indexing)
				S_MAT[temp_row_base] = newpivotcol - 1; //replacing entering variable index in leaving variable index (1-based indexing)
				//			}
				__syncthreads();

				//*******************************************************************************************************

				tid = threadIdx.x;
				//		__syncthreads();		//testing here
				if (threadIdx.x >= 2 && threadIdx.x < S_col) {
					//for (int j = 2; j < S_col; j++){		//Data Parallel section 3
					unsigned long row_base = row + base;	//avoiding re-computation
					temp_index = row_base + (tid * S_row);//avoiding re-computation
					S_MAT[temp_index] = (double)S_MAT[temp_index] / (double)col1[row];//updating pivot row by dividing, current pivot row value by (pivot element)
				}		//Data Parallel section 3 done
				__syncthreads();
				//*******************************************************************************************************
				tid = threadIdx.x;
				//if ((threadIdx.x >= 0) && (threadIdx.x < (S_row - 1))) { //updating all rows
				if (threadIdx.x < (S_row - 1)) { //updating all rows
					if (tid != row) {	//this if outside of for-loop give better performance
						//for (int i = 0; i < S_row; i++) {	//Data parallel section 4
						for (i = 2; i < S_col; i++) {
							temp_index1 = (unsigned long)(i * S_row + base);
							temp_index = (unsigned long)(tid + temp_index1);
							double zeroTemp = 0.0;
							zeroTemp = col1[tid] * S_MAT[(unsigned long)(row + temp_index1)];
							S_MAT[temp_index] = S_MAT[temp_index] - zeroTemp;

							//S_MAT[(unsigned long)(tid + i * S_row + base)] = S_MAT[(unsigned long)(tid + i * S_row + base)] - (col1[tid] * S_MAT[(unsigned long)(row + i * S_row + base)]);

						}
					}
				}	//Data Parallel section 4 done
				__syncthreads();
				//*******************************************************************************************************

				/*//tid = threadIdx.x;
				//if (threadIdx.x == 0) {
				//for (i = 2; i < (S_col - 1); i++) {
				i=threadIdx.x;
				if (threadIdx.x >= 2 && threadIdx.x < (S_col - 1)){
				if ((double)S_MAT[(Last_row + i * S_row) + base] < (double)0.0){//working
				c = false;
				//break;
				}
				}
				//}*/

				if (threadIdx.x == 0) {
					for (i = 2; i < (S_col - 1); i++) {
						if ((double)S_MAT[(Last_row + i * S_row) + base] < (double)0.0){ //working
							c = false;
							break;
						}
					}
				}
				__syncthreads();

			} /*else if (row == -1) { //No candidate leaving row have been found so remember this pivot column (and try next iteration although not mentioned in Algorithm)
				//This is actually the situation of UNBOUNDEDNESS
				c = true;
				}*/

			//Debugging after one iteration -----------------------
			/*if (threadIdx.x==0){
				printf("\n ----- After ONe Iteration ----- \n");
				for (int n=0;n<S_N;n++){
				int some = n * S_row * S_col * sizeof(double);
				for (int r=0;r<S_row;r++){
				for (int cl=0;cl<S_col;cl++){
				printf("%f  ",S_MAT[(int)(r + cl*S_row + some)]);
				}
				printf("\n");
				}
				printf("\n");
				}
				}*/
			/*				if (threadIdx.x==0){
								printf("\n ----- After ONe Iteration ----- \n");
								int some = 0 * S_row * S_col * sizeof(double);
								//for (int r=0;r<S_row;r++){
								int r=Last_row;
								for (int cl=0;cl<S_col;cl++){
								printf("%f  ",S_MAT[(int)(r + cl*S_row + some)]);
								}
								printf("\n");
								//}
								}*/
			//-----------------------------------------------------
			//__syncthreads(); //testing here has no effect
		} //end of while
		__syncthreads();

		if (threadIdx.x == 0) {
			Result[index] = S_MAT[(Last_row + (S_col - 1) * S_row) + base];
			// printf("Result[index] =%f\n",Result[index]);
			//printf("Size of double =%d\n",sizeof(double));
		}
	}
}

__host__ Simplex::Simplex(unsigned int N_S, unsigned int multipleBatch) {
	//Does NOT require multiple batch solving

	//std::cout << "\n*****GPU class instaintiated*****\n";


	number_of_LPs = N_S;
	M = 0;
	N = 0;
	c = 0;
	No_c = 0;
	multipleBatchOn = multipleBatch; //Flag to determine if multiple batched is called
	//R = (double*) calloc(N_S,sizeof(double));
}

//get status of particular simplex
__host__ int Simplex::getStatus(int n) {
	int s;
	for (int i = 0; i < C.size1(); i++) {
		if (i == (n - 1)) {
			if (R[i] == -1) {
				s = 6;	// 6 = Simplex Is Unbounded
			}
			else if (R[i] > 0) {
				s = 2;	// 2= Simplex has feasible and Optimal solution
			}
		}
	}
	return s;

}	//get status of particular simplex

//get the No of simplex the object is ruuning on GPU
__host__ int Simplex::getNo_OF_Simplx() {
	return C.size1();
}	//get the No of simplex the object is ruuning on GPU

//get the result of all simplex
__host__ std::vector<double> Simplex::getResultAll() {

	std::vector<double> Res(C.size1());
	if (multipleBatchOn >= 1){
		return allResult;
	}
	else {
		for (int i = 0; i < C.size1(); i++) {
			Res[i] = R[i];
			//std::cout<<"  R[i]="<<R[i];
		}
	}
	return Res;
}

//get the result of all simplex

__host__ double Simplex::getResult(int n) {
	// get result of particular simplex
	float r;
	for (int i = 0; i < C.size1(); i++) {
		if (i == (n - 1)) {
			r = R[i];
		}
	}
	return r;
}	// get result of particular simplex

__host__ std::vector<int> Simplex::getStatusAll() {

	std::vector<int> Status(C.size1());
	for (int i = 0; i < C.size1(); i++) {
		if (R[i] == -1)
			Status[i] = 6;
		else
			Status[i] = 2;
	}
	return Status;
}	//get the status of all simplex

__host__ void Simplex::setConstratint(math::matrix<double> A, std::vector<double> B) {
	unsigned int N_S = number_of_LPs;
	orig_CoefficientMatrix = A;
	BoundValue = B;
	unsigned int No_O = A.size2();
	unsigned int No_C = A.size1();
	//M = No_C + 1;
	M = No_C + 2;//Extra row for coefficient of objective function
	N = No_O + 3 + No_C;//original variables + slack + 3 extra(index,pivot-col,b_i); artificial is not included now/here
	c = 1 + No_O;

	//	MAT = (double *) calloc(N_S * M * N, sizeof(double));


	unsigned long long memSize = N_S * M * N * sizeof(double);
	cudaError_t err;
	err = cudaMallocHost(&MAT, memSize);//Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
	//	printf("CUDA cudaMallocHost-- MAT: %s\n", cudaGetErrorString(err));
	err = cudaMemset(MAT, 0, memSize);	//initializing all elements to zero ... Now try with Asynchronous Command
	//	err = cudaMemsetAsync(MAT, 0, memSize);// by default the last argument stream=0
	//	printf("CUDA MemsetAsync to 0-- MAT: %s\n", cudaGetErrorString(err));
	//cudaDeviceSynchronize();
	/*
	 * Simplex tableau Re-Structure Amit :: Note The variables are implemented as 1-based indexing
	 * row-size= (m-constraints + 1-row Z and Optimal Solution value + 1-row for copy of coefficient of Objective function) = m+2
	 * column-size = (n-original-variables + m-slack-variables + a-artificial-variables+ 3 (1 for index, 1 for coefficient of the basic variables, 1 for bounds b_i's)
	 * column 0: index of basic/slack variables
	 * column 1: coefficient of basic variables
	 * column 2 to n: coefficients of variables (non-basic) --implemented as 2 to No_O+2 where No_O is the size of original variables
	 * column (n+1) to (n+m): includes slack variables --implemented as (No_O+2 +1) to (No_O+2+No_C) where No_C is the number of constraints comprising slack variables
	 * column (n+m+1) to (n+m+a): includes artificial variables --implemented after (No_O+2+No_C + 1) to <(N-1) where N is the total size of columns
	 * column last column (N-1) : bounds b_i
	 *
	 * row 0 through m: contains the coefficient for Simplex method algorithm
	 * row (m+1), the 2nd last row: contains the values of the operations (Cj - Zj) of each iterations starting from column 2 through (last - 1) and the last column contains the
	 * value of the optimal solution of each iterations.
	 * Last row (m+2): contains a copy of the coefficients of the objective function required for Simplex Algorithm. Column 2 through (N-1) is used to store these values
	 * ** NB: In phase-I : Artificial variables contains -1 and all values are 0. But in phase-II artificial variables are eliminated and the original coefficients are replaced
	 * 		  with original variables having their respective values and 0 for slack variables.
	 *
	 */

	//printf("\nBefor Setting A, b multipleBatchOn=%d",multipleBatchOn);
#pragma omp parallel for
	for (long s = 0; s < N_S; s++) {
		unsigned long long some = M * N * s;
		for (unsigned int i = 0; i < (M - 2); i++) {
			for (unsigned int j = 0; j < N; j++) {
				if (j == 0) {	//index of basic/slack variables
					/*//if (multipleBatchOn > 1){
						std::cout<<"(i,j)=("<<i<<","<<j<<") index="<<(unsigned long long)(i  + some)<<"---";
						std::cout<<"c+i="<<c+i<<" --- "<<std::endl;
						MAT[(i  + some)] = c + i;
						//				std::cout<<"MAT="<<MAT[(i  + some)]<<" ";
						//}*/
					MAT[(unsigned long long)(i + some)] = c + i;
					//MAT[(unsigned long long)((i + j * M) + some)] = c + i;
				}
				else if (j > 1) { //excluding 'column 1' from loop
					if (j < (No_O + 2)) {	// coefficients of the variables (a.k.a. non-basic)
						//std::cout<<"(i,j)=("<<i<<","<<j<<") index="<<((i + j * M) + some)<<"---";
						//std::cout<<"A="<<A(i,j - 2)<<" ";
						//std::cout<<"MAT="<<MAT[(unsigned long long) ((i + j * M) + some)]<<std::endl;
						MAT[(unsigned long long) ((i + j * M) + some)] = A(i, j - 2);
					}
					else if (j == (N - 1)) { //last column stores the bounds b_i
						MAT[(unsigned long long) ((i + j * M) + some)] = B[i];
					}
					else if (j < (N - 1)) { //includes slack variables
						MAT[(unsigned long long) ((i + (No_O + 2 + i) * M) + some)] = 1;
					}
				}
			}
		}
	}

	//Debugging
	/*		printf("\n");
		for (int s = 0; s < N_S; s++) {
		unsigned int some=M * N * s;
		for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
		std::cout<<MAT[(int)(i+ j*M +some)]<<"  ";
		}
		printf("\n");
		}
		printf("\n");
		}*/

	//	printf("\nSetting A, b done\n");
}

__host__ void Simplex::ComputeLP(math::matrix<double> &C1, unsigned int number_of_streams) {

	cudaError_t err;
	unsigned int threads_per_block;	//Maximum threads depends on CC 1.x =512 2.x and > = 1024
	unsigned int number_of_blocks;//depends on our requirements (better to be much more than the number of SMs)
	int device;
	cudaDeviceProp props;
	cudaGetDevice(&device);
	cudaGetDeviceProperties(&props, device);
	int No_C = orig_CoefficientMatrix.size1();
	C = math::matrix<double>(C1);
	unsigned int N_S = C.size1();
	int No_O = C.size2();
	M = No_C + 2, N = No_O + 3 + No_C;// M is now + 2 instead of + 1
	int N_C = No_C;


	unsigned long memSize = N_S * sizeof(double);
	err = cudaMallocHost((void**)&R, memSize);	//PINNED Memory	 //cudaMallocHost((void**)&a, bytes) );      // host pinned
	//printf("CUDA cudaMallocHost-- R: %s\n", cudaGetErrorString(err));


	std::vector<int> rem;
	for (int i = 0; i < N_C; i++) {
		//std::cout<<BoundValue[i]<<"\n";
		if (BoundValue[i] < 0) {
			rem.push_back(i);
		}
	}

	//	std::cout<<"Number of Artificial Variables = "<< rem.size()<<"\n";
	int nc = N + rem.size();

	threads_per_block = 32 * (nc / 32) + 32; //if count equal 0 than nc=N so works for all Model
	if (threads_per_block > props.maxThreadsPerBlock) //Assuming maximum threads supported by CC is 1024
		threads_per_block = props.maxThreadsPerBlock;
	std::cout << "Number of threads per LP = " << threads_per_block << std::endl;
	int *R_index;	//reduction data
	double *R_data;	//reduction index
	err = cudaMalloc((void **)&R_data, C1.size1() * threads_per_block * sizeof(double));//C1.size1() * 96 being the maximum threads
	//printf("CUDA malloc R_data: %s\n", cudaGetErrorString(err));
	err = cudaMalloc((void **)&R_index, C1.size1() * threads_per_block * sizeof(int));//C1.size1() being the number of LPs
	//printf("CUDA malloc R_index: %s\n", cudaGetErrorString(err));
	err = cudaMalloc((void **)&G_R, N_S * sizeof(double));//Doing it here for the First Time
	//printf("CUDA malloc G_R: %s\n", cudaGetErrorString(err));

	//std::cout<<"Number of Artificial Variables = "<< rem.size()<<"\n";

	//std::cout << "Number of threads per block = " << threads_per_block << "\n";

	if (rem.size() > 0) {	//Two-Phase Simplex Algorithm for LPs with infeasible Initial Basic Solution\n";

		//	std::cout << "Two-Phase Simplex Algorithm Running\n";

		unsigned long long memSize = N_S * M * nc * sizeof(double);
		cudaError_t err;
		err = cudaMallocHost(&N_MAT, memSize);//Pinned memory Syntax:: cudaMallocHost(&h_ptr,bytes);
		//printf("CUDA cudaMallocHost-- MAT: %s\n", cudaGetErrorString(err));
		cudaMemset(N_MAT, 0, memSize);	//initializing all elements to zero
		//N_MAT = (double *) calloc(N_S * M * nc, sizeof(double)); //initialized to zero (this tableau include artificial variables)

		c = 1 + No_O;	//for correct assigning of slack variable index
		//int s;

#pragma omp parallel for
		for (int s = 0; s < N_S; s++) {
			unsigned long long some = M * N * s; //base address for each LP
			for (int j = 2; j < (No_O + 2); j++) {
				MAT[(unsigned long long) (((M - 2) + j * M) + some)] = C(s, j - 2); //Last row (M-1)  Amit::removed negative Now modified to M-2 the 2nd last row
				//Now keep a copy of coefficients of the objective function
				MAT[(unsigned long long) (((M - 1) + j * M) + some)] = C(s, j - 2); //slack is already zero as initialized
				//}
			}
		}

#pragma omp parallel for
		for (int i = 0; i < N_S; i++) {
			unsigned long long base = i * M * N;//base of every LP in MAT tableau
			unsigned long long basen = i * M * nc;//base of every LP in N_MAT tableau (this include artificial variables)
			for (unsigned int j = 0; j < M; j++) {	//from every row/constraints
				//base=i*M*N;
				N_MAT[(unsigned long long)(j + ((nc - 1) * M) + basen)] = MAT[(unsigned long long)(j + ((N - 1) * M) + base)]; // N_MAT[lastCol] = MAT[lastCol]
			}
			for (unsigned int j = 2; j < (nc - 1); j++) {	//from every column
				if (j >= ((No_O + 3 + No_C) - 1))//for Phase-I //original and slack variables will remain 0//N_MAT[(int) ((M-1) + j * M + some)] = 0;
					N_MAT[(unsigned long long)((M - 1) + j * M + basen)] = -1; //artificial variable
			}
		}

		//Creating Artificial Variables
#pragma omp parallel for
		for (int k = 0; k < N_S; k++) {
			bool once = false;
			int artif = 0, ch;
			unsigned long long base = k * M * N;//base of every LP in MAT tableau
			unsigned long long basen = k * M * nc;//base of every LP in N_MAT tableau (this include artificial variables)
			for (unsigned int i = 0; i < (M - 1); i++) { //for every row including 2nd last row of the tableau :: leave last row
				ch = 0;
				for (unsigned int j = 0; j < nc; j++) {  //for every column of the MAT and N_MAT tableau
					if (MAT[(unsigned long long)(i + ((N - 1) * M) + base)] < 0) { //this indicates the negative b_i from the MAT tableau
						if ((j >= (N - 1)) && (j < (nc - 1))) { //this indicate all columns that represent artificial variables
							if (!ch) {
								float v = N_MAT[(unsigned long long)((i - 1) + (j * M) + basen)]; //why (i-1)?
								if ((once) && (v == 1)){ // computing v is meaningless since once is false and not made true anywhere so this will always be false
									N_MAT[(unsigned long long)((i + (j + 1) * M) + basen)] = 1;//so this block will never be executed ToDo:: can be skipped
								}
								else {
									N_MAT[(unsigned long long)((i + j * M) + basen)] = 1;
								}
								ch = 1;	//this will allow populating 1, diagonally in artificial variables
							}
						}
						else if (j == (nc - 1)) { //this indicate the last column of N_MAT tableau which is b_i's
							N_MAT[(unsigned long long)((i + j * M) + basen)] = -1 * N_MAT[(unsigned long long)((i + j * M) + basen)]; //negating b_i's
						}
						else if (j == 1) { //the extra temporary working column
							N_MAT[(unsigned long long)((i + j * M) + basen)] = -1; //why populate -1 only for negative b_i's in the extra column? May be used as coefficient of artificial variables
						}
						else if (j == 0) { //first index column
							//NOTE: Binayak used table index as variable indexing so it is 1-based Indexing
							//ToDo:: Amit detected Bugged here in index computation for Artificial variables
							//N_MAT[((i + j * M)) + basen] = (N + i)-2;//computes the index of artificial variables as n+m+a where size of vars, slacks and artificial are n,m and a respectively
							N_MAT[(unsigned long long)((i + j * M) + basen)] = (N + artif) - 2;//increase index only when found artificial variable and not for every row i.
							//	std::cout<<" artif = " <<artif;
							artif++;//increase for next artificial variable found
							//std::cout<<" (N + i)-2 = " <<(N + i)-2;
						}
						else if (j > 1) { //negated all variables and slacks (only non-basic excluding artificial)
							N_MAT[(unsigned long long)((i + j * M) + basen)] = -1 * (MAT[(unsigned long long)((i + j * M) + base)]);
						}
					}
					else if ((i != (M - 2)) && (j < (N - 1))) { //except last row and last column of MAT i.e. b_i's
						N_MAT[(unsigned long long)((i + j * M) + basen)] = MAT[(unsigned long long)((i + j * M) + base)];//copy into N_MAT as it is
					}
					else if (i == (M - 2)) {
						if ((j >= (N - 1)) && (j < (nc - 1))) {
							N_MAT[(unsigned long long)((i + j * M) + basen)] = -1; //ALL artificial variable coefficient is assigned -1
						}
					}
				}
			}
		}

		//Creation of Last Row or Z-Value(Zj-Cj)
#pragma omp parallel for
		for (int k = 0; k < N_S; k++) {
			//int sum = 0;
			//base = k * M * N;
			unsigned long long basen = k * M * nc;
			for (unsigned int k1 = 2; k1 < nc; k1++) {//for all columns upto b_i from column 2
				double sum = 0.0; //reset for every column k1 (objective function value)
				for (unsigned int j = 0; j < (M - 2); j++) { //for all rows except the 2nd last row for which this computation is performed and also last row
					sum = sum + (N_MAT[(unsigned long long)((j + k1 * M) + basen)] * N_MAT[(unsigned long long)((j + 1 * M) + basen)]); // column 1 currently contains -1, the coefficient of artificial variables
				}
				//std::cout << sum << "-"	<< N_MAT[((M - 1) + k1 * M) + basen];
				N_MAT[(unsigned long long)(((M - 2) + k1 * M) + basen)] = sum - N_MAT[(unsigned long long)(((M - 2) + k1 * M) + basen)]; //formula Zj - Cj
			}
		}
		//cudaEvent_t start, stop;
		//		std::cout << "Before Kernel Called 1\n";
		err = cudaMalloc((void **)&G_MAT, (N_S * M * nc * sizeof(double)));
		//printf("CUDA Malloc G_MAT : %s\n", cudaGetErrorString(err));



		/*
			//Debugging ---------------------------------------
			printf("\n******************* N_MAT tableau *******************\n");
			for (int s = 0; s < N_S; s++) {
			unsigned int some=M * nc * s;
			for (int i = 0; i < M; i++) {
			for (int j = 0; j < nc; j++) {
			std::cout<<N_MAT[(int)(i+ j*M +some)]<<"  ";
			}
			printf("\n");
			}
			printf("\n");
			}
			// ---------------------------------------
			*/

		// *********************** Stream Processing Begins *******************
		//Using Asynchronous Memory copy:: needs //MAT to be a PINNED memory
		int num_streams = number_of_streams;//number of streams desired to create ::Note check for odd numbers
		unsigned long Each_LP_size = M * nc;	//  sizeof(float);
		unsigned long num_LPs_perStream;
		bool equal_stream = true;
		if (N_S % num_streams == 0) {
			num_LPs_perStream = (N_S / num_streams);
			equal_stream = true;
		}
		else {
			num_LPs_perStream = (N_S / num_streams);//last stream will not be of the same size
			num_streams = num_streams + 1;//one extra streams.where nos of LPs to be solved will be less;
			equal_stream = false;
		}
		unsigned long lastBlock_size;
		if (equal_stream == false) {
			lastBlock_size = N_S - (N_S / (num_streams - 1)) * (num_streams - 1);//LAST Stream Size
			//std::cout << "\nAmit Last Block size (LPs is )= " << lastBlock_size<< std::endl;
		}

		//std::cout << "\nNumber of LPs_perStream = " << num_LPs_perStream <<"  Total Streams="<<num_streams<< std::endl;
		// WINDOWS FIX allocate streams
		cudaStream_t *stream;
		stream = (cudaStream_t*)malloc(sizeof(cudaStream_t)*num_streams);

		cudaError_t result;
		// **************** Creation of Streams ****************
		for (int i = 0; i < num_streams; i++) {
			result = cudaStreamCreate(&stream[i]);	//Creation of Streams
			unsigned long long offset = i * Each_LP_size * num_LPs_perStream;//for memory copy the starting address
			unsigned long offset_res = i * num_LPs_perStream;	//for result here offset_res is a pointer to the LP number
			//Stream -- memcopy Host to Device
			if (equal_stream == false && i == (num_streams - 1)) {//last stream
				cudaMemcpyAsync(&G_MAT[offset], &N_MAT[offset], (lastBlock_size * M * nc * sizeof(double)), cudaMemcpyHostToDevice, stream[i]);
				//printf("CUDA memcopyAsync G_MAT last-stream: %s\n", cudaGetErrorString(err));
				mykernel << <lastBlock_size, threads_per_block, 0, stream[i] >> >(G_MAT, M, nc, G_R, lastBlock_size, R_data, R_index, offset_res);
				cudaMemcpyAsync(&R[offset_res], &G_R[offset_res], (lastBlock_size * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
				cudaMemcpyAsync(&N_MAT[offset], &G_MAT[offset], (lastBlock_size * M * nc * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
			}
			else {
				cudaMemcpyAsync(&G_MAT[offset], &N_MAT[offset], (num_LPs_perStream * M * nc * sizeof(double)), cudaMemcpyHostToDevice, stream[i]);
				//printf("CUDA memcopyAsync G_MAT: %s\n", cudaGetErrorString(err));
				mykernel << <num_LPs_perStream, threads_per_block, 0, stream[i] >> >(G_MAT, M, nc, G_R, num_LPs_perStream, R_data, R_index, offset_res);
				cudaMemcpyAsync(&R[offset_res], &G_R[offset_res], (num_LPs_perStream * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
				cudaMemcpyAsync(&N_MAT[offset], &G_MAT[offset], (num_LPs_perStream * M * nc * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
			}
			//Stream -- Kernel
			//Stream -- memcopy Device to Host
		}

		cudaDeviceSynchronize(); //Required before N_MAT is accessed for further computation

#pragma omp parallel for
		for (int i = 0; i < N_S; i++) {
			unsigned long long base = i * M * N;
			unsigned long long basen = i * M * nc;
			for (unsigned int j = 0; j < M; j++) { //for every row
				for (unsigned int k = 0; k < N; k++) { //for each column in MAT
					if (N_MAT[(unsigned long long)(j + 0 * M + basen)] == (k + 1)) { //column 0 i.e. (row,0) has index of variables starting from (1 to n+m+a) so k starts from (1 to N, N=n+m+a+3)
						//So this condition will be met every value of k. if variable indexing is correctly assigned (Note Artificial is in-correct)
						N_MAT[(unsigned long long)(j + 1 * M + basen)] = MAT[(unsigned long long)((M - 1) + (2 + k) * M + base)]; //in column 1 of N_MAT replacing original problem's objective coefficients
					}
				}
			}
		}

		bool feasible = true;
#pragma omp parallel for
		for (int s = 0; s < N_S; s++) {
			if (round(R[s]) == 0.000000) {
				//if ((roundf(R[s]/10000)*10000) == 0) {
				//std::cout<<"\nInside ** IF*** Result = "<<R[0]<<"\n";
				//int sum = 0;
				unsigned long long base = s * M * N;
				unsigned long long basen = s * M * nc;
				for (unsigned int i = 0; i < N; i++) { //for each column i
					double sum = 0;
					for (unsigned int j = 0; j < (M - 1); j++) { //for every row j except the new last row of obj. coefficients
						if ((j < (M - 2))) { //except the last row
							if (i != (N - 1)) { //except the last column ie b_i's
								//std::cout<<N_MAT[(j+(i*M))+basen]<<"*"<<N_MAT[(j+(1*M))+basen]<<std::endl;
								sum = sum + (N_MAT[(unsigned long long)((j + (i * M)) + basen)] * N_MAT[(unsigned long long)((j + (1 * M)) + basen)]); //sum = sum + N_MAT[row,1] * N_MAT[row,i]
								MAT[(unsigned long long)((j + (i * M)) + base)] = N_MAT[(unsigned long long)((j + (i * M)) + basen)]; //copy data from N_MAT to MAT for every column i, as row by row (row is j)
							}
							else if (i == N - 1) { //for the last column ie b_i's
								sum = sum + (N_MAT[(unsigned long long)((j + (nc - 1) * M) + basen)] * N_MAT[(unsigned long long)((j + (1 * M)) + basen)]); //sum is for objective value column, the product of b_i's and coefficient's in column 1
								MAT[(unsigned long long)((j + (i * M)) + base)] = N_MAT[(unsigned long long)((j + (nc - 1) * M) + basen)]; //copy data of b_i's from N_MAT to MAT
							}
						}
						//if (j == (M - 1)) { //for the last row
						if (j == (M - 2)) { //for the 2nd last row
							if (i > 1) { // for all column from variables to slack variables excluding artificial variables
								//std::cout<<sum<<" And "<<MAT[(j+(i*M))+base]<<std::endl;
								//MAT[(j + (i * M)) + base] = MAT[(j + (i * M)) + base] + (-1) * sum; // Zj = Zj - Cj
								MAT[(unsigned long long)((j + (i * M)) + base)] = sum + (-1 * MAT[(unsigned long long)((j + (i * M)) + base)]); // Zj = Zj - Cj  Amit:: Corrected
							}
						}
					}
				}
			}
			else{
				//std::cout << "\n\nThe problem is Infeasible !!!\n";
				//std::cout << "\n Value = " << R[s] << std::endl;
				feasible = false; //any one instance
			}
		}
		cudaFree(G_MAT);
		//	std::cout<<"\n\nThe problem is Feasible !!!\n";
		if (feasible){	//Second Section of Streaming for Two-phase Simplex Method

			//now stream this section
			cudaMalloc((void **)&G_MAT, (N_S * M * N * sizeof(double)));
			// *********************** Stream Processing Begins *******************
			//Using Asynchronous Memory copy:: needs //MAT to be a PINNED memory
			int num_streams = number_of_streams;//number of streams desired to create ::Note check for odd numbers
			unsigned int Each_LP_size = M * N;	//  sizeof(float);
			unsigned int num_LPs_perStream;
			bool equal_stream = true;
			if (N_S % num_streams == 0) {
				num_LPs_perStream = (N_S / num_streams);
				equal_stream = true;
			}
			else {
				num_LPs_perStream = (N_S / num_streams);//last stream will not be of the same size
				num_streams = num_streams + 1;//one extra streams.where nos of LPs to be solved will be less;
				equal_stream = false;
			}
			unsigned long lastBlock_size = 0;
			if (equal_stream == false) {
				lastBlock_size = N_S - (N_S / (num_streams - 1)) * (num_streams - 1);//LAST Stream Size
				//std::cout << "\nAmit Last Block size (LPs is )= " << lastBlock_size<< std::endl;
			}
			//std::cout << "\nNumber of LPs_perStream = " << num_LPs_perStream <<"  Total Streams="<<num_streams<< std::endl;

			//WINDOWS FIX re-allocate treams
			free(stream);
			stream = (cudaStream_t*)malloc(sizeof(cudaStream_t)*num_streams);
			
			cudaError_t result;
			// **************** Creation of Streams ****************
			for (unsigned int i = 0; i < num_streams; i++) {
				result = cudaStreamCreate(&stream[i]);
				unsigned long long offset = i * Each_LP_size * num_LPs_perStream;//for memory copy the starting address
				unsigned long offset_res = i * num_LPs_perStream;	//for result here offset_res is a pointer to the LP number
				if (equal_stream == false && i == (num_streams - 1)) { //last stream
					cudaMemcpyAsync(&G_MAT[offset], &MAT[offset], (lastBlock_size * M * N * sizeof(double)), cudaMemcpyHostToDevice, stream[i]);
					//printf("CUDA memcopyAsync G_MAT last-stream: %s\n", cudaGetErrorString(err));
					mykernel << <lastBlock_size, threads_per_block, 0, stream[i] >> >(G_MAT, M, N, G_R, lastBlock_size, R_data, R_index, offset_res);
					cudaMemcpyAsync(&R[offset_res], &G_R[offset_res], (lastBlock_size * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
					//printf("CUDA memcopyAsync G_R last-stream: %s\n", cudaGetErrorString(err));
				}
				else {
					cudaMemcpyAsync(&G_MAT[offset], &MAT[offset], (num_LPs_perStream * M * N * sizeof(double)), cudaMemcpyHostToDevice, stream[i]);
					//printf("CUDA memcopyAsync G_MAT: %s\n", cudaGetErrorString(err));
					mykernel << <num_LPs_perStream, threads_per_block, 0, stream[i] >> >(G_MAT, M, N, G_R, num_LPs_perStream, R_data, R_index, offset_res);
					cudaMemcpyAsync(&R[offset_res], &G_R[offset_res], (num_LPs_perStream * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
					//printf("CUDA memcopyAsync G_R : %s\n", cudaGetErrorString(err));
				}
			}

		}//end of two-phase
		free(stream);

	}
	else {	//Simplex Algorithm for LPs with feasible Initial Basic Solution\n";


#pragma omp parallel for
		for (int s = 0; s < N_S; s++) {
			unsigned long some = M * N * s; //base address for each LP
			for (unsigned int j = 2; j < (No_O + 2); j++) { //Amit::Infact can be < (No_O+2) for Optimization
				//if (j < 2 + No_O) { //assigning objective coefficients of variables only
				//MAT[(int) (((M-1) + j * M) + some)] = -C(s, j - 2); //Last row (M-1)
				MAT[(unsigned long long)(((M - 2) + j * M) + some)] = -C(s, j - 2); //2nd last row (Zj-Cj) initially Zj = 0
				//Now keep a copy of coefficients of the objective function
				MAT[(unsigned long long)(((M - 1) + j * M) + some)] = C(s, j - 2); //original coefficients of the objective function
			}
		}

		//std::cout<<"Number of Artificial Variables = "<< rem.size()<<"\n";
		//cudaEvent_t start, stop;

		err = cudaMalloc((void **)&G_MAT, (N_S * M * N * sizeof(double)));
		//printf("CUDA Malloc G_MAT : %s\n", cudaGetErrorString(err));

		/*
				//Debugging ---------------------------------------
				printf("\n******************* MAT tableau *******************\n");
				for (int s = 0; s < N_S; s++) {
				unsigned int some=M * N * s;
				for (int i = 0; i < M; i++) {
				for (int j = 0; j < N; j++) {
				std::cout<<MAT[(int)(i+ j*M +some)]<<"  ";
				}
				printf("\n");
				}
				printf("\n");
				}
				// ---------------------------------------
				*/

		// *********************** Stream Processing Begins *******************
		//Using Asynchronous Memory copy:: needs //MAT to be a PINNED memory
		int num_streams = number_of_streams;//number of streams desired to create ::Note check for odd numbers
		unsigned long Each_LP_size = M * N;	//  sizeof(float);
		unsigned int num_LPs_perStream;
		bool equal_stream = true;
		if (N_S % num_streams == 0) {
			num_LPs_perStream = (N_S / num_streams);
			equal_stream = true;
		}
		else {
			num_LPs_perStream = (N_S / num_streams);//last stream will not be of the same size
			num_streams = num_streams + 1;//one extra streams.where nos of LPs to be solved will be less;
			equal_stream = false;
		}
		unsigned long lastBlock_size = 0;
		if (equal_stream == false) {
			lastBlock_size = N_S - (N_S / (num_streams - 1)) * (num_streams - 1);//LAST Stream Size
			//std::cout << "\nAmit Last Block size (LPs is )= " << lastBlock_size<< std::endl;
		}
		std::cout << "\nNumber of LPs_perStream = " << num_LPs_perStream << "  Total Streams=" << num_streams << " Last Block size=" << lastBlock_size << std::endl;

		// WINDOWS FIX allocate streams
		cudaStream_t *stream;
		stream = (cudaStream_t*)malloc(sizeof(cudaStream_t)*num_streams);

		cudaError_t result;
		// **************** Creation of Streams ****************
		for (int i = 0; i < num_streams; i++) {
			result = cudaStreamCreate(&stream[i]);
			unsigned long long offset = i * Each_LP_size * num_LPs_perStream;//for memory copy the starting address
			unsigned long offset_res = i * num_LPs_perStream;	//for result here offset_res is a pointer to the LP number
			if (equal_stream == false && i == (num_streams - 1)) {//last stream
				//int offset = i * Each_LP_size * lastBlock_size;	//for memory copy //Todo:: I think this offset computation is WRONG
				cudaMemcpyAsync(&G_MAT[offset], &MAT[offset], (lastBlock_size * M * N * sizeof(double)), cudaMemcpyHostToDevice, stream[i]);
				//printf("CUDA memcopyAsync G_MAT last-stream: %s\n", cudaGetErrorString(err));
				mykernel << <lastBlock_size, threads_per_block, 0, stream[i] >> >(G_MAT, M, N, G_R, lastBlock_size, R_data, R_index, offset_res);
				cudaMemcpyAsync(&R[offset_res], &G_R[offset_res], (lastBlock_size * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
			}
			else {
				//int offset = i * Each_LP_size * num_LPs_perStream;//for memory copy the starting address
				cudaMemcpyAsync(&G_MAT[offset], &MAT[offset], (num_LPs_perStream * M * N * sizeof(double)), cudaMemcpyHostToDevice, stream[i]);
				//printf("CUDA memcopyAsync G_MAT: %s\n", cudaGetErrorString(err));
				mykernel << <num_LPs_perStream, threads_per_block, 0, stream[i] >> >(G_MAT, M, N, G_R, num_LPs_perStream, R_data, R_index, offset_res);
				cudaMemcpyAsync(&R[offset_res], &G_R[offset_res], (num_LPs_perStream * sizeof(double)), cudaMemcpyDeviceToHost, stream[i]);
			}
		}

		// WINDOWS FIX free dynamically allocated stream array
		free(stream);
	}

	//	std::cout<<"Number of Artificial Variables = "<< rem.size()<<"\n";
	cudaFree(G_MAT);
	//cudaFree(G_Sel);
	//cudaFree(G_R);
	//	cudaFree(R_index);	//Only to synchronize with the cudamemcpy
	//	cudaFreeHost(MAT);
	//cudaFree(R_data);	//Only to synchronize with the cudamemcpy
	cudaDeviceSynchronize();		//removed as hopping that cudaFree will handle it
	if (multipleBatchOn >= 1){
		//std::cout<<"multipleBatch="<<multipleBatchOn<<std::endl;
		//std::cout<<"Before Reset R[0]="<<R[0]<<std::endl;
		//allResult = (double *) malloc(N_S * sizeof(double));
		allResult.resize(C.size1());
#pragma omp parallel for
		for (int i = 0; i < C.size1(); i++) {
			allResult[i] = R[i];
		}
		cudaDeviceReset(); //Reset the GPU. This solves multiple batch running but has problem in retaining values in Streaming approach due to the use of cudaMallocHost
		//	std::cout<<"R[0]="<<R[0]<<std::endl;
	}
	//cudaDeviceReset();
}

