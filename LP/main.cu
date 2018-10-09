#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include "cuda_profiler_api.h"
#include <glm/glm.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "FileIO.h"
#include "Auxilary.h"

//Whether Objective function should should minimise distance to point, comment out for minimising linear function
//#define OBJECTIVE_DISTANCE

#define maxSpeed_ 100000000000 // large artificial circular constraint
#define RVO_EPSILON 0.00001f //something close to zero

//#define printInfo
//number of threads in a block. Must be a multiple of 32!
#define BlockDimSize 128

//enum for min max
enum optimisation { MINIMISE, MAXIMISE };


////////////////////////////////////
//Numerical functions

//Determinant of 2 2d vectors
__device__ float det(const glm::vec2 v1, const glm::vec2 v2)
{
	return (v1.x * v2.y) - (v1.y * v2.x);
}

//Vector doted with itself
__device__ float absSq(const glm::vec2 &vector)
{
	return glm::dot(vector, vector);
}

//Square of value a
__device__ float sqr(float a)
{
	return a * a;
}



////////////////////////////////////
//Atomic operation functions

//Custom, unoptimized float atomic max
__device__ static float atomicMax(float* address, float val)
{
	int* address_as_i = (int*)address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_i, assumed, __float_as_int(fmaxf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}

//Custom, unoptimized float atomic min
__device__ static float atomicMin(float* address, float val)
{
	int* address_as_i = (int*)address;
	int old = *address_as_i, assumed;
	do {
		assumed = old;
		old = atomicCAS(address_as_i, assumed, __float_as_int(fminf(val, __int_as_float(assumed))));
	} while (assumed != old);
	return __int_as_float(old);
}



////////////////////////////////////
//Block level operations

/** \brief block level sum reduction. Result written to t0
* \param input_data thread value to be reduced over
* \returns sum reduction written to thread 0, else returns 0 for other threads.
*
* not valid for block size greater than 1024 (32*32)
*/
__device__ int reduce(int input_data) {
	__shared__ int s_ballot_results[BlockDimSize >> 5]; //shared results of the ballots
	int int_ret = 0; //value to return, only non zero for (threadIdx.x = 0)

	s_ballot_results[threadIdx.x >> 5] = ballot(input_data);
	__syncthreads();
	int blockCompNum;
	if (threadIdx.x < 32) { //0th warp and only threads that are within range - not valid for block size greater than 1024 (32*32)
		if (threadIdx.x >= (BlockDimSize >> 5))
			blockCompNum = 0;
		else
			blockCompNum = __popc(s_ballot_results[threadIdx.x]);
		for (int offset = 16; offset>0; offset >>= 1)
			blockCompNum += __shfl_down(blockCompNum, offset);
	}
	if (threadIdx.x == 0)
		int_ret = blockCompNum;
	return int_ret;
}

/** \brief compresses the thread varaible input_data using warp shuffles
* \param input_data thread level
* \param compArr shared memory array to store the compressed indices
* \param comp_num the size of the compressed array
*/
__device__ int compress(int input_data, int *compArr) {

	const int tid = threadIdx.x;
	__shared__ int temp[BlockDimSize >> 5]; //stores warp scan results shared
	int int_ret; //value to return

	int temp1 = input_data;
	//scan within warp
	for (int d = 1; d<32; d <<= 1) {
		int temp2 = __shfl_up(temp1, d);
		if (tid % 32 >= d) temp1 += temp2;
	}
	if (tid % 32 == 31) temp[tid >> 5] = temp1;
	__syncthreads();
	//scan of warp sums
	if (threadIdx.x < 32) {
		int temp2 = 0.0f;
		if (tid < blockDim.x / 32)
			temp2 = temp[threadIdx.x];
		for (int d = 1; d<32; d <<= 1) {
			int temp3 = __shfl_up(temp2, d);
			if (tid % 32 >= d) temp2 += temp3;
		}
		if (tid < blockDim.x / 32) temp[tid] = temp2;
	}
	__syncthreads();
	//add to previous warp sums
	if (tid >= 32) temp1 += temp[tid / 32 - 1];
	//compress
	if (input_data == 1) {
		compArr[temp1 - 1] = threadIdx.x;
	}

	//get total number - reduction
	int_ret = reduce(input_data);

	return int_ret;
}





////////////////////////////////////
//Linear program auxillaries

__device__ bool linearProgram1Fractions(const glm::vec4 lines_lineNo, const glm::vec4 lines_i, const float2 t2, float* tnew, bool *tLeftb)
{
	const glm::vec2 lines_direction_lineNo = glm::vec2(lines_lineNo.x, lines_lineNo.y);
	const glm::vec2 lines_point_lineNo = glm::vec2(lines_lineNo.z, lines_lineNo.w);

	const glm::vec2 lines_direction_i = glm::vec2(lines_i.x, lines_i.y);
	const glm::vec2 lines_point_i = glm::vec2(lines_i.z, lines_i.w);

	const float denominator = det(lines_direction_lineNo, lines_direction_i);
	const float numerator = det(lines_direction_i, lines_point_lineNo - lines_point_i);

	if (fabsf(denominator) <= RVO_EPSILON) {
		// Lines lineNo and i are (almost) parallel.
		if (numerator < 0.0f) {
			return false;
		}
		else {
			//continue, i.e. save a value that is guarateed to not affect the results
			*tnew = -INT_MAX; //an arbitary large value that is larger than tright
			*tLeftb = true;
			return true;
		}
	}

	const float t = numerator / denominator;

	if (denominator >= 0.0f) {
		// Line i bounds line lineNo on the right.
		*tLeftb = false;
	}
	else {
		// Line i bounds line lineNo on the left.
		*tLeftb = true;
	}
	*tnew = t;

	return true;
}

/*
* Solve constrainsts subject to a maximisation/minimisation function
* lines: glm::vec4 array of constraint lines of form {x gradient, y gradient, x point, y point }
* output: the 2 dimensional result of the linear program
* batches: number of batches to solve
* size: size of each batch
* objective_function: variables to minimise.
*/
__global__ void lpsolve(const glm::vec4 * const lines, glm::vec2 *output, const int batches, const int size, const glm::vec2* const objective_function) {
	//thread index
	const int index = (blockIdx.x * blockDim.x) + threadIdx.x;
	const int tid = threadIdx.x;

	//initialize SM
	__shared__ int compArr[BlockDimSize]; //shared compressed array working list
	__shared__ int active_agents[1]; //shared number of threads in block that are in the active compression
	__shared__ glm::vec4 s_line[BlockDimSize]; //shared current orca line of interest x:direction.x y:direction.y z:point.x w:point.y
	__shared__ float2 s_t[BlockDimSize]; //tleft and tright shared
	__shared__ int s_lineFail[BlockDimSize]; //on which neighbour number the lp has failed shared. -1 if succeeded
	__shared__ glm::vec2 s_newv[BlockDimSize]; //solution
	__shared__ glm::vec2 s_desv[BlockDimSize]; //position to be closest to/value that optimises objective function
	//whether to minimise or maximise objective function
	enum optimisation optimiseFunc = MAXIMISE;

	//Initialise variables
	s_desv[tid] = objective_function[index]; //initialise position to be closest to / arguments of the objective function
	s_newv[tid] = (float)INT_MAX * glm::normalize(objective_function[index]) * ((optimiseFunc == MAXIMISE) ? 1.0f : -1.0f); //starting value guess. Choosen as to maximise/minimise objective function
	s_lineFail[tid] = -1; //initialise success integer
	s_t[tid] = make_float2(-INT_MAX, INT_MAX);

	//number of starting agents = all possible threads in block
	if (tid == 0) {
		(*active_agents) = (batches < (blockIdx.x * blockDim.x)) ? (batches % blockDim.x) : blockDim.x;
	}
    //sync sm data
	//__syncthreads();


	//loop through all lines in the batch and set per loop variables
	for (int i = 0; i < size; i++, s_t[tid] = make_float2(-INT_MAX, INT_MAX)) {

		//Whether this thread is performing computation (1). Needs to be integer due to block level operations, and check early out condition
		int bthread_data = (index < batches && s_lineFail[tid] == -1) ? 1 : 0;

		if (bthread_data == 1) {
			//Write line info for this iteration into SM
			s_line[tid] = lines[(index*size) + i];

			//check if newVel is satisfied by the constraint line. 1 if not satisfied and requires work, otherwise 0.
			bthread_data = (int)(det(glm::vec2(s_line[tid].x, s_line[tid].y), glm::vec2(s_line[tid].z, s_line[tid].w) - s_newv[tid]) > 0.0f); //issues on this line
		}


		//compress through exclusive scan
		int result = compress(bthread_data, compArr);
		//result written to thread0
		if (tid == 0) {
			(*active_agents) = result;
			//if (result == 0)
			//	break;
		}

		//For compArr to be filled properly
		__syncthreads();

#ifdef printInfo
		//index to examine
#define ite 0
		if (index == ite) {
			printf("---------------------------\niteration %d\n", i);
		}
		//if (active_agents != 0 && tid == 0){
		//	printf("Block %i active threads %d\n", blockIdx.x, result);
		//}
		if (ite == index) {
			printf("Current optimal value = (%5.3f, %5.3f)\n", s_newv[threadIdx.x].x, s_newv[threadIdx.x].y);
			printf("Line being considered = (%5.3f, %5.3f)dir (%5.3f, %5.3f)point\n", s_line[threadIdx.x].x, s_line[threadIdx.x].y, s_line[threadIdx.x].z, s_line[threadIdx.x].w );
			if (bthread_data) {
				printf("Performing recomputation for this line\n");
			}
			else {
				printf("No computation for this line needed\n");
			}
		}
#endif

		//calculate the total number of work unit items (where a work unit is a line read and calculation for a unqiue agent line index). i.e.
		int wu_count = ((*active_agents) * i);

		//divide work unit items between threads. i.e.
		for (int j = 0; j < wu_count; j += blockDim.x) {

			//calculate unique work unit index
			int wu_index = j + tid;

			//do work if there are still wu to complete
			if (wu_index < wu_count) {

				//for each thread work out which agent it is associated with
				int n_tid = compArr[wu_index / i];

				//for each thread work out which line index it should read
				int line_index = wu_index % i;

				//read in the unique agent line combination using the calculated indices
				int newIndex = n_tid + blockIdx.x*blockDim.x;
				glm::vec4 lines_i = lines[(newIndex*size) + line_index];

#ifdef printInfo
				if (n_tid + blockDim.x * blockIdx.x == ite) {
					printf("sub line %i considered = (%5.3f, %5.3f)dir (%5.3f, %5.3f)point\n", line_index, lines_i.x, lines_i.y, lines_i.z, lines_i.w );
				}
#endif

				//calculate denominator and numerator
				bool btleft;//whether the t value is left (or right if false)
				float t;//value of t
				if (!linearProgram1Fractions(s_line[n_tid], lines_i, s_t[n_tid], &t, &btleft)) {
					//operation failed
					s_lineFail[n_tid] = i;
#ifdef printInfo
					if (n_tid + blockDim.x * blockIdx.x == ite) {
						printf("Computation failed, half-planes are parallel in wrong direction\n");
					}
#endif
				}
#ifdef printInfo
				if (n_tid + blockDim.x * blockIdx.x == ite) {
					printf("sub line %i writing to bleft? %i, value = %5.3f\n", line_index, btleft, t);
				}
#endif
				//atomic write tleft and tright to shared memory using an atomic min and max
				if (btleft) {
					atomicMax(&s_t[n_tid].x, t);
				}
				else {
					atomicMin(&s_t[n_tid].y, t);
				}
			}
		}

		//sync to ensure all atomic writes are complete
		__syncthreads();

		//update the new velocity for each active agent
		if (tid < (*active_agents)) {

			//New index
			int n_tid = compArr[tid];

			//failure condition if no region of validity
			if (s_t[n_tid].x > s_t[n_tid].y) {
				s_lineFail[n_tid] = i;
#ifdef printInfo
				if (n_tid + blockDim.x * blockIdx.x == ite) {
					printf("Computation failed, no valid region\n");
				}
#endif
			}

			//If not failed up to this point
			if (s_lineFail[n_tid] == -1) {

				// Optimize closest point
				glm::vec2 lineDir = glm::vec2(s_line[n_tid].x, s_line[n_tid].y);
				glm::vec2 linePoint = glm::vec2(s_line[n_tid].z, s_line[n_tid].w);
				//Change this line to alter the functional form of the optimisation function
#ifdef OBJECTIVE_DISTANCE
				//for case of minimising distance to point @s_desv
				const float t = glm::dot(lineDir, s_desv[n_tid] - linePoint);

				//best value is to the left of what is allowed
				if (t < s_t[n_tid].x) {
					s_newv[n_tid] = linePoint + s_t[n_tid].x * lineDir;
				}
				//best value is to the right of what is allowed
				else if (t > s_t[n_tid].y) {
					s_newv[n_tid] = linePoint + s_t[n_tid].y * lineDir;
				}
				//best value is not on a vertex
				else {
					s_newv[n_tid] = linePoint + t * lineDir;
				}
#else
				//for case of minimising linear function
				//solution is guaranteed to be either on vertex of t_left or right, or along the line joining them (in which case the precise choise does not matter)

				//the objective function
				glm::vec2 fct = s_desv[n_tid];


				glm::vec2 t_left_sln = linePoint + s_t[n_tid].x * glm::vec2(lineDir.x, lineDir.y); //x,y value at t_left
				float t_left_val = glm::dot(fct, t_left_sln); //value of objective function at t_left
				//value of function at t_right
				glm::vec2 t_right_sln = linePoint + s_t[n_tid].y * glm::vec2(lineDir.x, lineDir.y); //x,y value at t_right
				float t_right_val = glm::dot(fct, t_right_sln); //value of objective function at t_right
				//assign answer from correct t.
				s_newv[n_tid] = ((t_left_val > t_right_val) != (optimiseFunc == (MINIMISE))) ? t_left_sln : t_right_sln;

#endif // OBJECTIVE_DISTANCE

#ifdef printInfo
				if (n_tid + blockDim.x * blockIdx.x == ite) {
					printf("New optimal value = (%5.3f, %5.3f)\n", s_newv[threadIdx.x].x, s_newv[threadIdx.x].y);
				}
#endif
			}
		}


#ifdef printInfo
#ifdef tte
		if (index == tte && thread_data == 1) {
			printf("Current iter %d soltn x: %5.3f y: %5.3f \n\n", i, s_newv[index].x, s_newv[index].y);
		}
#endif
#endif

		//sync to ensure all shared mem writes are complete
		__syncthreads();

	}

	//write to output
	if (index < batches) {
		output[index] = s_newv[tid];
	}

#ifdef printInfo
	//print results
	if (tid == 0) {
		for (int i = 0; i < size; i++) {
			printf("------------------------------\nFinal index %d soltn x: %5.3f y: %5.3f \n", i+(blockIdx.x*blockDim.x), s_newv[i].x, s_newv[i].y);
		}
	}
#endif

}



////////////////////////////////////
//Main

int main(int argc, const char* argv[])
{
	int batches = 0; //number of LPs
	int size = 0; //size of each LP
	glm::vec2* x = NULL;  //< solution to randomly generated LP
	glm::vec4* constraintsSingle = NULL; //array of constraints for 1 lp
	glm::vec2 optimiseSingle; // variables to minimise in optimisation function for 1 lp
	glm::vec4* constraints = NULL; // array of constraints

	glm::vec2* optimise = NULL; //variable to minimise in optimisation function
	glm::vec2* output = NULL; // array of optimal results in x & y.
	//float* outputVals = NULL; //array of optimal values, i.e. plug x & y into optimisation function


	//------------------------------------------
	//handle args
	if (argc != 3) {
		printf("\nIncorrect Number of Arguments!\n");
		printf("Correct Usages/Syntax:\n");
		printf("Argument 1) File-Name -- Name of input file in benchmarks folder. Cannot contain spaces\n");
		printf("Argument 2) Batch-size -- Number of LPs to be solved\n");
		return 1;
	}
	batches = atoi(argv[2]);

	//------------------------------------------
	//handle input

	//Input is from file
	printf("Parsing input files... ");
	if(!parseBenchmark4(argv[1], &constraintsSingle, &optimiseSingle, &size)){
		return 1;
	}
	printf("Done\n");




	//------------------------------------------
	//initialize memory initialsation timings
	cudaEvent_t start, stop;
	float memory_milliseconds = 0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	//start time
	cudaEventRecord(start);

	//------------------------------------------
	//memory allocation

	gpuErrchk(cudaMallocManaged(&output, sizeof(glm::vec2) * batches));
	gpuErrchk(cudaMallocManaged(&constraints, sizeof(glm::vec4) * batches * size));
	gpuErrchk(cudaMallocManaged(&optimise, sizeof(glm::vec2) * batches));



	//------------------------------------------
	//Tile LP data to all other LPs
	for (int i = 0; i < batches; i++) {
		for (int j = 0; j < size; j++) {
			constraints[i*size + j] = constraintsSingle[j]; //deep copy constraint data
		}
		optimise[i] = optimiseSingle;
	}

    //end time
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&memory_milliseconds, start, stop);
	//printf("Memory processing time: %f (ms)\n", memory_milliseconds);


	//------------------------------------------
	//Initialize kernel
	int blockSize = BlockDimSize;
	int gridSize = (batches + blockSize - 1) / blockSize;
	dim3 b, g;
	b.x = blockSize;
	g.x = gridSize;

    //------------------------------------------
	//initialize timings
	float milliseconds = 0;
	//start time
	cudaEventRecord(start);

	//------------------------------------------
	//Kernel execution


	//kernel
	lpsolve <<< g, b >>>(constraints, output, batches, size, optimise);

	//end time
	cudaEventRecord(stop);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&milliseconds, start, stop);
	//printf("Kernel Processing time: %f (ms)\n", milliseconds);

	cudaDeviceSynchronize();



	//------------------------------------------
	//print results to console
	int numResultsToPrint = (batches < 1) ? batches : 1;
	for (int i = 0; i < numResultsToPrint; i++) {
		printf("Batch %i \t Optimal location is x: %f y: %f \t value of %f\n", i, output[i].x, output[i].y, glm::dot(optimise[i], output[i]) );
	}

    //------------------------------------------
	//write timing to file
	writeTimingtoFile("timings/RGBMemtimings.txt", size, batches, memory_milliseconds );
	writeTimingtoFile("timings/RGBCalctimings.txt", size, batches, milliseconds);
	writeTimingtoFile("timings/RGBtimings.txt", size, batches, memory_milliseconds+milliseconds);

	//------------------------------------------
	//cleanup

	//memory
	gpuErrchk(cudaFree(constraints));
	free(x);


	//cuda device reset
	gpuErrchk(cudaDeviceReset());

  return 0;
}
