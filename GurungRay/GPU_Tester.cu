/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

/*
 * This Project contains a comparison of GLPK solver running multiple LPs sequentially one after another with our batched LP solver in GPU.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "simplex.cuh"
#include "math/glpk_lp_solver/glpk_lp_solver.h"
#include <vector>
#include "math/matrix.h"
#include <climits>
#include <iostream>
 //whether to use boost or chrono, chrono is consistent across tested CPU methods
#define CHRONO
#ifdef CHRONO
#include <chrono>
#else
#include "boost/timer/timer.hpp"
#endif
#include <fstream>
#ifdef WIN32
	
#else
	#include "sys/time.h"
#endif
#include "parseBenchmark.h"
#include "dataStructure.h"
#include <omp.h>


int main(int argc, char *argv[]) {
	unsigned int LP_size = 1, avg;
	unsigned int MaxMinFlag=1; //1 for Min and 2 for Max
	math::matrix<double> A;
	math::matrix<double> C, newC;
	std::vector<double> b, c;
	std::vector<double> result;
	std::vector<int> status_val;
	std::string afile, bfile, cfile, efile;



	if (argc > 1) {
		if (argc != 6) { //1(ApplicationName) + 5 (Input Arguments)
		//if (argc != 5) { //1(ApplicationName) + 4 (Input Arguments)
			std::cout << "\nInsufficient Number of Arguments!!!\n";
			std::cout << "Correct Usages/Syntax:\n";
			std::cout << "./ProjName   'Benchmark-NO'   'Average'  'Batch-Size'\n";
			std::cout << "Argument 1) Benchmark-Name -- The Benchmark Name to use, e.g. 'benchmark/256'\n";
			std::cout << "Argument 2) Average -- Number of Average reading to be taken\n";
			std::cout << "Argument 3) Batch-size -- Number of LPs to be solved \n";
			std::cout << "Argument 4) Stream-size -- Number of Streams to be used \n";
			std::cout << "Argument 5) 1/2 -- 1 for Benchmarks and 2 for dense-Random LP Experiment \n";
			/*
			 * Note if the last Argument is selected as 1 then the first Argument values determines the Benchmarks to be executed
			 * but if the last  Argument selected is 2 then the value entered in the first Argument is the dimension of the dense-Random LP
			 */
			return 0;
		}
		int randOrBenchmark;
		randOrBenchmark=atoi(argv[5]);

		if (randOrBenchmark==1){
			//selectBenchmark(atoi(argv[1]), afile, bfile, cfile);//Selecting the required benchmark files having Matrix A, vector b and c.
			selectBenchmarkNamed(argv[1], afile, bfile, cfile);//Selecting the required benchmark files having Matrix A, vector b and c.
			parseLP(afile.c_str(), bfile.c_str(), cfile.c_str(), A, b, c, MaxMinFlag);//parse the selected Benchmark files to convert into Matrix A along with vector b and c.
		}else{
			generateRandomLP(atoi(argv[1]),A,b,c,MaxMinFlag);
			/*std::cout << "\nA Matrix is \n";
			std::cout << A << std::endl;
			std::cout << "\nb Vector is \n";
			for (int i = 0; i < b.size(); i++)
				std::cout << b[i] << "\t";
			std::cout << "\nc Vector is \n";
			for (int i = 0; i < c.size(); i++)
				std::cout << c[i] << "\t";*/

			return 0;
		}
		// ************************ Smaller Test LPs ******************************************
		//parseSmallLP2(A,b,c, MaxMinFlag); //has solution 39.44444
		//parseSmallLP3(A,b,c, MaxMinFlag); //has solution 21.33333
		//parseSmallLP4(A,b,c, MaxMinFlag); //has solution 9.33333
		//parse5Dim(A,b,c, MaxMinFlag); //for solution see inside the function the value of c
		//parse28DimHelicopter(A,b,c, MaxMinFlag); //
		//parseSmallLP(A,b,c, MaxMinFlag); //Unbounded problem has many solutions at 2
		// ************************ Smaller Test LPs ******************************************



		avg = atoi(argv[2]);
		LP_size = atoi(argv[3]);
		unsigned int numberOfStreams=atoi(argv[4]);

		/*
		 * Experimental Note:
		 *  To record the average time we ignore the first reading taken in GPU as it does not
		 *  reflect the correct computation time because for the first GPU call it also include an
		 *  extra overhead of GPU initialization time.
		 *
		 *  Note: The NetLib benchmarks are taken as
		 *  		Minimize cx
		 *  		sub to   Ax <= b
		 *  		For All  x>=0
		 *
		 *  so, we convert it as cx to -cx and display the result as -1 * (the computed result)
		 *  		-1 * (Maximize -cx)
		 *  		sub to   Ax <= b
		 *  		For All  x>=0
		 *
		 */

#ifdef CHRONO
		std::chrono::duration<double, std::milli> sum;
#else
		double sum = 0.0;
		double wall_clock, return_Time;
		boost::timer::cpu_timer tt1, tt2, tt3;	//tt1 -- Variable declaration
#endif
		//Since our GPU LP Solver consider Maximize by default so we multiply -1 to the objective function

		C.resize(LP_size, c.size());	//objective function
		for (unsigned int i = 0; i < C.size1(); i++) {
			for (unsigned int j = 0; j < C.size2(); j++) {
				if (MaxMinFlag == 1)
					C(i, j) = -1 * c[j]; //Given problem is Minimize so converting to work on Maximization
				else
					C(i, j) = c[j]; //Given problem is Maximize
			}
		}

		/*for (unsigned int i = 0; i < C.size1(); i++) {
			for (unsigned int j = 0; j < C.size2(); j++) {
				std::cout<<C(i, j)<<"\t"; //Given problem is Maximize
			}
			std::cout<<"\n";
		}*/


		std::cout << std::fixed; //to assign precision on the std::output stream
		std::cout.precision(17); //cout << setprecision(17);

		//Computation for CPU ie GLPK
		//	sum=0.0;
		//std::cout << "\n*****GLPK RESULT*****\n";
		std::vector<double> dir(c.size());

		//***** MODEL SELECTION *****
		double res = 0.0;
		double batchTime = 0.0, AvgBatchTime = 0.0;
		std::vector<double> resul(C.size1());

/*
 * Before Solving the given batch-size of LPs
 * 1) Let me first find out the memory size requirement of one (1) LP
 * 2) Then find if possible to solve in a single batch.
 * 3) If not possible in one batch then call it multiple times.
 *
 */
// ************ Calculating size of One LP ******************************
		unsigned int rowSize= A.size1(), colSize=A.size2();
		unsigned long long oneLPsize;
		std::vector<int> rem;
		for (int i = 0; i < b.size(); i++) {
			if (b[i] < 0) {
				rem.push_back(i);
			}
		}
		colSize =colSize + rowSize + rem.size() + 3; //tableau column size
		oneLPsize = (rowSize+2) * colSize * sizeof(double); //size in Bytes

		/*
		 * kernel(double *S_MAT, int S_row, int S_col, double *Result, int S_N, double *R_data, int *R_index)
		 * To compute precise memory requirement (excluding shared memory and registers)
		 *
		 */
		unsigned int LPVarySize, fixForAll;
		fixForAll = sizeof(int) + sizeof(int) + sizeof(int); //4 bytes each for S_row, S_col and S_N:: This is fix for any number of LPs
		LPVarySize = sizeof(double);// * LP_size; //LP_size number of result of double DataType is returned
		unsigned int R_data_size = sizeof(double) * colSize;// * LP_size;
		unsigned int R_index_size = sizeof(int) * colSize;// * LP_size;
		LPVarySize = LPVarySize + oneLPsize + R_data_size + R_index_size;
		oneLPsize = fixForAll + LPVarySize; //This size of oneLPsize also include fixForAll size which take into account for any number of LPs. Thus (totalGMemSize / oneLPsize) is overApproximation
//		std::cout<<"oneLPsize = "<<oneLPsize<<std::endl;
		oneLPsize = oneLPsize + 1024 * 10;//testing with extra 1 KB per LP
//		std::cout<<"oneLPsize = "<<oneLPsize<<std::endl;
// ************ Size of One LP ******************************
		unsigned long long totalGMemSize;
		int device;
		cudaDeviceProp props;
		cudaGetDevice(&device);
		cudaGetDeviceProperties(&props, device);
		totalGMemSize = props.totalGlobalMem;
//		std::cout<<"totalGMemSize = "<<totalGMemSize<<std::endl;
		//totalGMemSize = props.totalGlobalMem / 1024; //converting into KiloBytes
// ************ Size of One LP ******************************

		// ***************  How many LPs can be solved? ********************************
		unsigned int numberSolved = totalGMemSize / oneLPsize; // returns Maximum LPs that can be solved in one batch
//		std::cout<<"numberSolved = "<<numberSolved<<std::endl;

		int multipleBatch=0;//Multiple batch NOT required

		if (LP_size > numberSolved)
			multipleBatch=1;//Require multiple batch solving
		// ***************  How many LPs can be solved? ********************************
int max=10;//Decide to view max number of LPs as verification for correctness of GPU compared to GLPK

		if (multipleBatch ==0){ //Does NOT require multiple batch solving
			for (unsigned int i = 0; i <= avg; i++) {
				Simplex s(C.size1(), multipleBatch);
				s.setConstratint(A, b);	//setting constraints also recorded like in GLPK for independent LP
#ifdef CHRONO
				auto tt1 = std::chrono::high_resolution_clock::now();
#else
				tt1.start();
#endif
				s.ComputeLP(C, numberOfStreams);
#ifdef CHRONO
				auto tt2 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> return_Time = tt2 - tt1;
				sum += return_Time;

#else
				tt1.stop();
				if (i == 0) {
					wall_clock = tt1.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
					return_Time = wall_clock / (double)1000; //convert milliseconds to seconds
					std::cout << "Iter 0: " << "Time = " << return_Time << " Seconds" << std::endl;
				}
				else {
					wall_clock = tt1.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
					return_Time = wall_clock / (double)1000; //convert milliseconds to seconds
					std::cout << "Iter " << i << ": Time = " << return_Time << " Seconds" << std::endl;
					sum = sum + return_Time; //convert nanoseconds to milliseconds
				}
#endif
				result = s.getResultAll();
			//	std::cout<<"results[0]="<<result[0]<< "  = "<<s.getResultAll()[1]<<std::endl;
			}

#ifdef CHRONO
			std::chrono::duration<double, std::milli> Final_Time2 = sum / (double)avg;
			std::cout << "\n*****GPU RESULT*****\n";
			std::cout << "\nCHRONO Time taken:Wall  (in ms):: GPU= " << Final_Time2.count() << std::endl;

#else
			double	Final_Time2 = sum / avg;
			std::cout << "\n*****GPU RESULT*****\n";
			std::cout << "\nBoost Time taken:Wall  (in Seconds):: GPU= " << (double) Final_Time2 << std::endl;

			//std::cout << "\nNumber of Simplex Solved = " << C.size1() << std::endl;
#endif
			std::cout << "\n**Answer_Of_All_Simplex**\n";
			//int max = 5;	//LP_size;	//Verifying results of only first 5 LPs.
			if (LP_size < max)
				max = LP_size;

			std::cout << "\nVERIFICATION FOR CORRECTNESS\n";
			for (int i = 0; i < max; i++) {
				if (MaxMinFlag == 1)
					std::cout << "GLPK: " << -1 * resul[i] << " || GPU: " << -1 * result[i] << " || GPU: " << -1 * result[LP_size - 1 - i] << std::endl;
				else
					std::cout << "GLPK: " << resul[i] << " || GPU: " << result[i] << " || GPU: " << result[LP_size - 1 - i] << std::endl;
			}


		} else {//Requires multiple batch solving
			/*
			 * Find out how many batches need to be solved?
			 * To be solved = LP_size
			 * Can be solved = numberSolved
			 */

			unsigned int timesSolved = (LP_size / numberSolved);
			if ((LP_size % numberSolved)==0)
				;//exact and equal LPs per batch times
			else
				timesSolved = timesSolved + 1; //unequal LPs per batch times OR last batch has less LPs to be solved
			std::cout << "Total number of Batches " << timesSolved << std::endl;
			// **************** Breaking into batches ***************************
			std::list<block_lp> bulk_lps;// size is created on every push_back// (timesSolved);	//list of sub-division of LPs
			struct block_lp myLPList;
			myLPList.block_obj_coeff.resize(numberSolved, C.size2());
			math::matrix<double> block_obj_coeff(numberSolved, C.size2());
			unsigned long long index = 0;
			for (unsigned long long lp_number = 0; lp_number < LP_size; lp_number++) {
				for (unsigned long long i = 0; i < C.size2(); i++) {
					myLPList.block_obj_coeff(index, i) = C(lp_number, i);
				}
				index++;
				if (index == numberSolved) {
					index = 0;
					bulk_lps.push_back(myLPList);
				}
			}	//end of all LPs

			if (bulk_lps.size() != timesSolved){
				myLPList.block_obj_coeff.resize(index, C.size2(),true);
				//std::cout<<"myLPList.block_obj_coeff\n"<<myLPList.block_obj_coeff;
				bulk_lps.push_back(myLPList);//the last created data
				//std::cout<<"\myLPList.block_obj_coeff.size1()=" <<myLPList.block_obj_coeff.size1()<<std::endl;
			}
//			std::cout<<"\nbulk_lps.size()=" <<bulk_lps.size()<<std::endl;
//			std::list<block_lp>::iterator it = bulk_lps.begin();
//			std::cout<<"\nbulk_lps.size1()=" <<(*it).block_obj_coeff.size1()<<std::endl;
//			std::cout<<"\nbulk_lps.size2()=" <<(*it).block_obj_coeff.size2()<<std::endl;

			// **************** Breaking into batches ***************************
			// **************** Solving each batches ***************************
			std::list<block_lp_result> bulk_result;// size is created on every push_back// (timesSolved);
			//struct block_lp_result eachBlock;
			//eachBlock.results.resize(numberSolved);	//last block will be less

			//sum = 0.0;
			for (std::list<block_lp>::iterator it = bulk_lps.begin(); it != bulk_lps.end(); it++) {
				struct block_lp_result eachBlock;
				math::matrix<double> obj_coeff;
				//std::cout<<"(*it).block_obj_coeff.size1()="<<(*it).block_obj_coeff.size1()<<std::endl;
				Simplex s((*it).block_obj_coeff.size1(), multipleBatch);
				s.setConstratint(A, b);	//setting constraints also recorded like in GLPK for independent LP
				obj_coeff = (*it).block_obj_coeff;
				//std::cout<<"obj_coeff\n"<<obj_coeff;
#ifdef CHRONO
				auto tt3 = std::chrono::high_resolution_clock::now();
#else
				tt3.start();
#endif
				s.ComputeLP(obj_coeff, numberOfStreams);
#ifdef CHRONO
				auto tt4 = std::chrono::high_resolution_clock::now();
				std::chrono::duration<double, std::milli> return_Time = tt4 - tt3;
				sum += return_Time;
				std::cout << "\nBatch Time (in Seconds):: GPU= " << return_Time.count() << std::endl;
#else
				tt3.stop();

				wall_clock = tt3.elapsed().wall / 1000; //convert nanoseconds to milliseconds
				return_Time = wall_clock / (double) 1000; //convert milliseconds to seconds
				sum = sum + return_Time; //convert nanoseconds to milliseconds
				std::cout << "\nBatch Time (in Seconds):: GPU= " << (double) return_Time<< std::endl;
#endif
				eachBlock.results = s.getResultAll();
				bulk_result.push_back(eachBlock);

				multipleBatch++;
			//	std::cout<<"eachBlock.results[0]="<<eachBlock.results[0]<< "  = "<<s.getResultAll()[1]<<std::endl;
			//	std::cout<<"eachBlock.results.size="<<eachBlock.results.size()<<std::endl;
			}
			// **************** Solving each batches ***************************
			// **************** Obtaining Result of each batches (if Required) ***************************
			std::vector<double> res(LP_size); //Final Result is here
			unsigned long long index_res = 0;
			for (std::list<block_lp_result>::iterator it = bulk_result.begin();
					it != bulk_result.end(); it++) {
				unsigned long block_result_size = (*it).results.size();
				for (unsigned long i = 0; i < block_result_size; i++) {
					res[index_res] = (*it).results[i];
					index_res++;
				}
			}
			// **************** Obtaining Result of each batches ***************************
			std::cout << "\n*****GPU RESULT*****\n";
#ifdef CHRONO
			std::cout << "\nChrono Time taken:Wall  (in Seconds):: GPU= " << sum.count() << std::endl;
#else
			std::cout << "\nBoost Time taken:Wall  (in Seconds):: GPU= " << (double) sum << std::endl;
#endif
			std::cout << "\n**Answer_Of_All_Simplex**\n";
			//int max = 5;	//LP_size;	//Verifying results of only first 5 LPs.
			if (LP_size < max)
				max = LP_size;

			std::cout << "\nVERIFICATION FOR CORRECTNESS\n";
			for (int i = 0; i < max; i++) {
				if (MaxMinFlag == 1)
					//std::cout << "GLPK: " << -1 * resul[i] << " || GPU: " << -1 * res[i] << std::endl;
					std::cout << "GLPK: " << -1 * resul[i] << " || GPU: " << -1 * res[i]<< std::endl;
				else
					std::cout << "GLPK: " << resul[i] << " || GPU: " << res[i] << std::endl;
			}
		}

		std::cout << "Writing timing results to files...";
		//write timings to file
		std::ofstream timingFile;
		int size = A.size1();
		int batches = LP_size;


		//GPU
		timingFile.open("timings/GRtimings.txt", std::ios::app);
#ifdef CHRONO
		timingFile << size << "\t" << batches << "\t" << sum.count() << std::endl;
#else
		timingFile << size << "\t" << batches << "\t" << sum << std::endl;
#endif

		timingFile.close();
		std::cout << "Done!" << std::endl;
	}
	return 0;

}


