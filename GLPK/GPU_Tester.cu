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
//#include "simplex.cuh"
//#include "math/glpk_lp_solver/glpk_lp_solver.h"
#include <vector>
//#include "math/matrix.h"
#include <climits>
#include <iostream>
#include "boost/timer/timer.hpp"
#include <fstream>
#ifdef WIN32
	
#else
	#include "sys/time.h"
#endif
#include "parseBenchmark.h"
//#include "dataStructure.h"
#include <omp.h>


int main(int argc, char *argv[]) {
	unsigned int LP_size = 1, avg;
	unsigned int MaxMinFlag=1; //1 for Min and 2 for Max
	math::matrix<double> A;
	math::matrix<double> C, newC;
	std::vector<double> b, c;
	std::vector<double> result;
	std::vector<int> status_val;
	double Final_Time;
	std::string afile, bfile, cfile, efile;



	if (argc > 1) {
		if (argc != 6) { //1(ApplicationName) + 5 (Input Arguments)
		//if (argc != 5) { //1(ApplicationName) + 4 (Input Arguments)
			std::cout << "\nInsufficient Number of Arguments!!!\n";
			std::cout << "Correct Usages/Syntax:\n";
			std::cout << "./ProjName   'Benchmark-NO'   'Average'  'Batch-Size'\n";
			std::cout << "Argument 1) Benchmark-NO -- The Benchmark No or Random LP dimension\n";
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
			selectBenchmark(atoi(argv[1]), afile, bfile, cfile);//Selecting the required benchmark files having Matrix A, vector b and c.
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

		double sum = 0.0;
		double wall_clock, return_Time;
		boost::timer::cpu_timer tt1, tt2, tt3;	//tt1 -- Variable declaration

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


		for (int k = 1; k <= avg; k++) {
			//batchTime = 0.0;
			tt1.start();
#pragma omp parallel for num_threads(4)
			for (int i = 0; i < C.size1(); i++) {
				glpk_lp_solver mylp;
				mylp.setMin_Or_Max(2);	//2 for Maximize and 1 for Minimize
				for (int j = 0; j < c.size(); j++) {
					dir[j] = C(i, j);
				}
				mylp.setConstraints(A, b, 1); //this function actually determines independent LP in GLPK
				res = mylp.Compute_LLP(dir); //We consider every dir an independent LP problem
				//res = mylp.Compute_LLP(c); //We consider every dir an independent LP problem
				resul[i] = res;
			}
			tt1.stop();
			wall_clock = tt1.elapsed().wall / 1000000; //convert nanoseconds to milliseconds
			return_Time = wall_clock / (double) 1000; //convert milliseconds to seconds
			batchTime = return_Time; //convert nanoseconds to milliseconds
			AvgBatchTime = AvgBatchTime + batchTime;
		}


		Final_Time = AvgBatchTime / avg;

		std::cout << "\n*****GLPK RESULT*****\n";
		std::cout << "\nBoost Time taken:Wall  (in Seconds):: GLPK= " << (double) Final_Time << std::endl;

		/*std::cout << "\nVERIFICATION FOR CORRECTNESS\n";
		for (int i = 0; i < LP_size; i++) {
			if (MaxMinFlag == 1)
				std::cout << "GLPK: " << -1 * resul[i]<< std::endl;
			else
				std::cout << "GLPK: " << resul[i] << std::endl;
		}*/



		std::cout << "Writing timing results to files...";
		//write timings to file
		std::ofstream timingFile;
		int size = A.size1();
		int batches = LP_size;
		//CPU
		double millisecondsCPU = Final_Time*1000.0;
		timingFile.open("BLPG_CPU_OMP_timing.txt", std::ios::app);
		timingFile << size << "\t" << batches << "\t" << millisecondsCPU << std::endl;
		timingFile.close();

	}
	return 0;

}


