/*
 * This Project contains a comparison of GLPK solver running multiple LPs sequentially one after another with our batched LP solver in GPU.
 *
 */
 //whether to use boost or chrono, chrono is consistent across tested CPU methods
#define CHRONO
 //whether to use multiple core (OMP), or single core
#define OMP

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <climits>
#include <iostream>
#include <fstream>
#ifdef CHRONO
#include <chrono>
#else
#include "boost/timer/timer.hpp"
#endif
#ifdef WIN32
	
#else
	#include "sys/time.h"
#endif
#include "parseBenchmark.h"

#ifdef OMP
#include <omp.h>
#endif


void generateRandomLP2D(int dimension, math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag) {
	unsigned int status = 0, res1 = 0, res2 = 0;
	// *************************************************************************************************************************************
	unsigned int art = 0; //Decides the number of artificial variables (or make the initial basic solution of LP as Feasible or Infeasible)
	// *************************************************************************************************************************************
	time_t t;
	/* Intializes random number generator */
	srand((unsigned)time(&t));

	int sign = -1, sum1;
	glpk_lp_solver lp;

	while (status != 5) {
		//We are assuming that LP is in the form A<=b
		A.resize(dimension, 2);
		b.resize(dimension);

		for (unsigned int j = 0; j < dimension; j++) {
			for (unsigned int k = 0; k < 2; k++) {
				//A(j, k) = rand() % (k + 10) + 1;
				A(j, k) = (rand() % 999) + 1; // Range 1 to 1000
			}

			if (j < art) {//for LPs with initial basic solution is Infeasible
				res1 = (rand() % 999);
				res2 = (res1 + 1);
				sum1 = res2 * sign;
				b[j] = (double)sum1;
				A(j, j) = A(j, j) * sign; //This is required to make the LP feasible -A(j,j)<=-b[j]
			}
			else	//for LPs with initial basic solution is  Feasible
					//b[j] = (rand() % (j + 1) + (10 + j));
				b[j] = (rand() % 999) + 1;  // Range 1 to 1000
		}

		//** Setting current Lp to GLPK
		lp.setMin_Or_Max(MaxMinFlag);
		lp.setConstraints(A, b, 1);
		status = lp.TestConstraints(); //The objective function is all set ot zero. //std::cout << "Status = " << status << "\n";
	}

	/*
	* Experimental Note:
	*  We generate randomly the first feasible LP problem with
	*  Matrix A (in the range 1 to 1000)
	*  Vector b (in the range 1 to 1000) and
	*  The objective function c (in the range 1 to 500).
	*  We use the system timestamp as the seed for the random generator
	*/

	c.resize(2);
	unsigned int i;
	for (unsigned int j = 0; j < c.size(); j++) {
		//c[j] = rand() % (j + 1) + 1;
		c[j] = (rand() % 499) + 1; //Range 1 to 500
	}
	double res = lp.Compute_LLP(c);

	std::cout << "GLPK::Solution = " << res << std::endl;
}


void WriteToFile(int dimension, math::matrix<double> A, std::vector<double> b, std::vector<double> c) {
	std::string fileToWrite = "benchmarks/" + std::to_string(dimension);
	std::string fileA = fileToWrite + "_A.txt";
	std::string fileB = fileToWrite + "_B.txt";
	std::string fileC = fileToWrite + "_C.txt";

	FILE *f = fopen(fileA.c_str(), "w");
	if (f == NULL) {
		printf("Error opening file for writing\n");
		return;
	}
	for (int i = 0; i < dimension; i++) {
		fprintf(f, "%f %f\n", A(i, 0), A(i, 1));
	}
	fclose(f);

	f = fopen(fileB.c_str(), "w");
	if (f == NULL) {
		printf("Error opening file for writing\n");
		return;
	}
	for (int i = 0; i < dimension; i++) {
		fprintf(f, "%f\n", b[i]);
	}
	fclose(f);

	f = fopen(fileC.c_str(), "w");
	if (f == NULL) {
		printf("Error opening file for writing\n");
		return;
	}
	fprintf(f, "%f %f\n", c[0], c[1]);
	fclose(f);
}

int main(int argc, char *argv[]) {
	unsigned int LP_size = 1, avg;
	unsigned int MaxMinFlag=2; //1 for Min and 2 for Max
	math::matrix<double> A;
	math::matrix<double> C, newC;
	std::vector<double> b, c;
	std::vector<double> result;
	std::vector<int> status_val;
	std::string afile, bfile, cfile, efile;



	if (argc != 2) { //1(ApplicationName) + 5 (Input Arguments)
		std::cout << "\nInsufficient Number of Arguments!!!\n";
		std::cout << "Correct Usages/Syntax:\n";
		std::cout << "./ProjName   'Benchmark-NO'\n";
		std::cout << "Argument 1) Benchmark-Name -- The size of problem to generate, e.g. '256'\n";
		return 0;
	}

	printf("Generating...\n");
	generateRandomLP2D(atoi(argv[1]), A, b, c, MaxMinFlag);

	WriteToFile(atoi(argv[1]), A, b, c);
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


