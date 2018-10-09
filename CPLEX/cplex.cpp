#include <ilcplex/ilocplex.h>
#include "../LP/FileIO.h"
#include <glm/glm.hpp>
#include <chrono>
#include <stdlib.h>
#include <string>
#include <iostream>

//uncomment for useful prints to console
//#define PRINTINFO

//Use MPS file to run CPLEX. Comment out to run code below
//#define USEMPS

ILOSTLBEGIN

int main(int argc, const char* argv[]) {
	//int batches = 0; //number of LPs
	int size = 0; //size of each LP
	glm::vec2* x = NULL;  //< solution to randomly generated LP
	std::vector<glm::vec2> A; //array of constraints for 1 lp
	std::vector<float> b; //array of constraints for 1 lp
	glm::vec2 optimiseSingle; // variables to minimise in optimisation function for 1 lp
	std::vector<std::vector<glm::vec2>> Aall; // array of constraints
	std::vector<std::vector<float>> ball; // array of constraints
	std::vector<glm::vec2> optimise; //variable to minimise in optimisation function

	ofstream logfile("nul"); //logfile of CPLEX solver


	//------------------------------------------
	//handle args
	if (argc != 3) {
		printf("\nIncorrect Number of Arguments!\n");
		printf("Correct Usages/Syntax:\n");
		printf("Argument 1) File-Name -- Name of input file in benchmarks folder. Cannot contain spaces\n");
		printf("Argument 2) Batch-size -- Number of LPs to be solved\n");
		return 1;
	}
	const unsigned int batches = atoi(argv[2]); //number of LPs

#ifdef USEMPS
	//create and MPS of the file
	convertToMPS(argv[1], argv[1], size);

	//generate a run file
	genrunfile(string(argv[1]) + ".mps;");

	//the cplex runfile command - suppress cplex output
	std::string runcommand = "cplex<CPLEXrunfile > nul";

#else
	Aall.reserve(batches);
	ball.reserve(batches);

	//------------------------------------------
	//handle input

	//Input is from file
	printf("Parsing input files... ");
	if (!parseBenchmark3(argv[1], A, b, &optimiseSingle, &size)) {
		return 1;
	}
	printf("Done\n");


	//------------------------------------------
	//Tile LP data to all other LPs
	for (unsigned int i = 0; i < batches; i++) {
		Aall.push_back(A);
		ball.push_back(b);
		for (int j = 0; j < size; j++) {
			//Aall[i].push_back( A[j] ); //deep copy constraint data
			//ball[i].push_back( b[j] ); //deep copy constraint data
		}
		optimise.push_back( optimiseSingle );
	}
#endif

	//------------------------------------------
	//initialize memory initialsation timings
	auto t1 = std::chrono::high_resolution_clock::now();


	//------------------------------------------
	//repeat runs
//#pragma omp parallel for num_threads(4)
	for (unsigned int n = 0; n < batches; n++) {
		try {

#ifdef USEMPS
			//run it in CPLEX
			system(runcommand.c_str());
#else
			//set up problem 
			IloEnv env;
			IloNumVarArray vars(env);
			vars.add(IloNumVar(env));
			vars.add(IloNumVar(env));

			IloModel model(env);
			model.add(IloMaximize(env, optimise[n].x * vars[0] + optimise[n].y * vars[1]));
			//add constraints
			for (int i = 0; i < size; i++) {
				model.add(Aall[n][i].x * vars[0] + Aall[n][i].y * vars[1] <= ball[n][i]);
			}

			//run problem
			IloCplex cplex(model);
			cplex.setOut(logfile);
			if (!cplex.solve()) {
				env.error() << "Failed to optimize LP." << endl;
				throw(-1);
			}

#ifdef PRINTINFO
			printf("iter %i\n", n);
			IloNumArray vals(env);
			env.out() << "Solution status = " << cplex.getStatus() << endl;
			env.out() << "Solution value = " << cplex.getObjValue() << endl;
			cplex.getValues(vals, vars);
			env.out() << "Values = " << vals << endl;
#endif // PRINTINFO

			env.end();
#endif
		}

		catch (IloException& e) {
			cerr << "Concert exception caught: " << e << endl;
		}
		catch (...) {
			cerr << "Unknown exception caught" << endl;
		}
	}
	

	//end time
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
	std::cout << "process took " << fp_ms.count() << " ms\n";

	//------------------------------------------
	//write timing to file
#ifdef USEMPS
	writeTimingtoFile("timings/CPLEXtimingsv2.txt", size, batches, fp_ms.count());
#else
	writeTimingtoFile("timings/CPLEXtimings.txt", size, batches, fp_ms.count());
#endif

	return 0;
}
