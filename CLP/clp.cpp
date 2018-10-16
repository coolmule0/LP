#include "ClpSimplex.hpp"
#include "../LP/FileIO.h"
#include <chrono>
#ifdef OMP
#include <omp.h>
#endif

int main(int argc, const char *argv[])
{
	int size = 0; //size of each LP

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

	//create and MPS of the file
	convertToMPS(argv[1], argv[1], size);
	std::string filename = std::string(argv[1]) + ".mps";

	ClpSimplex  model;
	int status;
	status = model.readMps(filename.c_str());
	if (status) {
		return 1;
	}

	//------------------------------------------
	//initialize memory initialsation timings
	auto t1 = std::chrono::high_resolution_clock::now();

	std::vector<ClpSimplex> models(batches);

	for (int n = 0; n < batches; n++) {

		models[n] = model;


		models[n].primal();
	}

	//end time
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
	std::cout << "process took " << fp_ms.count() << " ms\n";

	//------------------------------------------
	//write timing to file
	writeTimingtoFile("timings/CLP.txt", size, batches, fp_ms.count());


	return 0;
}
