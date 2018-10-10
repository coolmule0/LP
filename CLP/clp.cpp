
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.

#include "ClpSimplex.hpp"
#include "../LP/FileIO.h"

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

	//create and MPS of the file
	convertToMPS(argv[1], argv[1], size);
	std::string filename = std::string(argv[1]) + ".mps";

	ClpSimplex  model;
	int status;
	status = model.readMps(filename.c_str());
	if (!status) {
		model.primal();
		model.dual();
	}
	return 0;
}
