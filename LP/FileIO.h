#pragma once
//#include <cuda_runtime.h>
//#include "device_launch_parameters.h"
//#include "cuda_profiler_api.h"
#include <glm/glm.hpp>

#include <vector>
#include <string>
//#include "Auxilary.h"


/*
* Functions to read/write LP data
*/

//Whether to print info for FIleIO functions
//#define FILEIO_PRINT

/*
* Read in line constraint information from file and write data to newly malloced @ORCA_line
* Reads in data for 1 LP
*
* @filename: name of file to read in data from, from benchmarks folder
* @ORCA_line: array which will be allocated CPU memory and contain constraint data
* @optimisation: 2D objective function variables to minimise over.
* @size: number of constraints
*/
bool parseBenchmark4(std::string filename, glm::vec4** ORCA_line, glm::vec2 *optimisation, int *size);

/*
* Read in line constraint information from file and write data to newly allocated A and b
* Reads in data for 1 LP
*
* @filename: name of file to read in data from, from benchmarks folder
* @A: array which will be allocated CPU memory and contain constraint data
* @b: array which will be allocated CPU memory and contain constraint data
* @optimisation: 2D objective function variables to minimise over.
* @size: number of constraints
*/
bool parseBenchmark3(std::string filename, std::vector<glm::vec2> &A, std::vector<float> &b, glm::vec2 *optimisation, int *size);

/*
* Converts line of form Ax < b to lineDir & linePoint (in @fourVar)
*/
void convertLine3to4(glm::vec4 * fourVar, glm::vec2 A, float b);




/*
* Write data of h_lines to appropriate text file
*
* @h_lines: array to write to file
* @size: number of constraints
* @name: path and file name. Will be postfixed with _A, _B and _C .txt
*/
void writeLPtoFiles(glm::vec4 *h_lines, const int size, const char* const name);


////////////////////////////////////
//Linear program generation

/* generates an LP with feasible solution
* Line is visually represented as:
		 ^
		/* * * * *
Valid  / *Invalid*
	  /* * * * * *
	 / * * * * * *
	/* * * * * * *
*/
//void generateRandomLP(glm::vec4** ORCA_line, glm::vec2* optimisation, const int size);

/* Write time take in ms from @time with batch size @size and number of batches @batches to file @name
 *
 */
 int writeTimingtoFile(const char* const name, const int size, const int batches, const float time);


 /* Convert file_A _B and _C into a MPS file type of the same name, with .mps extension
  */
 void convertToMPS(std::string inname, std::string outname, int &size);

 /*Create a run file, which can be easily passed to CPLEX to solve
  *filename of file to optimise, with extension
  *creates a CPLEX run file called "CPLEXrunfile"
  */
 void genrunfile(std::string filename);