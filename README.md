# Random Batch LP Solver

This is a batch linear program solver for Nvidia GPUs. Works for GPUs greater than compute 2.1 (i.e. newer than Fermi architecture). Include glm library.

## How to build

### Windows

Open LP.sln solution file in visual studio. Using Visual Studio 2015 and Cuda 8.0.

### Linux

Run "make" from ./LP. By default this creates an executable called "RGBLP" in "./.". Tested on Ubuntu 16 for cuda 8.0.

## Running examples

Takes 2 executable arguments: 1) The name of constraint data to read in. 2) The number of batches of these to solve.

An example would be `LP.exe benchmarks/64 1024` to use the file containing 64 randomly generated constraints for 1024 batches.

## Program Output

The program will output timing to the "timings" folder in a timings.txt file. This file contains lines of the form `%1 %2 %3\n` where %1 is the number of constraints (e.g. 64), %2 is the number of batches (e.g. 1024), %3 is the time in milliseconds that the program took.

The solutions to the linear program are printed to console. Access to these values is accessible by both GPU and CPU in `output[i]` where i ranges over all batches.

## Custom Data
Create 3 files of the form *_A.txt, *_B.txt, *_C.txt. Where * is a consistent name for all 3 files.
File *_A.txt starts with a line of 2 numbers, the first is the number of constraints within the file (n) and the dimension of the problem (d) (should always be 2). The remainder of the file is (n) lines of (d) numbers. These represent the left hand side of the constraints A1x + A2y <= B.

File *_B.txt contains (n) lines of 1 number. These represent the right hand side of the constraints A1x + A2y <= B

File *_C.txt contains (d) lines of 1 number. These represent the objective function to maximise.

Ensure data is written with respects to maximising objective function, and inequalities are of form Ax <= B.

## Reproducing Results

Included are a bash and a batch script to run on Linux/Windows respectively. The script will run the executable numerous times. Currently a work in progress

## CPLEX

The same models can be ran using the CPLEX library. Currently only tested on windows. To build, ensure you set enviromental paths for: `CPLEX_STUDIO_DIR1271` to the base folder location of CPLEX.