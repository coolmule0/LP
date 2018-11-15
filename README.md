# Random Incremental Batch LP Solver for GPU

This branch is for running the comparison tests from the paper "Two-Dimensional Batch Linear Programming on the GPU"

## How to build

### CMAKE

Ensure the configuration is set for x64 architecture. Different sub-projects have different dependencies. In total these are BOOST, CLP and GLPK.

Ensure you set relevant CUDA compute version by setting the flag e.g. `CMAKE_CUDA_FLAGS=-arch=sm_52`. This sets it to compute 52. Requirement is at least 30. Please ensure the computer this is running on has a supported Nvidia GPU.

BOOST: set the flags `BOOST_INC` and `BOOST_LIB` to the boost source directory (containing boost/boost.h as a sub file) and to the boost library directory respectively.

CLP: set `CLP_INC` to include directory, containing coin/CLPSimplex.hpp as a subfile. Set `CLP_LIB` to the folder containing CLP library files.

GLPK: set flags `GLPK_INC` to GLPK source / src folder. `GLPK_LIB` to the GLPK library / w64 folder.

Ensure release mode is built by setting flag `CMAKE_BUILD_TYPE=Release`

On windows, ensure CudaToolkitBinDir is an environmental variable so runtime .dll's can be passed to the generated executable location. Failing this, manually copy cudart.dll to the executable location.

Successful setting of CMAKE flags will generate the appropriate MAKE file/Visual Studio .sln (e.t.c.). Building this will generate the executables.



## Running examples

CLP: Takes 2 executable arguments: 1) The name of constraint data to read in. 2) The number of batches of these to solve.

An example would be `LP.exe benchmarks/64 1024` to use the file containing 64 randomly generated constraints for 1024 batches.

CPLEX: Takes 2 executable arguments: 1) The name of constraint data to read in. 2) The number of batches of these to solve.

LP: Takes 2 executable arguments: 1) The name of constraint data to read in. 2) The number of batches of these to solve.

mGLPK: Takes 5 executable arguments: 1) The name of constraint data to read in. 3) The number of batches of these to solve.



## Program Output

The program will output timing to the "timings" folder in a timings.txt file. This file contains lines of the form `%1 %2 %3\n` where %1 is the number of constraints (e.g. 64), %2 is the number of batches (e.g. 1024), %3 is the time in milliseconds that the program took.

The solutions to the linear program are printed to console. Access to these values is accessible by both GPU and CPU in `output[i]` where i ranges over all batches.

## Using Custom Data

To use a problem set to solve using the programs, create 3 files of the form *_A.txt, *_B.txt, *_C.txt. Where * is a consistent name for all 3 files.
File *_A.txt starts with a line of 2 numbers, the first is the number of constraints within the file (n) and the dimension of the problem (d) (should always be 2). The remainder of the file is (n) lines of (d) numbers. These represent the left hand side of the constraints A1x + A2y <= B.

File *_B.txt contains (n) lines of 1 number. These represent the right hand side of the constraints A1x + A2y <= B

File *_C.txt contains (d) lines of 1 number. These represent the objective function to maximise.

Ensure data is written with respects to maximising objective function, and inequalities are of form Ax <= B.

## Remaking Manuscript Figures

A bash .sh script is included in the base directory, called `ReproduceResults.sh`, which runs all the executables over a variety of problem and batch sizes. Results are written to the "timings" directory. 

These timings can be coalesced into a graph using the Matlab script found in the "matlab" subdirectory called `LPtimingScript.m`.

Running `ReproduceResults.sh` followed by `LPtimingScript.m` should generate the figures presented in the paper, appropriate to the device being tested on.

