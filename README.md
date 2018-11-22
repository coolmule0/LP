# Random Batch LP Solver

This is a batch linear program solver for Nvidia GPUs, suitable for solving many 2D problems efficiently. 

To test results from paper "Two-Dimensional Batch Linear Programming on the GPU" please see the "Improvements" branch.

## How to build

### CMake

Use CMake to generate the relevant makefile, VS solution, e.t.c. Useful CMake flags to set would be adding `-arch=sm_xx` to "CMAKE_CUDA_FLAGS". This ensures that the project is built using a suitable architecture. Requires `xx` to be set to at least 30. Therefor ensure the GPU this is ran on is an NVidia graphics card with compute capability of at least 30.

### Windows

Directory contains a visual studio solution not requiring CMake. Open this .sln file in visual studio.

## What this project contains

CMake capabilities to generate a build tool such as make or .sln. The GLM library in `include/glm`

## How to Run

Build and compile the project to generate an executable. The default location will be `x64/"build_mode"/LP` where "build_mode" is either Debug or Release.

Executable takes 2 arguments: 1) The name of constraint data to read in. 2) The number of batches of these to solve.

An example would be `./x64/Release/LP benchmarks/64 1024` to use the pregenerated file containing 64 random constraints with feasible solution, for 1024 batches.

## Program Output

The program will output timing to the "timings" folder in a timings.txt file. This file contains lines of the form `%1 %2 %3\n` where %1 is the number of constraints (e.g. 64), %2 is the number of batches (e.g. 1024), %3 is the time in milliseconds that the program took.

The solutions to the linear program are printed to console. Access to these values is accessible within the code by both GPU and CPU code in `output[i]` where i ranges over all batches.

## Custom Data

To pass custom data to the executable create 3 files of the form *_A.txt, *_B.txt, *_C.txt. Where * is a consistent name for all 3 files. File *_A.txt starts with a line of 2 numbers, the first is the number of constraints within the file (n) and the dimension of the problem (d) (should always be 2). The remainder of the file is (n) lines of (d) numbers. These represent the left hand side of the constraints A1x + A2y <= B.
File *_B.txt contains (n) lines of 1 number. These represent the right hand side of the constraints A1x + A2y <= B. File *_C.txt contains (d) lines of 1 number. These represent the objective function to maximise. Ensure data is written with respects to maximising objective function, and inequalities are of form Ax <= B.