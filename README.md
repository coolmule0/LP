# Random Batch LP Solver

This is a batch linear program solver for Nvidia GPUs.

## How to build

### Windows

Open LP.sln solution file in visual studio. Using Visual Studio 2015 and Cuda 8.0.

## Linux

Run "make".

## Running examples

Takes 2 executable arguments: 1) The name of constraint data to read in. 2) The number of batches of these to solve.

An example would be `benchmarks/64 1024` to use the file containing 64 randomly generated constraints for 1024 batches.

### Custom Data
Create 3 files of the form *_A.txt, *_B.txt, *_C.txt. Where * is a consistent name for all 3 files. 
File *_A.txt starts with a line of 2 numbers, the first is the number of constraints within the file (n) and the dimension of the problem (d) (should always be 2). The remainder of the file is (n) lines of (d) numbers. These represent the left hand side of the constraints A1x + A2y <= B.

File *_B.txt contains (n) lines of 1 number. These represent the right hand side of the constraints A1x + A2y <= B

File *_C.txt contains (d) lines of 1 number. These represent the objective function to maximise.

Ensure data is written with respects to maximising objective function, and inequalities are of form Ax <= B.

