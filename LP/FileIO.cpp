#include "FileIO.hpp"
#include <iostream>
#include <fstream>
#include <string>

#define PI 3.14159265358979323846264338327950

using namespace std;

bool parseBenchmark(float4** constraints, glm::vec2 *optimisation, int *size) {

	//expects files of form:
	//_A.txt	first line contain ("%i %i", number of constraints, dimensions = 2)
	//			rest of file containing equations of form A_i*x_i < B. i sums over dimensions
	//B file	containing the appropriate value for B for above equation
	//C file	containing objective function parameters x & y to MINIMIZE

	char *fileA = "benchmarks/Test_A.txt"; char *fileB = "benchmarks/Test_B.txt"; char *fileC = "benchmarks/Test_C.txt";
	int dimTest; //check dimensions of input file
	vector<glm::vec2> A; //values in A file
	vector<float> b; //values in B file
	//glm::vec2 c; //2 values in C file
	int counter = 0; //checks loop iterations

#ifdef FILEIO_PRINT
	//print absolute path
	char absPath[100];
	_fullpath(absPath, fileA, 100);
	printf("Absolute path expected for file A is %s\n", absPath);
#endif // FILEIO_PRINT


	///////////////////////////////////////
	//Read in from files
	//File A
	ifstream f(fileA);
	if (!f.is_open()) {
		cout << "unable to open file A" << endl;
		return false;
	}
	//get number of constraints (in size) and double check file is in 2D
	f >> *size >> dimTest;
	if (dimTest != 2) {
		cout<< "Expected 2 dimensions for input file!" << endl;
	}

	//extend memory
	*constraints = (float4*)malloc(sizeof(float4) * (*size));
	A.reserve(*size);
	b.reserve(*size);


	//read remainder of file
	float xelem = 0, yelem = 0;
	while (f >> xelem >> yelem) {
		A.push_back(glm::vec2(xelem, yelem));
		counter++;
	}
	f.close();
	if (counter != *size) {
		std::cout << "Warning: given number of constraints and expected number of constraints do not match! " << counter << " vs " << *size << endl;
	}
	counter = 0; //reset counter;

	//File B
	f.open(fileB);
	if (!f.is_open()) {
		cout << "unable to open file B" << endl;
		return false;
	}
	while (f >> xelem) {
		b.push_back(xelem);
		counter++;
	}
	f.close();
	if (counter != *size) {
		std::cout << "Warning: given number of constraints and expected number of constraints do not match! " << counter << " vs " << *size << endl;
	}
	counter = 0; //reset counter;

	//File C
	f.open(fileC);
	if (!f.is_open()) {
		cout << "unable to open file C" << endl;
		return false;
	}
	f >> xelem >> yelem;
	*optimisation = glm::vec2(xelem, yelem);
	f.close();

	////////////////////////////
	//Convert A and b from Ax < b to lineDir, linePoint
	for (int i = 0; i < (*size); i++) {
		convertLine3to4(&((*constraints)[i]), A[i], b[i]);
	}

	return true;
}

/*
* Converts line of form Ax < b to lineDir & linePoint (in @fourVar)
*/
void convertLine3to4(float4 * fourVar, glm::vec2 A, float b) {
	float ycomp = A.y;
	float xcomp = A.x;
	float x_intc = 0; //does not exist to be read
	float y_intc = b;

	float dirx, diry, pointx, pointy;
	//vertical line
	if (ycomp == 0) {
		diry = xcomp;
		dirx = 0;
		pointx = y_intc / xcomp;
		pointy = 0;
	}
	//horizontal line
	else if (xcomp == 0) {
		diry = 0;
		dirx = -ycomp;
		pointx = 0;
		pointy = y_intc / ycomp;
	}
	//usual case
	else {
		diry = 1;
		dirx = -xcomp / ycomp;
		pointx = 0;
		pointy = y_intc / ycomp;
	}

	//point dirx diry in right direction
	if (ycomp < 0) {
		//dirx is positive
		dirx = (dirx < 0) ? -dirx : dirx;
	}
	else {
		//dirx negative
		dirx = (dirx < 0) ? dirx : -dirx;
	}
	if (xcomp < 0) {
		//diry is negative
		diry = (diry < 0) ? diry : -diry;
	}
	else {
		//diry is positive
		diry = (diry < 0) ? -diry : diry;
	}

	//save
	*fourVar = make_float4(dirx, diry, pointx, pointy);

}

/*
* Converts line of form lineDir & linePoint (in @fourVar) to Ax < b
*/
void convertLine4to3(float4 fourVar, glm::vec2 *A, float *b) {
	float m = fourVar.y / fourVar.x;
	float c = fourVar.w - (fourVar.z * m);

	*A = glm::vec2(m, 1);
	*b = c;
}


void writeLPtoFiles(float4 *h_lines, glm::vec2 optimisation, int size, const char* const name){
	//name in string form
	string s_name = string(name);

	printf("Writing lp to file...\n");
	glm::vec2 *A = (glm::vec2*)malloc(sizeof(glm::vec2) * size);
	float *B = (float*)malloc(sizeof(float) * size);

	//convert dirxy and pointxy to a11 x1 + a12 x2 <= b1
	for (int i = 0; i < size; i++) {
		convertLine4to3(h_lines[i], &A[i], &B[i]);
	}

	//write to file

	FILE *f = fopen((s_name + string("_A.txt")).c_str(), "w");
	if (f == NULL) {
		printf("Error opening file for writing\n");
		return;
	}
	fprintf(f, "%i %i\n", size, 2);
	for (int i = 0; i < size; i++) {
		fprintf(f, "%f %f\n", A[i].x, A[i].y);
	}
	fclose(f);

	f = fopen((s_name + string("_B.txt")).c_str(), "w");
	if (f == NULL) {
		printf("Error opening file for writing\n");
		return;
	}
	for (int i = 0; i < size; i++) {
		fprintf(f, "%f\n", B[i]);
	}
	fclose(f);

	f = fopen((s_name + string("_C.txt")).c_str(), "w");
	if (f == NULL) {
		printf("Error opening file for writing\n");
		return;
	}
	fprintf(f, "%f\n%f\n", optimisation.x, optimisation.y);
	fclose(f);
	printf("Writing of LP to file completed\n");
}




void generateRandomLP(float4** lines, glm::vec2* optimisation, const int size){
	//y-offset of line to solution approximated
	float yOffset = 50;

	//choose arbitary x,y as solution
	float xsol = randF() * 100, ysol = randF() * 100;
	*optimisation = glm::vec2(xsol, ysol);

#ifdef FILEIO_PRINT
	//print solution
	printf("Sln x: %f y:%f\n", xsol, ysol);
#endif

	//memory allocation
	*lines = (float4*)malloc(sizeof(float4) * size);

	//loop over all lps to generate
	for (int i = 0; i < size; i++) {

		//generate dir between -1 and 1
		float angle = randF() * PI;
		float dirx = cos(angle);
		float diry = sin(angle);

		//generate y-offset of line to solution
		float offset = -randF() * yOffset ;
		//if line passes in x direction, ensure line has smaller y-intercept value
		if (dirx > 0)
			offset *= -1;

		//a point on the line
		float pointx = xsol;
		float pointy = ysol + offset;

		//store it
		(*lines)[i].x = dirx;
		(*lines)[i].y = diry;
		(*lines)[i].z = pointx;
		(*lines)[i].w = pointy;

#ifdef FILEIO_PRINT
		//print ORCA
		printf("\ndirx: %f diry: %f \t px: %f py: %f\n", dirx, diry, pointx, pointy);
#endif

		////check solution
		//float res = m * xsol + c;
		//if (dirx > 0) {
		//	if (ysol < res)
		//		printf("%i calculated incorrect\n", i);
		//}
		//else {
		//	if (ysol > res)
		//		printf("%i calculated incorrect\n", i);
		//}
	}
}


int writeTimingtoFile(const char* const name, const int size, const int batches, const float time){
	FILE *f = fopen(name, "a");
	if (f == NULL) {
		printf("Error opening timing file for appending\n");
		return 1;
	}
	fprintf(f, "%i\t%i\t%f\n", size, batches, time);
	fclose(f);
	printf("Wrote to timing file data successfully\n");
	return 0;
}
