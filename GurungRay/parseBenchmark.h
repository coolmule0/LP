/*
 * parseBenchmark.h
 *
 *  Created on: 09-Jun-2017
 *      Author: amit
 */

#ifndef PARSEBENCHMARK_H_
#define PARSEBENCHMARK_H_

#include <vector>
#include <string>
#include "math/matrix.h"
#include "math/glpk_lp_solver/glpk_lp_solver.h"

//void selectBenchmark(int benchmarkNo, char* argvA, char* argvB, char* argvC);
void selectBenchmark(int benchmarkNo, std::string& argvA, std::string&  argvB, std::string& argvC);

void parseLP(const char* argvA,const char* argvB,const char* argvC, math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);

void parseSmallLP(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);
void parseSmallLP2(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);
void parseSmallLP3(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);
void parseSmallLP4(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);
void parseSmallLP5(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);

void parse5Dim(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);
void parse28DimHelicopter(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);

void generateRandomLP(int dimension, math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag);

#endif /* PARSEBENCHMARK_H_ */
