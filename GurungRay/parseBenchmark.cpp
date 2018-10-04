#include "parseBenchmark.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
using namespace std;

/*
 * Minimize cx
 * sub to   Ax <=b
 *    x>= 0
 *
 * NetLib Benchmarks
 *                       PROBLEM SUMMARY TABLE
   Name       Rows   Cols   Nonzeros    Bytes        Optimal Value
*************************** Smaller Dimension LP **************************
***ADLITTLE     57     97      465       3690        2.2549496316E+05
***AFIRO        28     32       88        794       -4.6475314286E+02
***BLEND        75     83      521       3227       -3.0812149846E+01
***ISRAEL      175    142     2358      12109       -8.9664482186E+05
***SC105       106    103      281       3307       -5.2202061212E+01
***SC205       206    203      552       6380       -5.2202061212E+01
***SC50A        51     48      131       1615       -6.4575077059E+01
***SC50B        51     48      119       1567       -7.0000000000E+01
*
*/

void selectBenchmark(int benchmarkNo, std::string& argvA, std::string&  argvB, std::string& argvC){
	if (benchmarkNo==1){//ADLITTLE  has solution as 2.2549496316E+05
		argvA = "benchmarks/ADLITTLE_A.txt"; argvB = "benchmarks/ADLITTLE_B.txt";	argvC = "benchmarks/ADLITTLE_C.txt";
		std::cout<<"Benchmark = ADLITTLE"<<std::endl;
	} else if (benchmarkNo==2){ //AFIRO has solution as -4.6475314286E+02
		argvA="benchmarks/AFIRO_A.txt";  argvB="benchmarks/AFIRO_B.txt"; argvC="benchmarks/AFIRO_C.txt";
		std::cout<<"Benchmark = AFIRO"<<std::endl;
	} else if (benchmarkNo==3){//BLEND  has solution as -3.0812149846E+01
		argvA="benchmarks/BLEND_A.txt";  argvB="benchmarks/BLEND_B.txt"; argvC="benchmarks/BLEND_C.txt";
		std::cout<<"Benchmark = BLEND"<<std::endl;
	} else if (benchmarkNo==4){//ISRAEL  has solution as -8.9664482186E+05
		argvA="benchmarks/ISRAEL_A.txt";  argvB="benchmarks/ISRAEL_B.txt"; argvC="benchmarks/ISRAEL_C.txt";
		std::cout<<"Benchmark = ISRAEL"<<std::endl;
	} else if (benchmarkNo==5){//SC105  has solution as -5.2202061212E+01
		argvA="benchmarks/SC105_A.txt";  argvB="benchmarks/SC105_B.txt"; argvC="benchmarks/SC105_C.txt";
		std::cout<<"Benchmark = SC105"<<std::endl;
	} else if (benchmarkNo==6){//SC205  has solution as -5.2202061212E+01
		argvA="benchmarks/SC205_A.txt";  argvB="benchmarks/SC205_B.txt"; argvC="benchmarks/SC205_C.txt";
		std::cout<<"Benchmark = SC205"<<std::endl;
	} else if (benchmarkNo==7){//SC50A  has solution as -6.4575077059E+01
		argvA="benchmarks/SC50A_A.txt";  argvB="benchmarks/SC50A_B.txt"; argvC="benchmarks/SC50A_C.txt";
		std::cout<<"Benchmark = SC50A"<<std::endl;
	} else if (benchmarkNo==8){//SC50B  has solution as -7.0000000000E+01
		argvA="benchmarks/SC50B_A.txt";  argvB="benchmarks/SC50B_B.txt"; argvC="benchmarks/SC50B_C.txt";
		std::cout<<"Benchmark = SC50B"<<std::endl;
	} else if (benchmarkNo == 9) {
		argvA = "../LP/LP/benchmarks/Test_A.txt";  argvB = "../LP/LP/benchmarks/Test_B.txt"; argvC = "../LP/LP/benchmarks/Test_C.txt";
		std::cout << "Benchmark = Test" << std::endl;
	} else if (benchmarkNo == 10) {
		argvA = "../../LP/LP/benchmarks/Test_A.txt";  argvB = "../../LP/LP/benchmarks/Test_B.txt"; argvC = "../../LP/LP/benchmarks/Test_C.txt";
		std::cout << "Benchmark = Test" << std::endl;
	}
	else {
		std::cout << "Invalid benchmark no" << std::endl;
	}

	//get absolute path
	char absPath[100];
	_fullpath(absPath, argvA.c_str(), 100);
	std::cout << "Absolute path is " << absPath << "\n";
}

void selectBenchmarkNamed(const std::string inName, std::string& argvA, std::string&  argvB, std::string& argvC) {
	argvA = inName + "_A.txt";
	argvB = inName + "_B.txt";
	argvC = inName + "_C.txt";

	////get absolute path
	//char absPath[100];
	//_fullpath(absPath, argvA.c_str(), 100);
	//std::cout << "Absolute path is " << absPath << "\n";
}


void parseLP(const char* argvA,const  char* argvB,const  char* argvC, math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag ){
	/*std::cout << "Argument 1) Afile -- "<<argvA <<std::endl;
	std::cout << "Argument 2) Bfile -- "<<argvB <<std::endl;
	std::cout << "Argument 3) Cfile -- "<<argvC <<std::endl;*/

	int row = 0, col = 0;
	double field;
	int tot_constraints, tot_variables;
	std::ifstream in(argvA);
	if (in) {
		std::string line1;
		std::getline(in, line1); //Reading First Line for total Constraints and Variables
		std::istringstream row1(line1);
		row1 >> tot_constraints; //Reading First value
		row1 >> tot_variables; //Reading Second value
		std::cout << "Total Constraints =" << tot_constraints << "   ";
		std::cout << "Total Variables =" << tot_variables << std::endl;
		A.resize(tot_constraints, tot_variables);
		b.resize(tot_constraints);	//total rows or constraints same as A
		c.resize(tot_variables);	//total number of variables same as A
		//reading remaining Lines
		row = 0; //First Constraint
		while (std::getline(in, line1)) {
			col = 0;
			std::istringstream row1(line1);
			while (row1 >> field) {
				if (field==0 || field==-0){
					//std::cout<< field <<"  ";
					A(row, col) = (double)0.0;
				} else
					A(row, col) = field;
				//cout << field <<"\t";
				col++;
			}
			row++; //Next Constraints
			//cout<<"\n";
		}
	}
	else {
		cerr << "File " << argvA << " could not be opened!\n"; // Report error
		cerr << "Error code: " << strerror(errno) << "\n"; // Get some info as to why
	}//Reading A matrix is Over
	std::ifstream inB(argvB); //reading Bfile
	if (inB) {
		std::string line1;
		row = 0; //First Constraint
		while (std::getline(inB, line1)) {
			std::istringstream row1(line1);
			row1 >> field;
			//std::cout<< field <<"  ";
			b[row] = field;
			//cout << field <<"\t";
			row++; //Next Constraints
			//cout<<"\n";
		}
	}
	else {
		cerr << "File " << argvB << " could not be opened!\n"; // Report error
		cerr << "Error code: " << strerror(errno) << "\n"; // Get some info as to why
	}//Reading B matrix is Over
	std::ifstream inC(argvC); //reading Cfile
	if (inC) {
		std::string line1;
		row = 0; //First Constraint
		while (std::getline(inC, line1)) {
			std::istringstream row1(line1);
			row1 >> field;
			//std::cout<< field <<"  ";
			c[row] = field;
			row++; //Next Constraints
			//cout<<"\n";
		}
	}
	else {
		cerr << "File " << argvC << " could not be opened!\n"; // Report error
		cerr << "Error code: " << strerror(errno) << "\n"; // Get some info as to why
	}//Reading C matrix is Over

	MaxMinFlag = 1; //Minimize Netlib benchmarks are all Minimization problem. //1=min 2=max
}

/*
 * Max 		2x1 - x2
 * sub to   2x1 - x2 <= 2
 * 			x1 - 5x2 <= -4
 * 		x1, x2 >= 0
 *
 *
 *  //Unbounded problem has many solutions at 2
 */

void parseSmallLP(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){
	A.resize(2,2);
	b.resize(2);
	c.resize(2);
	A(0,0)=2;	A(0,1)=-1;
	A(1,0)=1;	A(1,1)=-5;
	b[0]=2; b[1]=-4;
	c[0]=2; c[1]=-1;

	MaxMinFlag = 2; //Maximize
}


/*
 * Max 		5x1 + 6x2
 * sub to   x1 + x2 <= 10
 * 			x1 - x2 >= 3 or -x1 +x2 <= -3
 * 			5x1 + 4x2 <= 35
 * 		x1, x2 >= 0
 *
 *
 * 		//has solution 39.44444
 */

void parseSmallLP2(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){
	A.resize(3,2);
	b.resize(3);
	c.resize(2);
	A(0,0)=1;	A(0,1)=1;
	A(1,0)=-1;	A(1,1)=1;
	A(2,0)=5;	A(2,1)=4;
	b[0]=10; b[1]=-3; b[2]=35;
	c[0]=5; c[1]=6;

	MaxMinFlag = 2; //Maximize
}

/*
 * Max 		x1 + 3x2
 * sub to   x1 - x2 <= 8
 * 			-x1 - x2 <= -3
 * 			-x1 + 4x2 <= 2
 * 		x1, x2 >= 0
 *
 *
 * 		//has solution 21.33333
 */

void parseSmallLP3(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){
	A.resize(3,2);
	b.resize(3);
	c.resize(2);
	A(0,0)=1;	A(0,1)=-1;
	A(1,0)=-1;	A(1,1)=-1;
	A(2,0)=-1;	A(2,1)=4;
	b[0]=8; b[1]=-3; b[2]=2;
	c[0]=1; c[1]=3;

	MaxMinFlag = 2; //Maximize
}

/*
 * Max 		2x1 - 6x3
 * sub to   x1 + x2 - x3 <= 7
 * 			3x1 - x2 >= 8        or  -3x1 + x2 <= -8
 * 			-x1 + 2x2 + 2x3 >= 0 or  x1 - 2x2 - 2x3 <= 0
 * 		x1, x2 >= 0
 *
 * 		//has solution 9.33333
 */

void parseSmallLP4(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){
	A.resize(3,3);
	b.resize(3);
	c.resize(3);
	A(0,0)=1;	A(0,1)=1;	A(0,2)=-1;
	A(1,0)=-3;	A(1,1)=1;	A(1,2)=0;
	A(2,0)=1;	A(2,1)=-2;	A(2,2)=-2;

	b[0]=7; b[1]=-8; b[2]=0;
	c[0]=2; c[1]=0;  c[2]=-6;

	MaxMinFlag = 2; //Maximize
}

/*
 * Max 		x1 + x2
 * sub to   4x1 - x2 <= 8
 * 			2x1 + x2  <= 10
 * 			5x1 - 2x2 >= -2    or  -5x1 + 2x2 <= 2
 * 		x1, x2 >= 0
 *
 * 		//has solution 39.44444
 */

void parseSmallLP5(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){

	A.resize(3,2);
	b.resize(3);
	c.resize(2);
	A(0,0)=4;	A(0,1)=-1;
	A(1,0)=2;	A(1,1)=1;
	A(2,0)=-5;	A(2,1)=2;
	b[0]=8; b[1]=10; b[2]=2;
	c[0]=1; c[1]=1;

	MaxMinFlag = 2; //Maximize
}


/*
 * Continuous System Example: A 5 Dimensional System
 *** We take the initial input convex set
 */
void parse5Dim(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){
	int row,col;
	row = 10;	col = 5;
	A.resize(row, col);
	A(0, 0) = 1;	A(0, 1) = 0;	A(0, 2) = 0;	A(0, 3) = 0;	A(0, 4) = 0;
	A(1, 0) = -1;	A(1, 1) = 0;	A(1, 2) = 0;	A(1, 3) = 0;	A(1, 4) = 0;
	A(2, 0) = 0;	A(2, 1) = 1;	A(2, 2) = 0;	A(2, 3) = 0;	A(2, 4) = 0;
	A(3, 0) = 0;	A(3, 1) = -1;	A(3, 2) = 0;	A(3, 3) = 0;	A(3, 4) = 0;
	A(4, 0) = 0;	A(4, 1) = 0;	A(4, 2) = 1;	A(4, 3) = 0;	A(4, 4) = 0;
	A(5, 0) = 0;	A(5, 1) = 0;	A(5, 2) = -1;	A(5, 3) = 0;	A(5, 4) = 0;
	A(6, 0) = 0;	A(6, 1) = 0;	A(6, 2) = 0;	A(6, 3) = 1;	A(6, 4) = 0;
	A(7, 0) = 0;	A(7, 1) = 0;	A(7, 2) = 0;	A(7, 3) = -1;	A(7, 4) = 0;
	A(8, 0) = 0;	A(8, 1) = 0;	A(8, 2) = 0;	A(8, 3) = 0;	A(8, 4) = 1;
	A(9, 0) = 0;	A(9, 1) = 0;	A(9, 2) = 0;	A(9, 3) = 0;	A(9, 4) = -1;

	b.resize(row);
	b[0] = 1.01;	b[1] = -0.99;	b[2] = 0.01;	b[3] = 0.01;	b[4] = 0.01;
	b[5] = 0.01;	b[6] = 0.01;	b[7] = 0.01;	b[8] = 0.01;	b[9] = 0.01;

	c.resize(col);
	//c[0]=1; c[1]=0;  c[2]=0;	c[3]=0;	c[4]=0;		//Solution is 1.01
	//c[0]=-1; c[1]=0;  c[2]=0;	c[3]=0;	c[4]=0;		//Solution is -0.99
	c[0]=0; c[1]=1;  c[2]=0;	c[3]=0;	c[4]=0;		//Solution is -0.99

	MaxMinFlag = 2; //Maximize
}

/*
 * Continuous System Example: Input for a 28 Dimensional Helicopter Controller System
 */
void parse28DimHelicopter(math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){
	int row,col;
	row = 56;
	col = 28;//constraints for each 28variables but excluding the time-constraint
	A.resize(row, col);
	b.resize(row);

	for (int i = 0; i < col; i++) {
		for (int j = 0; j < col; j++) {
			if (i == j) {	//for diagonal elements assign 1 and -1
				A(2 * i, j) = 1;			//for xi
				A(2 * i + 1, j) = -1;		//for -xi
			} else {	//for all other elements assign zeros
				A(2 * i, j) = 0;			//for xi
				A(2 * i + 1, j) = 0;		//for -xi
			}
		}
		if (i < 8) {		//x1  to x8
			b[2 * i] = 0.1;
			b[2 * i + 1] = 0;
		}
		if ((i >= 8) && (i < col)) {		//x9  to x28
			b[2 * i] = 0;
			b[2 * i + 1] = 0;
		}
	}

	c.resize(col);

	for (unsigned int j = 0; j < 28; j++) {
		if (j == 0)
			c[j] = 1;
		else
			c[j] = 0;
	}

	MaxMinFlag = 2; //Maximize
}

void generateRandomLP(int dimension, math::matrix<double>& A, std::vector<double>& b, std::vector<double>& c, unsigned int& MaxMinFlag){
	unsigned int status = 0, res1 = 0, res2 = 0;
	// *************************************************************************************************************************************
	unsigned int art=0; //Decides the number of artificial variables (or make the initial basic solution of LP as Feasible or Infeasible)
	// *************************************************************************************************************************************
	time_t t;
	/* Intializes random number generator */
	   srand((unsigned) time(&t));

	int sign = -1, sum1;
	glpk_lp_solver lp;

	while (status != 5) {
		//We are assuming that LP is in the form A<=b
		A.resize(dimension, dimension);
		b.resize(dimension);
		for (unsigned int j = 0; j < dimension; j++) {
			for (unsigned int k = 0; k < dimension; k++) {
				//A(j, k) = rand() % (k + 10) + 1;
				A(j, k) = (rand() % 999) + 1; // Range 1 to 1000
			}
			if (j < art) {//for LPs with initial basic solution is Infeasible
				res1 = (rand() % 999);
				res2 = (res1 + 1);
				sum1 = res2 * sign;
				b[j] = (double) sum1;
				A(j, j) = A(j, j) * sign; //This is required to make the LP feasible -A(j,j)<=-b[j]
			} else	//for LPs with initial basic solution is  Feasible
				//b[j] = (rand() % (j + 1) + (10 + j));
				b[j] = (rand() % 999) + 1;  // Range 1 to 1000
		}
		//** Setting current Lp to GLPK
		lp.setMin_Or_Max(2);
		lp.setConstraints(A, b, 1);
		status = lp.TestConstraints(); //The objective function is all set ot zero. //std::cout << "Status = " << status << "\n";
	}
	/*
	 * Experimental Note:
	 *  We generate randomly the first feasible LP problem with
	 *  Matrix A (in the range 1 to 1000)
	 *  Vector b (in the range 1 to 1000) and
	 *  The objective function c (in the range 1 to 500).
	 *  We use the system timestamp as the seed for the random generator
	 */

	MaxMinFlag = 2; //Maximize
	c.resize(dimension);
	unsigned int i;
	for (unsigned int j = 0; j < c.size(); j++) {
		//c[j] = rand() % (j + 1) + 1;
		c[j] = (rand() % 499) + 1; //Range 1 to 500
	}
	double res=lp.Compute_LLP(c);

	std::cout<<"GLPK::Solution = "<<res<<std::endl;
}
