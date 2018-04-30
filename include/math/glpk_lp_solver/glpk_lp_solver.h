/*
 * lp_solver.h
 *
 *  Created on: 23-Apr-2014
 *      Author: gurung
 */

#ifndef GLPK_SOLVER_H_
#define GLPK_SOLVER_H_

#include <iostream>
#include <vector>
#include <glpk.h>
#include <boost/shared_ptr.hpp>

#include "math/matrix.h"

/* solution status:
 *
 #define GLP_UNDEF          1  		// solution is undefined
 #define GLP_FEAS           2  		// solution is feasible
 #define GLP_INFEAS         3  		// solution is infeasible
 #define GLP_NOFEAS         4  		// no feasible solution exists
 #define GLP_OPT            5  		// solution is optimal
 #define GLP_UNBND          6  		// solution is unbounded
 */

class glpk_lp_solver {
public:
	typedef boost::shared_ptr<glpk_lp_solver> glpk_ptr;

	void setDefaultObject(); // same as calling a constructor
	/*
	 glpk_lp_solver(std::vector<std::vector <double> > coeff_constraints, std::vector <double> bounds,
	 std::vector <int> bound_signs);
	 void setConstraints(std::vector<std::vector <double> > coeff_constraints, std::vector <double> bounds,
	 int bound_signs);
	 */
	glpk_lp_solver(math::matrix<double> coeff_constraints,
			std::vector<double> bounds, std::vector<int> bound_signs);
	void setConstraints(math::matrix<double> coeff_constraints,
			std::vector<double> bounds, int bound_signs);	//Before calling setConstraints(), the member
															//function setMin_or_Max() must be called.

	/**
	 * Adds the variables and constraints of the polytope to the current LP problem.
	 * Note that this function is different from just adding constraints since joining
	 * a polytope implies adding its variables as new variables and also adding constraints
	 * on these new variables given as parameters.
	 */
	void join_poly_constraints(math::matrix<double> coeff_constraints,
			std::vector<double> bounds, int bound_signs);

	/**
	 * Adds n new column variable to the problem
	 */
	void add_new_cols(int n);
	/**
	 * Adds a new constrain with the problem, in the ax <= b form
	 * len is the number of non-zero entries in the constraint
	 * ind is the array of column indices of the coefficient in the lp coeff matrix
	 * val is the corresponding array of values of the coefficient matrix
	 */
	void add_new_constraint(int len, int ind[], double val[]);
	void setMin_Or_Max(int Min_Or_Max);
	/**
	 * Sets the coefficient of the objective function
	 */
	void set_obj_coeff(unsigned int j, double val);
	/*
	 * Returns the coefficient of the column variable with index j
	 */
	double get_obj_coeff(unsigned int j);
	/**
	 * Solve the LP w.r.t the objective function that is set already
	 */
	double solve();

	/*
	 * Gets the optimal value of the variable, after solving the lp, having the index passed as parameter
	 */
	double get_col_val(int i);

	int getMin_Or_Max();
	void setIteration_Limit(int limits);
	void setInitial_SimplexControlParameters();
	unsigned int getStatus();

	/*double Compute_LLP(std::vector<std::vector<double> > coeff_function, std::vector<std::vector <double> > coeff_constraints,
	 std::vector <double> bounds,std::vector <int> bound_signs);
	 double Compute_LLP(std::vector<std::vector<double> > coeff_function);
	 void setConstraints(std::vector<std::vector <double> > coeff_constraints, std::vector <double> bounds,
	 std::vector <int> bound_signs);
	 *
	 * Executes the simplex method with the given function.
	 */
	double Compute_LLP(std::vector<double> coeff_function) ;

	//returns the vector/point that gives the maximum value
	std::vector<double> getMaximizing_Variables();

	/*
	 * Executes the simplex method with the objective function set to zero. That is assigning all coefficient to zero
	 * This results in checking if the specified constraints has feasible solution or solution is infeasible or has no feasible solution
	 * This application can be used to test the intersection of polytope with guard or with another polytope.
	 */
	unsigned int TestConstraints();

	// ******* Functions to be Interfaced later to the common LP_Solver class **********
	void displayMaxVariables();
	/*
	 * Deletes the GLPK problem
	 */
	void delete_prob();

	void free_glpk_lp_solver();
	static void free_environment_glpk_lp_solver();
	static unsigned int lp_count;
	//   ******************************************************************************

	glpk_lp_solver();
	~glpk_lp_solver();
private:

	glp_prob *mylp;
	glp_smcp param; // simplex parameter structure
	//char *lp_name;
	int Min_Or_Max;		//1 for Min and 2 for Max
	/* optimization direction flag:
	 #define GLP_MIN            1  // minimization
	 #define GLP_MAX            2  // maximization
	 */

	int dimension;				//dimension of the system
	int number_of_constraints;	//rows in glpk or facets of a polytope
	std::vector<double> sv;// values of the variable that maximizes the result
	double result;				//The optimal objective function value

};

#endif /* GLPK_SOLVER_H_ */
