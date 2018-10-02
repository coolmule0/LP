#include <ilcplex/ilocplex.h>
#include "../LP/FileIO.h"
#include <glm/glm.hpp>

ILOSTLBEGIN



int main(int argc, const char* argv[]) {

	int batches = 0; //number of LPs
	int size = 0; //size of each LP
	glm::vec2* x = NULL;  //< solution to randomly generated LP
	std::vector<glm::vec2> A; //array of constraints for 1 lp
	std::vector<float> b; //array of constraints for 1 lp
	glm::vec2 optimiseSingle; // variables to minimise in optimisation function for 1 lp
	glm::vec4* constraints = NULL; // array of constraints

	glm::vec2* optimise = NULL; //variable to minimise in optimisation function
	glm::vec2* output = NULL; // array of optimal results in x & y.

	IloEnv env;							  

	//------------------------------------------
	//handle args
	if (argc != 3) {
		printf("\nIncorrect Number of Arguments!\n");
		printf("Correct Usages/Syntax:\n");
		printf("Argument 1) File-Name -- Name of input file in benchmarks folder. Cannot contain spaces\n");
		printf("Argument 2) Batch-size -- Number of LPs to be solved\n");
		return 1;
	}
	batches = atoi(argv[2]);

	//------------------------------------------
	//handle input

	//Input is from file
	printf("Parsing input files... ");
	if (!parseBenchmark3(argv[1], A, b, &optimiseSingle, &size)) {
		return 1;
	}
	printf("Done\n");

	//------------------------------------------
	//set up problem
	try {
		IloModel model(env);
		IloNumVarArray vars(env);
		vars.add(IloNumVar(env));
		vars.add(IloNumVar(env));
		model.add(IloMaximize(env, optimiseSingle.x * vars[0] + optimiseSingle.y * vars[1] ));
		//model.add(-vars[0] + vars[1]  <= 20);
		//model.add(vars[0] - 3 * vars[1]  <= 30);

		//IloModel model(env);
		//IloNumVarArray vars(env);
		//vars.add(IloNumVarArray(env)); // 0<x<inf
		//vars.add(IloNumVarArray(env)); // 0<y<inf

		////model.add(IloMinimize(env, optimiseSingle.x * vars[0] + optimiseSingle.y * vars[1]));	// Objective Function: Minimize distance

		////add constraints
		for (int i = 0; i < size; i++) {
			model.add(A[i].x * vars[0] + A[i].y * vars[1] <= b[i]);
		}

		//run problem
		IloCplex cplex(model);
		if (!cplex.solve()) {
			env.error() << "Failed to optimize LP." << endl;
			throw(-1);
		}

		IloNumArray vals(env);
		env.out() << "Solution status = " << cplex.getStatus() << endl;
		env.out() << "Solution value = " << cplex.getObjValue() << endl;
		cplex.getValues(vals, vars);
		env.out() << "Values = " << vals << endl;

	}
	catch (IloException& e) {
		cerr << "Concert exception caught: " << e << endl;
	}
	catch (...) {
		cerr << "Unknown exception caught" << endl;
	}

	env.end();

	return 0;
}


// -------------------------------------------------------------- -*- C++ -*-
// File: blend.cpp
// Version 12.7.1  
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
// 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
// Copyright IBM Corporation 2000, 2017. All Rights Reserved.
//
// US Government Users Restricted Rights - Use, duplication or
// disclosure restricted by GSA ADP Schedule Contract with
// IBM Corp.
// --------------------------------------------------------------------------
//
// blend.cpp -- A blending problem

/* ------------------------------------------------------------

Problem Description
-------------------

Goal is to blend four sources to produce an alloy: pure metal, raw
materials, scrap, and ingots.

Each source has a cost.
Each source is made up of elements in different proportions.
Alloy has minimum and maximum proportion of each element.

Minimize cost of producing a requested quantity of alloy.

------------------------------------------------------------ */


//IloInt nbElements, nbRaw, nbScrap, nbIngot;
//IloNum alloy;
//IloNumArray nm, nr, ns, ni, p, P;
//
//IloNumArray2 PRaw, PScrap, PIngot;
//
//
//void define_data(IloEnv env) {
//	nbElements = 3;
//	nbRaw = 2;
//	nbScrap = 2;
//	nbIngot = 1;
//	alloy = 71;
//	nm = IloNumArray(env, nbElements, 22.0, 10.0, 13.0);
//	nr = IloNumArray(env, nbRaw, 6.0, 5.0);
//	ns = IloNumArray(env, nbScrap, 7.0, 8.0);
//	ni = IloNumArray(env, nbIngot, 9.0);
//	p = IloNumArray(env, nbElements, 0.05, 0.30, 0.60);
//	P = IloNumArray(env, nbElements, 0.10, 0.40, 0.80);
//	PRaw = IloNumArray2(env, nbElements);
//	PScrap = IloNumArray2(env, nbElements);
//	PIngot = IloNumArray2(env, nbElements);
//	PRaw[0] = IloNumArray(env, nbRaw, 0.20, 0.01);
//	PRaw[1] = IloNumArray(env, nbRaw, 0.05, 0.00);
//	PRaw[2] = IloNumArray(env, nbRaw, 0.05, 0.30);
//
//	PScrap[0] = IloNumArray(env, nbScrap, 0.00, 0.01);
//	PScrap[1] = IloNumArray(env, nbScrap, 0.60, 0.00);
//	PScrap[2] = IloNumArray(env, nbScrap, 0.40, 0.70);
//
//	PIngot[0] = IloNumArray(env, nbIngot, 0.10);
//	PIngot[1] = IloNumArray(env, nbIngot, 0.45);
//	PIngot[2] = IloNumArray(env, nbIngot, 0.45);
//}
//
//int
//main(int, char**)
//{
//	IloEnv env;
//	try {
//		IloInt j;
//
//		define_data(env);
//
//		IloModel model(env);
//
//		IloNumVarArray m(env, nbElements, 0.0, IloInfinity);
//		IloNumVarArray r(env, nbRaw, 0.0, IloInfinity);
//		IloNumVarArray s(env, nbScrap, 0.0, IloInfinity);
//		IloNumVarArray i(env, nbIngot, 0.0, 100000);
//		IloNumVarArray e(env, nbElements);
//
//		// Objective Function: Minimize Cost
//		model.add(IloMinimize(env, IloScalProd(nm, m) + IloScalProd(nr, r) +
//			IloScalProd(ns, s) + IloScalProd(ni, i)));
//
//		// Min and max quantity of each element in alloy
//		for (j = 0; j < nbElements; j++) {
//			e[j] = IloNumVar(env, p[j] * alloy, P[j] * alloy);
//		}
//
//		// Constraint: produce requested quantity of alloy
//		model.add(IloSum(e) == alloy);
//
//		// Constraints: Satisfy element quantity requirements for alloy
//		for (j = 0; j < nbElements; j++) {
//			model.add(e[j] == m[j] + IloScalProd(PRaw[j], r)
//				+ IloScalProd(PScrap[j], s)
//				+ IloScalProd(PIngot[j], i));
//		}
//
//		// Optimize
//		IloCplex cplex(model);
//		cplex.setOut(env.getNullStream());
//		cplex.setWarning(env.getNullStream());
//		cplex.solve();
//
//		if (cplex.getStatus() == IloAlgorithm::Infeasible)
//			env.out() << "No Solution" << endl;
//
//		env.out() << "Solution status: " << cplex.getStatus() << endl;
//
//		// Print results
//		env.out() << "Cost:" << cplex.getObjValue() << endl;
//		env.out() << "Pure metal:" << endl;
//		for (j = 0; j < nbElements; j++)
//			env.out() << j << ") " << cplex.getValue(m[j]) << endl;
//		env.out() << "Raw material:" << endl;
//		for (j = 0; j < nbRaw; j++)
//			env.out() << j << ") " << cplex.getValue(r[j]) << endl;
//		env.out() << "Scrap:" << endl;
//		for (j = 0; j < nbScrap; j++)
//			env.out() << j << ") " << cplex.getValue(s[j]) << endl;
//		env.out() << "Ingots : " << endl;
//		for (j = 0; j < nbIngot; j++)
//			env.out() << j << ") " << cplex.getValue(i[j]) << endl;
//		env.out() << "Elements:" << endl;
//		for (j = 0; j < nbElements; j++)
//			env.out() << j << ") " << cplex.getValue(e[j]) << endl;
//	}
//	catch (IloException& ex) {
//		cerr << "Error: " << ex << endl;
//	}
//	catch (...) {
//		cerr << "Error" << endl;
//	}
//	env.end();
//	return 0;
//}
//
///*
//Cost:653.554
//Pure metal:
//0) 0
//1) 0
//2) 0
//Raw material:
//0) 0
//1) 0
//Scrap:
//0) 17.059
//1) 30.2311
//Ingots :
//0) 32.4769
//Elements:
//0) 3.55
//1) 24.85
//2) 42.6
//*/
