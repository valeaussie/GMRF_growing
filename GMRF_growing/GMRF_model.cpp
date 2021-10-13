#define _CRT_SECURE_NO_WARNINGS

#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "GMRF_header.h"
#include <cmath>
#include <Eigen/Dense>




const int dim_grid = 4; //dimension of the square grid
const int dim_prec = pow(dim_grid, 2); //dimension of the precision matrix
Eigen::MatrixXd Q; //the precision matrix
std::vector < double > X; //the GMRF vector



Eigen::MatrixXd Precision(int dim_grid);
std::vector < double > Chol_and_LTsol( int dim_grid, Eigen::MatrixXd Prec_Q );


int GMRF_model() {

	GMRF_func();//calls the file with all the functions
	
	//create the precision matrix
	Q = Precision(dim_grid);

	//Cholesky decomposition Q = LL^T ad solution of L^T x = z 
	X = Chol_and_LTsol(dim_grid, Q);


	//creates a dat file with the values of X called "GMRF_vector_X.dat"
	std::ofstream outFile("./GMRF_vector_X.dat");
	for (double n : X) {
		outFile << n << std::endl;
	}
	outFile.close();


	
	return 0;
}