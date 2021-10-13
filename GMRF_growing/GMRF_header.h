#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED



//DEFINITIONS

#include <Eigen/Dense>


extern const int dim_grid; //dimension of the square grid
extern const int dim_prec; //number of columns (and rows) of the precision matrix (nu*nu)
extern const int adj; //number of maximum adjacent elements in a grid
extern int steps; //number of steps of the random walk
extern const int max_steps; //number of max steps of the random walk
extern Eigen::MatrixXd Q; //Precision matrix
extern std::vector < double > X; //GMRF vector
extern std::vector < double > P; //the vector containing the index of the observed elements of the grid
extern std::vector < double > O; //the vector containing the index of the max observed elements of the grid
extern std::vector < std::vector < double > > mat_Z; //matrix of the state of the observations


int GMRF_func(); //contains all functions
int GMRF_model(); //draws from the GMRF
int GMRF_obs(); //draws the observations from a random walk
int GMRF_gold(); //calculates the analytical solutions


#endif
#pragma once
