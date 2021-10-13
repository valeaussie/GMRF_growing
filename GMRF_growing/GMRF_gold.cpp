#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <math.h>
#include "GMRF_header.h"
#include "numeric"
#include <Eigen/Dense>



Eigen::MatrixXd Q_AA(std::vector < double > A, Eigen::MatrixXd Prec_Q);
Eigen::MatrixXd Q_AB(std::vector < double > A, Eigen::MatrixXd Prec_Q);
void F_outExp(std::vector <double> exp);

int GMRF_gold() {
	
	//calculate the mean of the conditional posterior

	Eigen::MatrixXd QAB = Q_AB(P, Q);

	//crete a vector of size mu of non negative integers
	std::vector < double > integers(dim_prec); //vectors of size dim_prec
	std::iota(std::begin(integers), std::end(integers), 0); //fill with 0, 1, ..., dim_prec

	int size_of_obs = P.size();

	//calculate the matrix QAA
	Eigen::MatrixXd QAA = Q_AA(P, Q);


	//typecast P into eigen
	std::vector < double > obs;
	for (int i : P) {
		obs.push_back(X[i]);
	}

	int obs_size = obs.size();
	double* ptr4 = &obs[0];
	Eigen::Map<Eigen::VectorXd> obs_eigen(ptr4, obs_size);
	
	
	//calculate the mean
	Eigen::VectorXd mean = -(QAA.inverse() * QAB) * obs_eigen;

	//typecast Eigen to std
	std::vector < double > mean_std(mean.size(), 0.0);
	Eigen::Map<Eigen::VectorXd>(mean_std.data(), mean.size()) = mean;


	//Create a dat files for expectations and variances
	F_outExp(mean_std);

	
	return 0;
}

//Creates a dat file with the values of the Expectations
void F_outExp(std::vector <double> exp) {
	std::ofstream outFile("./GMRF_gold_000.csv");
	outFile << std::endl;
	for (double n : exp) {
		outFile << n << std::endl;
	}
	outFile.close();
}