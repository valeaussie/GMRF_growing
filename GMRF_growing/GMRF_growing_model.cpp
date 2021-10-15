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


/********FUNCTIONS*********/
void F_print_matrix(std::vector < std::vector < int > > m);
void F_print_matrix(std::vector < std::vector < double > > M);
void F_print_vector(std::vector < double > v);
void F_print_vector(std::vector < int > v);
void remove_intersection(std::vector<double>& a, std::vector<double>& b);
void remove_duplicates(std::vector<double> &v);

Eigen::MatrixXd Precision(int dim_grid);
std::vector < double > Chol_and_LTsol(int dim_grid, Eigen::MatrixXd Prec_Q);
Eigen::VectorXd Chol_and_LTsol_eigen(int dim_grid, Eigen::MatrixXd Prec_Q);
std::vector < double > Chol_and_Lsol(int dim_grid, Eigen::MatrixXd Prec_Q);
Eigen::VectorXd Chol_and_Lsol_eigen(int dim_grid, Eigen::MatrixXd Prec_Q);
Eigen::MatrixXd Q_AA(std::vector < double > A, Eigen::MatrixXd Prec_Q);
Eigen::MatrixXd Q_AB(std::vector < double > A, Eigen::MatrixXd Prec_Q);
Eigen::VectorXd XAcondXB(Eigen::MatrixXd QA, Eigen::MatrixXd QB, Eigen::VectorXd XB);



int GMRF_model_growing() {

	//for an evolving system:
			//add the observe element at this iteration
			//to the vector of observed elements
	neigh_obs.push_back(o);
	//generate vector of indexes for the neighbours of the newly observed element
	std::vector < double > neigh_unobs;
	for (int qi = 0; qi < sqrt(Q.size()); qi++) {
		if (Q(o, qi) == -1) {
			neigh_unobs.push_back(qi);
		}
	}
	//remove already observed or symulated points
	std::vector < double > toremove;
	toremove.push_back(o);
	for (int k = 0; k < neigh_obs.size(); k++) {
		toremove.push_back(neigh_obs[k]);
	}
	remove_intersection(neigh_unobs, toremove);
	//create a vector with all non zero elements (all the observed or simulated one)
	std::vector < double > all_non_zero;
	all_non_zero.insert(all_non_zero.end(), neigh_obs.begin(), neigh_obs.end());
	all_non_zero.insert(all_non_zero.end(), neigh_unobs.begin(), neigh_unobs.end());
	sort(all_non_zero.begin(), all_non_zero.end());
	remove_duplicates(all_non_zero);
	//create a vector with the idenxes of the observed (in the new matrix)
	std::vector < double > new_neigh_obs;
	for (double n = 0; n < all_non_zero.size(); n++) {
		if (std::find(neigh_obs.begin(), neigh_obs.end(), all_non_zero[n]) != neigh_obs.end()) {
			new_neigh_obs.push_back(n);
		}
	}
	//create a vector with the idenxes of the unobserved (in the new matrix)
	std::vector < double > new_neigh_unobs;
	for (double n = 0; n < all_non_zero.size(); n++) {
		if (std::find(neigh_unobs.begin(), neigh_unobs.end(), all_non_zero[n]) != neigh_unobs.end()) {
			new_neigh_unobs.push_back(n);
		}
	}

	std::cout << "all non zero" << std::endl;
	F_print_vector(all_non_zero);
	std::cout << "neigh obs" << std::endl;
	F_print_vector(neigh_obs);
	std::cout << "neigh unobs" << std::endl;
	F_print_vector(neigh_unobs);
	std::cout << "new neigh obs" << std::endl;
	F_print_vector(new_neigh_obs);
	std::cout << "new neigh unobs" << std::endl;
	F_print_vector(new_neigh_unobs);

	std::vector < double > all_elements;
	for (int i = 0; i < sqrt(Q.size()); i++) {
		all_elements.push_back(i);
	}
	//creates a vector a with the indexes of the elements yet to be observed
	std::vector < double > zero_elements;
	std::remove_copy_if(all_elements.begin(), all_elements.end(), std::back_inserter(zero_elements),
		[&all_non_zero](const int& arg)
	{ return (std::find(all_non_zero.begin(), all_non_zero.end(), arg) != all_non_zero.end()); });
	std::cout << "zero_elements" << std::endl;
	F_print_vector(zero_elements);
	//create the matrices
	Eigen::MatrixXd Qreduced = Q_AA(zero_elements, Q);
	Eigen::MatrixXd Qaa = Q_AA(new_neigh_unobs, Qreduced);
	Eigen::MatrixXd Qab = Q_AB(new_neigh_unobs, Qreduced);
	//create the vector with the observed elements
	std::vector < double > XB;
	for (int i : new_neigh_unobs) {
		XB.push_back(X[i]);
	}
	//typecast new_neigh_obs into eigen
	int new_neigh_unobs_size = new_neigh_unobs.size();
	double* ptr7 = &new_neigh_unobs[0];
	Eigen::Map<Eigen::VectorXd> new_neigh_unobs_eigen(ptr7, new_neigh_unobs_size);
	//sample from conditional
	Eigen::VectorXd samp_x = XAcondXB(Qaa, Qab, new_neigh_unobs_eigen);

	std::cout << "Qreduced" << std::endl;
	std::cout << Qreduced << std::endl;
	std::cout << "Qaa" << std::endl;
	std::cout << Qaa << std::endl;
	std::cout << "Qab" << std::endl;
	std::cout << Qab << std::endl;
	std::cout << "XB" << std::endl;
	F_print_vector(XB);
	std::cout << "samp_x" << std::endl;
	std::cout << samp_x << std::endl;

	//populate neigh_obs for the next step
	for (int k = 0; k < neigh_unobs.size(); k++) {
		neigh_obs.push_back(neigh_unobs[k]);
	}
	remove_duplicates(neigh_obs);

}