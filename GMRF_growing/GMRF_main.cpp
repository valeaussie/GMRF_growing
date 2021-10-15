#include <random>
#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include <algorithm>
#include <unordered_set>
#include <math.h>
#include "GMRF_header.h"
#include "numeric"
#include <Eigen/Dense>
#include <time.h>



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



int main() {

	clock_t tStart = clock();

	GMRF_func(); //calls the file with all the functions
	GMRF_model();  //simulates the Gaussian Markov Random Field
	GMRF_obs(); //simulate the state of the observations

	
	int n = 100; //number of particles


	//define the container for the sampled particles 
	std::vector < std::vector < double > >  sampled(dim_prec, std::vector < double >(n, 0.0));
	//define the container for the resampled particles 
	std::vector < std::vector < std::vector < double > > > resampled(steps, 
		std::vector < std::vector < double > >(dim_prec, std::vector < double >(n, 0.0)));
	//define and initialise the containter for the unnormalised weights
	std::vector < std::vector < double > > un_weights(n, std::vector < double >(steps, 0.0));
	//define and initialise the container for the normalised weights
	std::vector < std::vector < double > > weights(n, std::vector < double >(steps, 0.0));




	std::random_device rd;
	std::mt19937 generator(rd());


	
	//draw the process for n particles for not evolving systems
	Eigen::MatrixXd mat_sim(n, dim_prec);
	for (int i = 0; i < n; i++) {
		Eigen::VectorXd sim = Chol_and_LTsol_eigen(dim_grid, Q);
		Eigen::RowVectorXd row_sim = sim.transpose();
		mat_sim.row(i) = row_sim;
	}
	for(int i = 0; i < n; i++) {
		for (int j = 0; j < dim_prec; j++) {
			sampled[j][i] = mat_sim(i, j);
		}
	}
	
	
	//initialise for evolving system:
	for (int i = 0; i < n; i++) {
		un_weights[i][0] = 1; //set equal weights at time 1
		weights[i][0] = 1.0 / n; //normalise the weights
		resampled[0][0][i] = X[0];
	}

	//typecast X into eigen
	int X_size = X.size();
	double* ptr3 = &X[0];
	Eigen::Map<Eigen::VectorXd> X_eigen(ptr3, X_size);

	
	
	//Iterate:

	std::vector < double > sim_and_obs{ 0 };
	std::vector < double > neigh_obs{ 0 };
	for (int j = 1; j < steps; j++) {

		//populate the resample vector from previous time
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < dim_prec; k++) {
				resampled[j][k][i] = resampled[j - 1][k][i];
			}
		}


		int o = P[j]; // index of observed elements at time j

		//for an evolving system simulate next observed elements:

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
		remove_intersection(neigh_unobs, sim_and_obs);

		//populate the vector with all observed or simulated elements
		sim_and_obs.insert(sim_and_obs.end(), neigh_obs.begin(), neigh_obs.end());
		sim_and_obs.insert(sim_and_obs.end(), neigh_unobs.begin(), neigh_unobs.end());
		sort(sim_and_obs.begin(), sim_and_obs.end());
		remove_duplicates(sim_and_obs);

		//create a vector with the idenxes of the observed (in the new matrix)
		std::vector < double > new_neigh_obs;
		for (double n = 0; n < sim_and_obs.size(); n++) {
			if (std::find(neigh_obs.begin(), neigh_obs.end(), sim_and_obs[n]) != neigh_obs.end()) {
				new_neigh_obs.push_back(n);
			}
		}
		//create a vector with the idenxes of the unobserved (in the new matrix)
		std::vector < double > new_neigh_unobs;
		for (double n = 0; n < sim_and_obs.size(); n++) {
			if (std::find(neigh_unobs.begin(), neigh_unobs.end(), sim_and_obs[n]) != neigh_unobs.end()) {
				new_neigh_unobs.push_back(n);
			}
		}
		
		std::cout << "sim and obs" << std::endl;
		F_print_vector(sim_and_obs);
		std::cout << "neigh obs" << std::endl;
		F_print_vector(neigh_obs);
		std::cout << "neigh unobs" << std::endl;
		F_print_vector(neigh_unobs);
		std::cout << "new neigh obs" << std::endl;
		F_print_vector(new_neigh_obs);
		std::cout << "new neigh unobs" << std::endl;
		F_print_vector(new_neigh_unobs);
		
		if (!neigh_unobs.empty()) {
			std::vector < double > all_elements;
			for (int i = 0; i < sqrt(Q.size()); i++) {
				all_elements.push_back(i);
			}
			//creates a vector a with the indexes of the elements yet to be observed
			std::vector < double > zero_elements;
			std::remove_copy_if(all_elements.begin(), all_elements.end(), std::back_inserter(zero_elements),
				[&sim_and_obs](const int& arg)
			{ return (std::find(sim_and_obs.begin(), sim_and_obs.end(), arg) != sim_and_obs.end()); });
			std::cout << "zero_elements" << std::endl;
			F_print_vector(zero_elements);

			//create the matrices

			Eigen::MatrixXd Qreduced = Q_AA(zero_elements, Q);
			std::vector < double > count_of_reduced;
			for (double n = 0; n < sqrt(Qreduced.size()); n++) {
				count_of_reduced.push_back(n);
			}
			remove_intersection(count_of_reduced, new_neigh_unobs);
			std::cout << "count of reduced" << std::endl;
			F_print_vector(count_of_reduced);

			Eigen::MatrixXd Qaa = Q_AA(count_of_reduced, Qreduced);

			Eigen::MatrixXd Qab = Q_AB(count_of_reduced, Qreduced);

			//create the vector with the observed elements	
			std::vector < double > XB;
			for (int i : new_neigh_unobs) {
				XB.push_back(X[i]);
			}
			//typecast new_neigh_obs into eigen
			int count_of_reduced_size = count_of_reduced.size();
			double* ptr7 = &count_of_reduced[0];
			Eigen::Map<Eigen::VectorXd> count_of_reduced_eigen(ptr7, count_of_reduced_size);
			//sample from conditional

			Eigen::VectorXd samp_x = XAcondXB(Qaa, Qab, count_of_reduced_eigen);

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
				sim_and_obs.push_back(neigh_unobs[k]);
			}
			remove_duplicates(neigh_obs);

			//populated resampled with newly sampled and observed elements
			for (int i = 0; i < n; i++) {
				for (int k = 0; k < neigh_unobs.size(); k++) {
					resampled[j][neigh_unobs[k]][i] = samp_x(k);
				}
			}
		}

		
		//make correction
		for (int i = 0; i < n; i++) {
			resampled[j][o][i] = X[o];
		}

		//calculate weights
		for (int i = 0; i < n; i++) {
			//fill vector sim
			Eigen::VectorXd sim = Eigen::VectorXd::Zero(dim_prec);
			for (int k = 0; k < dim_prec; k++) {
				sim(k) = resampled[j][k][i]; //populate Eigen::vector sim
			}
			
			//calculate weights for particle i at time j
			double Xosqr = pow(X_eigen(o), 2);
			double Simosqr = pow(sim(o), 2);
			double Xo = X_eigen(o);
			double simo = sim(o);
			double w = exp(
				( ( Q.row(o) * sim - simo * Q(o,o) ) * ( simo - Xo)
				+  0.5 * Q(o,o) * ( Simosqr - Xosqr )  )
			);
			un_weights[i][j] = w; 
		}

		//normalise the weights
		double sum_of_weights{ 0 };
		for (int i = 0; i < n; i++) {
			sum_of_weights = sum_of_weights + un_weights[i][j];
		}

		for (int i = 0; i < n; i++) {
			weights[i][j] = un_weights[i][j] / sum_of_weights;
		}
				
		//resampling
		std::vector < double > drawing_vector(n, 0.0);
		for (int i = 0; i < n; i++) {
			drawing_vector[i] = un_weights[i][j];
		}			
		double index_resampled;
		std::vector < int > vec_index;
		for (int i = 0; i < n; i++) {
			std::discrete_distribution < int > discrete(drawing_vector.begin(), drawing_vector.end());
			index_resampled = discrete(generator);
			vec_index.push_back(index_resampled);
		}
		std::vector < std::vector < double > > newmatrix(dim_prec, std::vector < double >(n, 0.0));
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < dim_prec; k++) {
				newmatrix[k][i] = resampled[j][k][vec_index[i]];
			}
		}
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < dim_prec; k++) {
				resampled[j][k][i] = newmatrix[k][i];
			}
		}
		/*
		//add jitter
		std::normal_distribution< double > nor(0,0.1);
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < dim_prec; k++) {
				resampled[j][k][i] = resampled[j][k][i] + nor(generator);
			}
		}*/
	}

	
	std::cout << "P" << std::endl;
	F_print_vector(P);

	
	std::ofstream outFile("./resampled_000.csv");
	outFile << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < dim_prec; j++) {
			outFile << resampled[steps-1][j][i] << ",";
		}
		outFile << std::endl;
	}
	outFile.close();

	std::ofstream outFile0("./resampled_001.csv");
	outFile0 << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < dim_prec; j++) {
			outFile0 << resampled[1][j][i] << ",";
		}
		outFile0 << std::endl;
	}
	outFile0.close();

	std::ofstream outFile1("./sampled.csv");
	outFile1 << std::endl;
	for (int i = 0; i < dim_prec; i++) {
		for (int j = 0; j < n; j++) {
			outFile1 << sampled[i][j] << ",";
		}
		outFile1 << std::endl;
	}
	outFile1.close();


	std::ofstream outFile2("./weights.csv");
	outFile2 << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < steps; j++) {
			outFile2 << weights[i][j] << ",";
		}
		outFile2 << std::endl;
	}
	outFile2.close();

	std::ofstream outFile3("./un_weights.csv");
	outFile3 << std::endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < steps; j++) {
			outFile3 << un_weights[i][j] << ",";
		}
		outFile3 << std::endl;
	}
	outFile3.close();
	
	
	//create a .csv file containing the parameters
	std::ofstream outparam("./parameters.csv");
	outparam << "mu" << "," << "max_steps" << "," << "part" << std::endl;
	outparam << dim_prec << "," << max_steps << "," << n << "," << std::endl;
	outparam.close();


	GMRF_gold();

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	return 0;

}

