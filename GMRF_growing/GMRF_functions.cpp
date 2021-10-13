#define _CRT_SECURE_NO_WARNINGS

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


using namespace std;


int GMRF_func() {

	return 0;
}


//*********** GENERIC FUNCTIONS ***************

//Prints a matrix of int
void F_print_matrix(vector < vector < int > > m) {
	for (const vector < int > v : m) {
		for (int x : v) cout << x << ' ';
		cout << endl;
	}
}
//Prints a matrix of signed doubles
void F_print_matrix(vector < vector < double > > M) {
	for (const vector < double > v : M) {
		for (double x : v) cout << x << ' ';
		cout << endl;
	}
}
//Prints a vector of doubles
void F_print_vector(vector < double > v) {
	for (const double x : v) cout << x << ' ';
	cout << endl;
}

//Prints a vector of doubles
void F_print_vector(vector < int > v) {
	for (const size_t x : v) cout << x << ' ';
	cout << endl;
}

//Remove duplicates from a vector
void remove_duplicates(std::vector<double> &v)
{
	auto end = v.end();
	for (auto it = v.begin(); it != end; ++it) {
		end = std::remove(it + 1, end, *it);
	}

	v.erase(end, v.end());
}

//Remove elements of one vector from another vector
void remove_intersection(std::vector<double>& a, std::vector<double>& b) {
	std::unordered_multiset<double> st;
	st.insert(a.begin(), a.end());
	st.insert(b.begin(), b.end());
	auto predicate = [&st](const double& k) { return st.count(k) > 1; };
	a.erase(std::remove_if(a.begin(), a.end(), predicate), a.end());
}



//creates the precision matrix (sparse) for any value of the dimension of the grid
Eigen::MatrixXd Precision(int dim_grid) {

	float param_delta = 0.01;

	int dim_prec = dim_grid * dim_grid;
	Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(dim_prec, dim_prec);
	for (int i = 0; i < dim_prec; i++) {
		for (int j = 0; j < dim_prec; j++) {
			Q(i, i) = 4 + param_delta;
		}
	}
	// central elements of the grid
	for (int i = 0; i < dim_prec; i++) {
		for (int j = 0; j < dim_prec; j++) {
			for (int k = 1; k < dim_grid - 1; k++) {
				for (int i = (dim_grid*k) + 1; i < (k + 1) * dim_grid - 1; i++) {
					for (int j = 0; j < dim_prec; j++) {
						Q(i, i - 1) = -1;
						Q(i, i + 1) = -1;
						Q(i, i - dim_grid) = -1;
						Q(i, i + dim_grid) = -1;
					}
				}
			}
		}
	}
	// corners of the grid
	for (int i = 0; i < dim_prec; i++) {
		for (int j = 0; j < dim_prec; j++) {
			// top left corner
			if (i == 0) {
				Q(i, i + 1) = -1;
				Q(i, i + dim_grid) = -1;
				Q(i, i + dim_grid - 1) = -1; // circular
				Q(i, dim_prec - dim_grid) = -1; // circular
			}
			// top right corner
			if (i == dim_grid - 1) {
				Q(i, i - 1) = -1;
				Q(i, i + dim_grid) = -1;
				Q(i, dim_grid - 1 - i) = -1; // circular
				Q(i, dim_prec - 1) = -1; // circular
			}
			// bottom right corner
			if (i == dim_prec - 1) {
				Q(i, i - 1) = -1;
				Q(i, i - dim_grid) = -1;
				Q(i, i - dim_grid + 1) = -1; // circular
				Q(i, dim_grid - 1) = -1; // circular
			}
			//bottom left corner
			if (i == dim_prec - dim_grid) {
				Q(i, i + 1) = -1;
				Q(i, i - dim_grid) = -1;
				Q(i, i + dim_grid - 1) = -1; // circular
				Q(i, 0) = -1; // circular
			}
			for (int k = 1; k < dim_grid - 1; k++) {
				// top side
				if (i == k) {
					Q(i, i - 1) = -1;
					Q(i, i + 1) = -1;
					Q(i, i + dim_grid) = -1;
					Q(i, (dim_prec - dim_grid) + i) = -1; // circular
				}
				// right side
				if (i == (k + 1) * dim_grid - 1) {
					Q(i, i - 1) = -1;
					Q(i, i - dim_grid) = -1;
					Q(i, i - dim_grid + 1) = -1; // circular
					Q(i, i + dim_grid) = -1;
				}
				// bottom side
				if (i == dim_prec - k - 1) {
					Q(i, i - 1) = -1;
					Q(i, i + 1) = -1;
					Q(i, i - dim_grid) = -1;
					Q(i, dim_grid - k - 1) = -1; // circular
				}
				// left side
				if (i == dim_prec - (k + 1) * dim_grid) {
					Q(i, i + 1) = -1;
					Q(i, i - dim_grid) = -1;
					Q(i, i + dim_grid - 1) = -1; // circular
					Q(i, i + dim_grid) = -1;
				}
			}
		}
	}
	return Q;
}




//Cholesky decomposition Q = LL^T ad solution of L^T x = z. Outputs an Std::vector
std::vector < double > Chol_and_LTsol(int DimGrid, Eigen::MatrixXd Prec_Q) {

	int dim_prec = DimGrid * DimGrid;
	std::vector < double > GMRF_vec(dim_prec, 0);

	//Sampling z from a standard multivariate normal distribution of size mu
	//Eigen::VectorXd z = SampleMND(dim_prec);
	
	std::random_device rd;
	std::mt19937 generator(rd());

	Eigen::VectorXd z = Eigen::VectorXd::Zero(dim_prec);
	for (int ind = 0; ind < dim_prec; ind++) {
		std::normal_distribution < double > normal{0, 1};
		z(ind) = normal(generator);
	}


	Eigen::LLT<Eigen::MatrixXd> lltOfQ(Prec_Q);
	Eigen::MatrixXd U = lltOfQ.matrixU();
	Eigen::VectorXd x = U.colPivHouseholderQr().solve(z);

	//creating an std::vector from Eigen::VectorXd
	Eigen::Map<Eigen::VectorXd>(GMRF_vec.data(), dim_prec) = x;


	return GMRF_vec;
}

//Cholesky decomposition Q = LL^T ad solution of L^T x = z. Outputs an Eigen::vector
Eigen::VectorXd Chol_and_LTsol_eigen(int DimGrid, Eigen::MatrixXd Prec_Q) {

	int dim_prec = DimGrid * DimGrid;

	//Sampling z from a standard multivariate normal distribution of size mu
	//Eigen::VectorXd z = SampleMND(dim_prec);

	std::random_device rd;
	std::mt19937 generator(rd());

	Eigen::VectorXd z = Eigen::VectorXd::Zero(dim_prec);
	for (int ind = 0; ind < dim_prec; ind++) {
		std::normal_distribution < double > normal{ 0, 1 };
		z(ind) = normal(generator);
	}

	Eigen::LLT<Eigen::MatrixXd> lltOfQ(Prec_Q);
	Eigen::MatrixXd U = lltOfQ.matrixU();
	Eigen::VectorXd x = U.colPivHouseholderQr().solve(z);

	return x;
}



//Cholesky decomposition Q = LL^T ad solution of L x = z. Outputs an Std::vector
std::vector < double > Chol_and_Lsol( int DimGrid, Eigen::MatrixXd Prec_Q ) {

	int dim_prec = DimGrid * DimGrid;
	std::vector < double > GMRF_vec( dim_prec, 0.0 );

	//Sampling z from a standard multivariate normal distribution of size mu
	//Eigen::VectorXd z = SampleMND(dim_prec);

	std::random_device rd;
	std::mt19937 generator(rd());

	Eigen::VectorXd z = Eigen::VectorXd::Zero(dim_prec);
	for (int ind = 0; ind < dim_prec; ind++) {
		std::normal_distribution < double > normal{ 0, 1 };
		z(ind) = normal(generator);
	}

	Eigen::LLT<Eigen::MatrixXd> lltOfQ(Prec_Q);
	Eigen::MatrixXd L = lltOfQ.matrixL();
	Eigen::VectorXd x = L.colPivHouseholderQr().solve(z);

	//creating an std::vector from Eigen::VectorXd
	Eigen::Map<Eigen::VectorXd>(GMRF_vec.data(), dim_prec) = x;

	return GMRF_vec;
}


//Cholesky decomposition Q = LL^T ad solution of L x = z. Outputs an Eigen::vector
Eigen::VectorXd Chol_and_Lsol_eigen(int DimGrid, Eigen::MatrixXd Prec_Q) {

	int dim_prec = DimGrid * DimGrid;

	//Sampling z from a standard multivariate normal distribution of size mu
	//Eigen::VectorXd z = SampleMND(dim_prec);

	std::random_device rd;
	std::mt19937 generator(rd());

	Eigen::VectorXd z = Eigen::VectorXd::Zero(dim_prec);
	for (int ind = 0; ind < dim_prec; ind++) {
		std::normal_distribution < double > normal{ 0, 1 };
		z(ind) = normal(generator);
	}

	Eigen::LLT<Eigen::MatrixXd> lltOfQ(Prec_Q);
	Eigen::MatrixXd L = lltOfQ.matrixL();
	Eigen::VectorXd x = L.colPivHouseholderQr().solve(z);

	return x;
}



//creates a patition Q_AA of the precision matrix Q 
//for the unobserved elements. Vector b contains the indexes of the observed elements.
Eigen::MatrixXd Q_AA(std::vector < double > b, Eigen::MatrixXd Prec_Q) {
	
	int dim_b = b.size();
	//creates a vector a with the indexes of the unobserved elements	
	int dim_a = sqrt(Prec_Q.size()) - dim_b;
	std::vector < double > v;
	for (int i = 0; i < sqrt(Prec_Q.size()); i++) {
		v.push_back(i);
	}
	//creates a vector a with the elements of v that are not in b
	std::vector < double > a;
	std::remove_copy_if(v.begin(), v.end(), std::back_inserter(a),
		[&b](const int& arg)
	{ return (std::find(b.begin(), b.end(), arg) != b.end()); });

	Eigen::MatrixXd QAA = Eigen::MatrixXd::Zero(dim_a, dim_a);

	for (int i = 0; i < dim_a; i++) {
		for (int j = 0; j < dim_a; j++) {
			int l = a[i];
			int m = a[j];
			QAA(i, j) = Prec_Q(l, m);
		}
	}
	return QAA;
}




//creates a patition Q_AB of the precision matrix Q
//Vector b contains the indexes of the observed elements.
Eigen::MatrixXd Q_AB(std::vector < double > b, Eigen::MatrixXd Prec_Q) {

	int dim_b = b.size();
	//creates a vector a with the indexes of the unobserved elements	
	int dim_a = sqrt(Prec_Q.size()) - dim_b;
	std::vector < double > v;
	for (int i = 0; i < sqrt(Prec_Q.size()); i++) {
		v.push_back(i);
	}
	//creates a vector a with the elements of v that are not in b
	std::vector < double > a;
	std::remove_copy_if(v.begin(), v.end(), std::back_inserter(a),
		[&b](const int& arg)
	{ return (std::find(b.begin(), b.end(), arg) != b.end()); });
	//creates the partition Q_AB
	Eigen::MatrixXd QAB = Eigen::MatrixXd::Zero(dim_a, dim_b);
	for (int i = 0; i < dim_a; i++) {
		for (int j = 0; j < dim_b; j++) {
			int l = b[j];
			int m = a[i];
			QAB(i, j) = Prec_Q(m, l);
		}
	}
	return QAB;
}

//sample from conditional x_A | x_B
Eigen::VectorXd XAcondXB (Eigen::MatrixXd QA, Eigen::MatrixXd QB, Eigen::VectorXd XB) {
	std::random_device rd;
	std::mt19937 generator(rd());
	//calcualate mean
	Eigen::VectorXd mean = QB * XB;
	//solve L w = mean
	Eigen::LLT<Eigen::MatrixXd> lltOfQA(QA);
	Eigen::MatrixXd L = lltOfQA.matrixL();
	Eigen::VectorXd w = L.colPivHouseholderQr().solve(mean);
	//solve LT mu = w
	Eigen::MatrixXd U = lltOfQA.matrixU();
	Eigen::VectorXd mu = U.colPivHouseholderQr().solve(w);
	//sample from normal distribution
	int dim_QA = sqrt(QA.size());
	Eigen::VectorXd z = Eigen::VectorXd::Zero(dim_QA);
	for (int ind = 0; ind < dim_QA; ind++) {
		std::normal_distribution < double > normal{ 0, 1 };
		z(ind) = normal(generator);
	}
	//compute LT v = z
	Eigen::VectorXd v = U.colPivHouseholderQr().solve(z);
	//compute x = mu + v
	Eigen::VectorXd samp_x = mu + v;

	return samp_x;
}




