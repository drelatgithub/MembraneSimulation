#include"simulation_process.h"

#include<vector>
#include<math.h>
#include<iostream>
#include<fstream>

#include"surface_mesh.h"
#include"surface_free_energy.h"

// Parameteres in Armijo-Goldstein condition
const double c = 0.5; // Inequality relaxation
const double tau = 0.5; // Shrink size of alpha after each iteration

const double eps = 1e-5; // Maximum deviation for each coordinate

int minimization(std::vector<MS::vertex> &vertices);

int MS::simulation_start(std::vector<vertex> &vertices) {
	//double sum_a_1 = 0, sum_a_2 = 0;
	int len = vertices.size();

	for (int i = 0; i < len; i++) {
		vertices[i].update_geo();
		vertices[i].make_initial();
	}

	minimization(vertices);

	//std::cout << sum_a_1/len<<'\t'<<sqrt(sum_a_2/len-sum_a_1*sum_a_1/len/len)<< std::endl;
	
	return 0;
}

int minimization(std::vector<MS::vertex> &vertices) {
	bool finished = false;
	int N = vertices.size(); // Number of vertices
	double H = 0, H_new = 0;
	double *d_H = new double[3 * N];
	double *d_H_new = new double[3 * N];
	double *p = new double[3 * N]; // Search direction
	double alpha;
	double beta;

	// Initializing
	for (int i = 0; i < N; i++) {
		// Get H and d_H
		H += MS::h_all(&vertices[i]);
		for (int j = 0; j < 3; j++) {
			d_H[i * 3 + j] = MS::d_h_all(&vertices[i], j);
			// Initialize search direction
			p[i * 3 + j] = -d_H[i * 3 + j];
		}

		// Store vertices location as the last location
		vertices[i].make_last();
	}

	int k = 0; // Iteration counter. Only for debug use.

	while (!finished) {
		// Find alpha and update coordinates
		alpha = 1.0; // TODO initial alpha should not be too big (less iterations to find alpha) nor too small (significant decrease in H)
		double judge_lhs, judge_rhs, m = 0;
		for (int i = 0; i < 3*N; i++) {
			m += p[i] * d_H[i];
		}
		judge_rhs = c*m;

		int l = 0; // debug counter
		while(true) {
			for (int i = 0; i < N; i++) {
				vertices[i].point->x = vertices[i].point_last->x + alpha*p[i * 3];
				vertices[i].point->y = vertices[i].point_last->y + alpha*p[i * 3 + 1];
				vertices[i].point->z = vertices[i].point_last->z + alpha*p[i * 3 + 2];
			}
			for (int i = 0; i < N; i++) {
				vertices[i].update_geo();
			}
			H_new = 0;
			for (int i = 0; i < N; i++) {
				H_new += MS::h_all(&vertices[i]);
			}
			judge_lhs = (H_new - H) / alpha;

			l++;

			if (judge_lhs <= judge_rhs)break;
			else alpha *= tau;
		}

		// Find d_H_new (d_H is already calculated in finding alpha)
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < 3; j++) {
				d_H_new[i * 3 + j] = MS::d_h_all(&vertices[i], j);
			}
		}

		// Find beta (Fletcher-Reeves)
		double a = 0, b = 0;
		for (int i = 0; i < 3 * N; i++) {
			a += d_H_new[i] * d_H_new[i];
			b += d_H[i] * d_H[i];
		}
		beta = a / b;

		// Renew search direction
		for (int i = 0; i < 3 * N; i++) {
			p[i] = -d_H_new[i] + beta*p[i];
		}

		// Judge when we should exit loop
		finished = true;
		for (int i = 0; i < N; i++) {
			if (finished && (
				vertices[i].point->x - vertices[i].point_last->x > eps || vertices[i].point_last->x - vertices[i].point->x > eps ||
				vertices[i].point->y - vertices[i].point_last->y > eps || vertices[i].point_last->y - vertices[i].point->y > eps ||
				vertices[i].point->z - vertices[i].point_last->y > eps || vertices[i].point_last->z - vertices[i].point->z > eps)) {
				finished = false;
			}
			// Record coordinates as last-time coordinates
			vertices[i].make_last();
		}

		// Renew H and d_H
		H = H_new;
		for (int i = 0; i < 3 * N; i++) {
			d_H[i] = d_H_new[i];
		}

		k++;
		std::cout << "Iteration: " << k << "\tFree energy: " << H << std::endl;
	}

	delete[]d_H;
	delete[]d_H_new;
	delete[]p;

	return 0;
}