#include"simulation_process.h"

#include<vector>
#include<math.h>
#include<iostream>
#include<fstream>

#include"surface_mesh.h"
#include"surface_free_energy.h"

// Parameteres in Armijo-Goldstein condition
const double c = 0.4; // Inequality relaxation
const double tau = 0.5; // Shrink size of alpha after each iteration

const double eps = 1e-3; // Maximum deviation for each coordinate

int minimization(std::vector<MS::vertex> &vertices);
void test_derivatives(std::vector<MS::vertex> &vertices);

int MS::simulation_start(std::vector<vertex> &vertices) {
	int len = vertices.size();

	for (int i = 0; i < len; i++) {
		vertices[i].count_neighbours();
		if (i == 0) {
			int a = 1;
		}
		vertices[i].update_geo();
		vertices[i].make_initial();
	}

	std::ofstream p_out;
	p_out.open("F:\\p_out.txt");

	//minimization(vertices);

	
	for (double a = -1; a < 1.5; a += 0.05) {
		std::cout << update_len(a) << std::endl;
		minimization(vertices);
		for (int i = 0; i < len; i++) {
			p_out << vertices[i].point->x << '\t' << vertices[i].point->y << '\t' << vertices[i].point->z << '\t';
		}
		p_out << '\n';
	}
	

	p_out.close();

	//test_derivatives(vertices);
	
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

	MS::update_len(-1);

	//std::ofstream tp_out;
	//tp_out.open("F:\\tp_out.txt");
	
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
		alpha = 0.09; // TODO initial alpha should not be too big (less iterations to find alpha) nor too small (significant decrease in H)
		double judge_lhs, judge_rhs, m = 0;
		for (int i = 0; i < 3*N; i++) {
			m += p[i] * d_H[i];
		}
		if (m > 0)
			std::cout << "Warning: along search direction is increasing free energy\n";
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
		//double a = 0, b = 0;
		//for (int i = 0; i < 3 * N; i++) {
		//	a += d_H_new[i] * d_H_new[i];
		//	b += d_H[i] * d_H[i];
		//}
		//beta = a / b;
		// Find beta (Polak-Ribiere)
		double a = 0, b = 0;
		for (int i = 0; i < 3 * N; i++) {
			a += d_H_new[i] * (d_H_new[i] - d_H[i]);
			b += d_H[i] * d_H[i];
		}
		beta = (a >= 0) ? a / b : 0;

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
				vertices[i].point->z - vertices[i].point_last->z > eps || vertices[i].point_last->z - vertices[i].point->z > eps)) {
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
		std::cout << k << "\tm: " << m << "\talpha: " << alpha << "\tFree energy: " << H << std::endl;
		//for (int i = 0; i < N; i++) {
		//	tp_out << vertices[i].point->x << '\t' << vertices[i].point->y << '\t' << vertices[i].point->z << '\t';
		//}
		//tp_out << '\n';
	}

	delete[]d_H;
	delete[]d_H_new;
	delete[]p;

	return 0;
}

void test_derivatives(std::vector<MS::vertex> &vertices) {
	// Test local properties (main)
	/*
	double l1 = vertices[0].curv_h, l2 = vertices[0].curv_g, l3 = vertices[0].cot_theta2[0];
	double d1 = vertices[0].dz_curv_h, d2 = vertices[0].dz_curv_g, d3 = vertices[0].dz_cot_theta2[0];

	vertices[0].point->z += 0.0001;
	vertices[0].update_geo();

	double c1 = (vertices[0].curv_h - l1) / .0001,
		c2 = (vertices[0].curv_g - l2) / .0001,
		c3 = (vertices[0].cot_theta2[0] - l3) / .0001;
	std::cout << "dx h:\tactual: " << c1 << "\tsupposed: " << d1 << std::endl;
	std::cout << "dx g:\tactual: " << c2 << "\tsupposed: " << d2 << std::endl;
	std::cout << "dx t2:\tactual: " << c3 << "\tsupposed: " << d3 << std::endl;
	*/

	// Test local properties (neighbours)
	/*
	MS::vertex *n = vertices[0].n[0], *nn = vertices[0].n_next[0], *np = vertices[0].n_prev[0];
	//double l1 = vertices[0].curv_h, l2 = vertices[0].curv_g, l3 = vertices[0].area, l4 = vertices[0].cot_theta3[0];
	double l1 = vertices[0].cot_theta3[0], l2 = vertices[0].theta3[0];
	//double d1 = vertices[0].dyn_curv_h[0], d2 = vertices[0].dyn_curv_g[0], d3 = vertices[0].dyn_area[0], d4=vertices[0].dyn_cot_theta3[0];
	double d1 = vertices[0].dznn_cot_theta3[0], d2 = vertices[0].dznn_theta3[0];

	//n->point->y += 0.0001;
	nn->point->z += 0.0001;
	vertices[0].update_geo();

	//double c1 = (vertices[0].curv_h - l1) / .0001,
		//c2 = (vertices[0].curv_g - l2) / .0001,
		//c3 = (vertices[0].area - l3) / .0001,
		//c4 = (vertices[0].cot_theta3[0] - l4) / .0001;
	double c1 = (vertices[0].cot_theta3[0] - l1) / .0001, c2 = (vertices[0].theta3[0] - l2) / .0001;
	std::cout << "dz:\tactual: " << c1 << "\tsupposed: " << d1 << std::endl;
	std::cout << "dz:\tactual: " << c2 << "\tsupposed: " << d2 << std::endl;
	//std::cout << "dy a:\tactual: " << c3 << "\tsupposed: " << d3 << std::endl;
	//std::cout << "dy r:\tactual: " << c4 << "\tsupposed: " << d4 << std::endl;
	*/

	// Test local free energies with neighbours
	/*
	double l1 = MS::h_curv_h(&vertices[0]), l2 = MS::h_curv_g(&vertices[0]), l3 = MS::h_tension(&vertices[0]);
	for (int i = 0; i < vertices[0].neighbours; i++) {
		l1 += MS::h_curv_h(vertices[0].n[i]);
		l2 += MS::h_curv_g(vertices[0].n[i]);
		l3 += MS::h_tension(vertices[0].n[i]);
	}
	double d1 = MS::d_h_curv_h(&vertices[0], 1), d2 = MS::d_h_curv_g(&vertices[0], 1), d3 = MS::d_h_tension(&vertices[0], 1);


	vertices[0].point->y += 0.00001;
	vertices[0].update_geo();
	for (int i = 0; i < vertices[0].neighbours; i++) {
		vertices[0].n[i]->update_geo();
	}

	double n1 = MS::h_curv_h(&vertices[0]), n2 = MS::h_curv_g(&vertices[0]), n3 = MS::h_tension(&vertices[0]);
	for (int i = 0; i < vertices[0].neighbours; i++) {
		n1 += MS::h_curv_h(vertices[0].n[i]);
		n2 += MS::h_curv_g(vertices[0].n[i]);
		n3 += MS::h_tension(vertices[0].n[i]);
	}
	double c1 = (n1 - l1) / .00001,
		c2 = (n2 - l2) / .00001,
		c3 = (n3 - l3) / .00001;

	std::cout << "dy h:\tactual: " << c1 << "\tsupposed: " << d1 << std::endl;
	std::cout << "dy g:\tactual: " << c2 << "\tsupposed: " << d2 << std::endl;
	std::cout << "dy t:\tactual: " << c3 << "\tsupposed: " << d3 << std::endl;
	*/

	// Test Overall free energy
	int N = vertices.size();
	double H = 0, H_new = 0;
	double *d_H = new double[3 * N];
	double *d_H_new = new double[3 * N];

	for (int i = 0; i < N; i++) {
		H += MS::h_all(&vertices[i]);
		for (int j = 0; j < 3; j++) {
			d_H[i * 3 + j] = MS::d_h_all(&vertices[i], j);
		}
	}

	vertices[0].point->y += 0.000001;
	vertices[100].point->x += 0.000001;
	vertices[200].point->z += 0.000001;
	for (int i = 0; i < N; i++) {
		vertices[i].update_geo();
	}

	for (int i = 0; i < N; i++) {
		H_new += MS::h_all(&vertices[i]);
		for (int j = 0; j < 3; j++) {
			d_H_new[i * 3 + j] = MS::d_h_all(&vertices[i], j);
		}
	}

	double c1 = (H_new - H) / .000001;
	std::cout << "dy H:\tactual: " << c1 << "\tsupposed: " << d_H[1] + d_H[300] + d_H[602] << std::endl;

}