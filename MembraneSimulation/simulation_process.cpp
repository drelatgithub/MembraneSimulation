#include"simulation_process.h"

#include<vector>
#include<math.h>
#include<iostream>
#include<fstream>

#include"surface_mesh.h"
#include"surface_free_energy.h"

#define USE_STEEPEST_DESCENT true
#define USE_LINE_SEARCH true

/* RUN_MODE
	0: Normal simulation
	1: Derivative test
	2: Direction force profile of one vertex
*/
#define RUN_MODE 1

// Wolfe conditions
// Armijo Rule
const double c1 = 0.0001; // Inequality relaxation
const double tau = 0.5; // Shrink size of alpha after each iteration
// Curvature condition
const double c2 = 0.1;

const double eps = 5e-3; // Maximum deviation for each coordinate

int minimization(std::vector<MS::vertex*> &vertices);

double line_search(std::vector<MS::vertex*> &vertices, int N, double H, double &H_new, double* p, double *d_H_new, double m, double &m_new, double alpha0);
double zoom(std::vector<MS::vertex*> &vertices, int N, double H, double &H_new, double* p, double *d_H_new, double m, double &m_new, double alphal, double alphah);

void test_derivatives(std::vector<MS::vertex*> &vertices);
void force_profile(std::vector<MS::vertex*> &vertices);

int MS::simulation_start(std::vector<vertex*> &vertices) {
	int len = vertices.size();

	for (int i = 0; i < len; i++) {
		vertices[i]->count_neighbours();
		if (i == 0) {
			int a = 1;
		}
		vertices[i]->update_geo();
		vertices[i]->make_initial();
	}

	std::ofstream p_out, f_out;
	p_out.open("F:\\p_out.txt");
	f_out.open("F:\\f_out.txt");


	switch (RUN_MODE) {
	case 0:
		for (double a = -1.3; a < 1.5; a += 0.05) {
			std::cout << update_len(a) << std::endl;

			minimization(vertices);

			for (int i = 0; i < len; i++) {
				p_out << vertices[i]->point->x << '\t' << vertices[i]->point->y << '\t' << vertices[i]->point->z << '\t';
				for (int j = 0; j < 3; j++) {
					f_out << MS::d_h_all(vertices[i], j) << '\t';
				}
			}
			p_out << std::endl;
			f_out << std::endl;

			break;
		}

		std::cout << update_len(-3) << std::endl;

		minimization(vertices);

		for (int i = 0; i < len; i++) {
			p_out << vertices[i]->point->x << '\t' << vertices[i]->point->y << '\t' << vertices[i]->point->z << '\t';
			for (int j = 0; j < 3; j++) {
				f_out << MS::d_h_all(vertices[i], j) << '\t';
			}
		}
		p_out << std::endl;
		f_out << std::endl;
		break;

	case 1:
		test_derivatives(vertices);
		break;

	case 2:
		force_profile(vertices);
		break;
	}
	

	p_out.close();
	f_out.close();

	return 0;
}

int minimization(std::vector<MS::vertex*> &vertices) {
	bool finished = false;
	int N = vertices.size(); // Number of vertices
	double H = 0, H_new = 0;
	double *d_H = new double[3 * N];
	double *d_H_new = new double[3 * N];
	double *p = new double[3 * N]; // Search direction
	double alpha0;
	double alpha;
	double beta;

	std::ofstream p_min_out, f_min_out, sd_min_out;
	p_min_out.open("F:\\p_min_out.txt");
	f_min_out.open("F:\\f_min_out.txt");
	sd_min_out.open("F:\\sd_min_out.txt");
	
	// Initializing
	for (int i = 0; i < N; i++) {
		// Get H and d_H
		H += MS::h_all(vertices[i]);
		for (int j = 0; j < 3; j++) {
			d_H[i * 3 + j] = MS::d_h_all(vertices[i], j);
			// Initialize search direction
			p[i * 3 + j] = -d_H[i * 3 + j];
		}

		// Store vertices location as the last location
		vertices[i]->make_last();

		//temp_out << MS::d_h_curv_g(&vertices[i], 0) << '\t' << MS::d_h_curv_g(&vertices[i], 1) << '\t' << MS::d_h_curv_g(&vertices[i], 2) << '\t';
	}
	//temp_out.close();

	int k = 0; // Iteration counter. Only for debug use.

	while (!finished) {
		// Find alpha and update coordinates
		alpha0 = 0.001; // TODO initial alpha?

		double m = 0, m_new = 0;
		if (USE_STEEPEST_DESCENT) {
			for (int i = 0; i < 3 * N; i++) {
				p[i] = -d_H[i];
				m += p[i] * d_H[i];
			}
		}
		else {
			for (int i = 0; i < 3*N; i++) {
				m += p[i] * d_H[i];
			}
			if (m > 0) {
				std::cout << "Warning: along search direction is increasing energy. Reassigning search direction.\n";
				m = 0;
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < 3; j++) {
						// Initialize search direction
						p[i * 3 + j] = -d_H[i * 3 + j];
						m += p[i * 3 + j] * d_H[i * 3 + j];
					}
				}
			}
		}

		alpha = line_search(vertices, N, H, H_new, p, d_H_new, m, m_new, alpha0);
		//std::cout << "New! Hn-H-c1*a*m=" << H_new - H - c1*alpha*m << "\t|mn|+c2*m=" << abs(m_new) + c2*m << std::endl;
		// H_new and d_H_new are updated.

		double aaa = MS::d_h_all(vertices[1078], 0);
		aaa = MS::d_h_all(vertices[1078], 1);
		aaa = MS::d_h_all(vertices[1078], 2);
		MS::vertex* nnn = vertices[1078]->n[0];
		aaa = MS::d_h_all(nnn, 0);
		aaa = MS::d_h_all(nnn, 1);
		aaa = MS::d_h_all(nnn, 2);

		if (!USE_STEEPEST_DESCENT) {
			// Find beta (Fletcher-Reeves)
			double a = 0, b = 0;
			for (int i = 0; i < 3 * N; i++) {
				a += d_H_new[i] * d_H_new[i];
				b += d_H[i] * d_H[i];
			}
			beta = a / b;
			// Find beta (Polak-Ribiere)
			//double a = 0, b = 0;
			//for (int i = 0; i < 3 * N; i++) {
			//	a += d_H_new[i] * (d_H_new[i] - d_H[i]);
			//	b += d_H[i] * d_H[i];
			//}
			//beta = (a >= 0) ? a / b : 0;

			// Renew search direction
			for (int i = 0; i < 3 * N; i++) {
				p[i] = -d_H_new[i] + beta*p[i];
			}
		}

		// Judge when we should exit loop
		finished = true;
		for (int i = 0; i < N; i++) {
			if (finished && (
				abs(d_H_new[i * 3 + 0]) > eps || abs(d_H_new[i * 3 + 1]) > eps || abs(d_H_new[i * 3 + 2]) > eps
				//vertices[i]->point->x - vertices[i]->point_last->x > eps || vertices[i]->point_last->x - vertices[i]->point->x > eps ||
				//vertices[i]->point->y - vertices[i]->point_last->y > eps || vertices[i]->point_last->y - vertices[i]->point->y > eps ||
				//vertices[i]->point->z - vertices[i]->point_last->z > eps || vertices[i]->point_last->z - vertices[i]->point->z > eps
				)) {
				finished = false;
			}
			// Record coordinates as last-time coordinates
			vertices[i]->make_last();
		}

		// Renew H and d_H
		H = H_new;
		for (int i = 0; i < 3 * N; i++) {
			d_H[i] = d_H_new[i];
		}

		k++;
		std::cout << k << "\tm: " << m << "\talpha: " << alpha << "\tFree energy: " << H << std::endl;
		for (int i = 0; i < N; i++) {
			p_min_out << vertices[i]->point->x << '\t' << vertices[i]->point->y << '\t' << vertices[i]->point->z << '\t';
			for (int j = 0; j < 3; j++) {
				f_min_out << d_H[i * 3 + j] << '\t';
				sd_min_out << p[i * 3 + j] << '\t';
			}
		}
		p_min_out << std::endl;
		f_min_out << std::endl;
		sd_min_out << std::endl;
	}

	f_min_out.close();
	p_min_out.close();
	sd_min_out.close();

	delete[]d_H;
	delete[]d_H_new;
	delete[]p;

	return 0;
}
double line_search(std::vector<MS::vertex*> &vertices, int N, double H, double &H_new, double* p, double *d_H_new, double m, double &m_new, double alpha0) {

	double alpha = alpha0 * tau;
	double alphap = 0;
	double Hp = H;
	bool accepted = false;

	while (!accepted) {
		accepted = true;

		for (int i = 0; i < N; i++) {
			vertices[i]->point->x = vertices[i]->point_last->x + alpha*p[i * 3];
			vertices[i]->point->y = vertices[i]->point_last->y + alpha*p[i * 3 + 1];
			vertices[i]->point->z = vertices[i]->point_last->z + alpha*p[i * 3 + 2];
		}
		
		for (int i = 0; i < N; i++) {
			vertices[i]->update_geo();
			if (vertices[i]->area <= 0) {
				accepted = false;
			}
		}
		if (!accepted) {
			alpha = alphap + (alpha - alphap) * tau;
			continue;
		}

		// Renew energy
		H_new = 0;
		for (int i = 0; i < N; i++) {
			H_new += MS::h_all(vertices[i]);
		}

		if (USE_LINE_SEARCH && (H_new >= Hp || H_new - H > c1*alpha*m)) { // Armijo condition not satisfied.
			std::cout << "Too far!\talphap: " << alphap << "\talpha: " << alpha << std::endl;
			alpha = zoom(vertices, N, H, H_new, p, d_H_new, m, m_new, alphap, alpha);
			return alpha;
		} // Armijo condition satisfied.

		// Renew energy derivatives and m value
		m_new = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < 3; j++) {
				d_H_new[i * 3 + j] = MS::d_h_all(vertices[i], j);
				m_new += p[i * 3 + j] * d_H_new[i * 3 + j];
			}
		}
		
		if (!USE_LINE_SEARCH || abs(m_new) <= -c2 * m) { // Curvature condition satisfied. Good.
			return alpha;
		} // Curvature condition not satisfied

		if (m_new >= 0) { // Minimum is inside
			std::cout << "Minimum is inside!\talpha: " << alpha << "\talphap: " << alphap << std::endl;
			alpha = zoom(vertices, N, H, H_new, p, d_H_new, m, m_new, alpha, alphap);
			return alpha;
		}

		// Increase alpha and try again
		accepted = false;
		Hp = H_new;
		alphap = alpha;
		alpha += (alpha0 - alpha) * tau;

		std::cout << "alphap: " << alphap << "\talpha: " << alpha << std::endl;

	}
}
double zoom(std::vector<MS::vertex*> &vertices, int N, double H, double &H_new, double* p, double *d_H_new, double m, double &m_new, double alphal, double alphah) {
	double alpha = 0.5 * (alphal + alphah);
	bool accepted = false;
	while (!accepted) {

		accepted = true;

		for (int i = 0; i < N; i++) {
			vertices[i]->point->x = vertices[i]->point_last->x + alpha*p[i * 3];
			vertices[i]->point->y = vertices[i]->point_last->y + alpha*p[i * 3 + 1];
			vertices[i]->point->z = vertices[i]->point_last->z + alpha*p[i * 3 + 2];
		}
		for (int i = 0; i < N; i++) {
			vertices[i]->update_geo();
			if (vertices[i]->area <= 0) {
				accepted = false;
			}
		}
		if (!accepted) {
			double alpha_small = ((alphal < alphah) ? alphal : alphah);
			alpha = alpha_small + (alpha - alpha_small) * tau;
			continue;
		}
		accepted = false;

		// Renew energy
		H_new = 0;
		for (int i = 0; i < N; i++) {
			H_new += MS::h_all(vertices[i]);
		}

		if (H_new - H > c1*alpha*m) { // Armijo condition not satisfied.
			alphah = alpha;
		}
		else {

			// Renew energy derivatives and m value
			m_new = 0;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < 3; j++) {
					d_H_new[i * 3 + j] = MS::d_h_all(vertices[i], j);
					m_new += p[i * 3 + j] * d_H_new[i * 3 + j];
				}
			}

			if (abs(m_new) <= -c2*m) {
				std::cout << "m_new: " << m_new << "\tm:" << m << std::endl;
				return alpha;
			}

			if (m_new * (alphah - alphal) >= 0) {
				alphah = alphal;
			}
			alphal = alpha;
		}

		alpha = 0.5 * (alphal + alphah);

		std::cout << "alphal: " << alphal << "\talphah: " << alphah << std::endl;
	}
}

void test_derivatives(std::vector<MS::vertex*> &vertices) {
	// Test local properties (main)
	/*
	int vind = 1078;
	double increment = 0.0001;

	double l1 = vertices[vind]->curv_h, l2 = vertices[vind]->area, l3 = vertices[vind]->cot_theta2[0];
	double d1 = vertices[vind]->dx_curv_h, d2 = vertices[vind]->dx_area, d3 = vertices[vind]->dx_cot_theta2[0];

	vertices[vind]->point->x += increment;
	vertices[vind]->update_geo();

	double c1 = (vertices[vind]->curv_h - l1) / increment,
		c2 = (vertices[vind]->area - l2) / increment,
		c3 = (vertices[vind]->cot_theta2[0] - l3) / increment;
	std::cout << "dx h:\tactual: " << c1 << "\tsupposed: " << d1 << std::endl;
	std::cout << "dx a:\tactual: " << c2 << "\tsupposed: " << d2 << std::endl;
	std::cout << "dx t2:\tactual: " << c3 << "\tsupposed: " << d3 << std::endl;

	*/

	// Test local properties (neighbours)
	int vind = 1078;
	double increment = 0.0001;

	MS::vertex *n = vertices[vind]->n[0], *nn = vertices[vind]->n_next[0], *np = vertices[vind]->n_prev[0];
	//double l1 = vertices[0].curv_h, l2 = vertices[0].curv_g, l3 = vertices[0].area, l4 = vertices[0].cot_theta3[0];
	double l1 = vertices[vind]->cot_theta3[0], l2 = vertices[vind]->theta3[0];
	//double d1 = vertices[0].dyn_curv_h[0], d2 = vertices[0].dyn_curv_g[0], d3 = vertices[0].dyn_area[0], d4=vertices[0].dyn_cot_theta3[0];
	double d1 = vertices[vind]->dznn_cot_theta3[0], d2 = vertices[vind]->dznn_theta3[0];

	//n->point->y += increment;
	nn->point->z += increment;
	vertices[vind]->update_geo();

	//double c1 = (vertices[0].curv_h - l1) / increment,
		//c2 = (vertices[0].curv_g - l2) / increment,
		//c3 = (vertices[0].area - l3) / increment,
		//c4 = (vertices[0].cot_theta3[0] - l4) / increment;
	double c1 = (vertices[vind]->cot_theta3[0] - l1) / increment, c2 = (vertices[vind]->theta3[0] - l2) / increment;
	std::cout << "dz:\tactual: " << c1 << "\tsupposed: " << d1 << std::endl;
	std::cout << "dz:\tactual: " << c2 << "\tsupposed: " << d2 << std::endl;
	//std::cout << "dy a:\tactual: " << c3 << "\tsupposed: " << d3 << std::endl;
	//std::cout << "dy r:\tactual: " << c4 << "\tsupposed: " << d4 << std::endl;

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
	/*
	int N = vertices.size();
	double H = 0, H_new = 0;
	double *d_H = new double[3 * N];
	double *d_H_new = new double[3 * N];

	for (int i = 0; i < N; i++) {
		H += MS::h_all(vertices[i]);
		for (int j = 0; j < 3; j++) {
			d_H[i * 3 + j] = MS::d_h_all(vertices[i], j);
		}
	}

	vertices[0]->point->y += 0.000001;
	vertices[100]->point->x += 0.000001;
	vertices[200]->point->z += 0.000001;
	for (int i = 0; i < N; i++) {
		vertices[i]->update_geo();
	}

	for (int i = 0; i < N; i++) {
		H_new += MS::h_all(vertices[i]);
		for (int j = 0; j < 3; j++) {
			d_H_new[i * 3 + j] = MS::d_h_all(vertices[i], j);
		}
	}

	double c1 = (H_new - H) / .000001;
	std::cout << "dy H:\tactual: " << c1 << "\tsupposed: " << d_H[1] + d_H[300] + d_H[602] << std::endl;
	*/

}

void force_profile(std::vector<MS::vertex*> &vertices) {
	int v_index = 10;
	int N = vertices.size();
	int len = vertices[v_index]->n.size();

	std::ofstream nfp, lfp1;
	nfp.open("F:\\nfp.txt");
	lfp1.open("F:\\lfp1.txt");

	vertices[v_index]->update_geo();
	double H = MS::h_all(vertices[v_index]);
	for (int j = 0; j < len; j++) {
		vertices[v_index]->n[j]->update_geo();
		H += MS::h_all(vertices[v_index]->n[j]);
	}

	vertices[v_index]->make_last();
	double n_x = vertices[v_index]->n_x;
	double n_y = vertices[v_index]->n_y;
	double n_z = vertices[v_index]->n_z;
	double l1_x = -n_y;
	double l1_y = n_x;
	double l1_z = 0;

	double H_new;

	for (double move = -0.02; move < 0.02; move += 0.001) {
		vertices[v_index]->point->x = vertices[v_index]->point_last->x + n_x*move;
		vertices[v_index]->point->y = vertices[v_index]->point_last->y + n_y*move;
		vertices[v_index]->point->z = vertices[v_index]->point_last->z + n_z*move;
		vertices[v_index]->update_geo();
		H_new = 0;
		H_new = MS::h_all(vertices[v_index]);
		for (int j = 0; j < len; j++) {
			vertices[v_index]->n[j]->update_geo();
			H_new += MS::h_all(vertices[v_index]->n[j]);
		}
		//std::cout << "move: " << move << "\tdH: " << H_new - H << std::endl;
		nfp << move << '\t' << H_new - H << std::endl;
	}

	for (double move = -0.08; move < 0.04; move += 0.001) {
		vertices[v_index]->point->x = vertices[v_index]->point_last->x + l1_x*move;
		vertices[v_index]->point->y = vertices[v_index]->point_last->y + l1_y*move;
		vertices[v_index]->point->z = vertices[v_index]->point_last->z + l1_z*move;
		vertices[v_index]->update_geo();
		H_new = 0;
		H_new = MS::h_all(vertices[v_index]);
		for (int j = 0; j < len; j++) {
			vertices[v_index]->n[j]->update_geo();
			H_new += MS::h_all(vertices[v_index]->n[j]);
		}
		//std::cout << "move: " << move << "\tdH: " << H_new - H << std::endl;
		lfp1 << move << '\t' << H_new - H << std::endl;
	}

	nfp.close();
	lfp1.close();
}