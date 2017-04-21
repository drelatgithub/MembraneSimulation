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
#define RUN_MODE 0

// Wolfe conditions
// Armijo Rule
const double c1 = 0.0001; // Inequality relaxation
const double tau = 0.1; // Shrink size of alpha after each iteration
// Curvature condition
const double c2 = 0.1;

const double h_eps = 1e-8; // Maximum tolerance for forces
const double d_eps = 1e-8; // Maximum tolerance for coordinates
const double max_move = 5e-8; // Maximum displacement for each step in any direction


int minimization(std::vector<MS::vertex*> &vertices);

double line_search(std::vector<MS::vertex*> &vertices, int N, double H, double &H_new, double* p, double d_H_max, double *d_H_new, double m, double &m_new, double alpha0);

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

	std::ofstream p_out, f_out, a_out;
	p_out.open("F:\\p_out.txt");
	f_out.open("F:\\f_out.txt");
	a_out.open("F:\\a_out.txt");


	switch (RUN_MODE) {
	case 0:
		for (double a = 0.9e-6; a < 1.2e-6; a += 0.02e-6) {
			std::cout << update_len(a) << std::endl;

			minimization(vertices);

			for (int i = 0; i < len; i++) {
				p_out << vertices[i]->point->x << '\t' << vertices[i]->point->y << '\t' << vertices[i]->point->z << '\t';
				for (int j = 0; j < 3; j++) {
					f_out << MS::d_h_all(vertices[i], j) << '\t';
				}
				a_out << vertices[i]->area << '\t' << vertices[i]->area0<<'\t';
			}
			p_out << std::endl;
			f_out << std::endl;
			a_out << std::endl;

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
	a_out.close();

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

	while (true) {
		// Find alpha and update coordinates
		double d_H_max = 0.0;
		for (int i=0; i<3*N; i++){
			if(d_H_max < abs(d_H[i])) d_H_max = abs(d_H[i]);
		}
		if(d_H_max < h_eps) break; // Force is almost zero
		alpha0 = max_move / d_H_max; // This ensures that no vertex would have greater step than max_move
		std::cout << "[MIN] Max gradient: " << d_H_max <<"\talpha0: "<<alpha0<< std::endl;

		// m is the inner product of the gradient and the search direcion, and must be non-negative.
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
				std::cout << "[MIN] Warning: along search direction is increasing energy. Reassigning search direction.\n";
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

		std::cout << "[MIN] H: " << H << "\tm: " << m << std::endl;


		if (false) { // Data verification
			int vinds[] = { 0,100,200,1078 };
			for (int i = 0; i < 4; i++) {
				MS::vertex* v = vertices[vinds[i]];
				std::cout << vinds[i] << "\tarea: " << v->area<< std::endl;
				std::cout << std::endl;
			}
		}

		alpha = line_search(vertices, N, H, H_new, p, d_H_max, d_H_new, m, m_new, alpha0);
		std::ofstream t1;
		t1.open("F:\\t1.txt");
		for (int i = 0; i < 3 * N; i++) {
			t1 << p[i] << '\t' << d_H[i] << '\t' << d_H_new[i] << std::endl;
		}
		std::cout << "[MIN] t1 complete." << std::endl;
		t1.close();
		//std::cout << "New! Hn-H-c1*a*m=" << H_new - H - c1*alpha*m << "\t|mn|+c2*m=" << abs(m_new) + c2*m << std::endl;
		// H_new and d_H_new are updated.

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

		for (int i = 0; i < N; i++) {
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
double line_search(std::vector<MS::vertex*> &vertices, int N, double H, double &H_new, double* p, double d_H_max, double *d_H_new, double m, double &m_new, double alpha0) {

	double MIN_D_ALPHA_FAC = 1e-15; // Minimum delta

	double alpha_init = -0.1 * abs(H) / m;
	double alpha = 0;
	double d_alpha = std::fmin(alpha_init, alpha0 * 0.5); // To ensure that 1st alpha is not larger than alpha0
	double H_p = H, m_p = m;
	bool accepted = false;
	bool backtrack = true;

	while (true) {
		accepted = true;
		backtrack = false;

		alpha += d_alpha;
		if(alpha > alpha0){
			alpha -= d_alpha;
			return alpha; // Ensure this won't happen for the 1st iteration.
		}

		for (int i = 0; i < N; i++) {
			vertices[i]->point->x = vertices[i]->point_last->x + alpha * p[i * 3];
			vertices[i]->point->y = vertices[i]->point_last->y + alpha * p[i * 3 + 1];
			vertices[i]->point->z = vertices[i]->point_last->z + alpha * p[i * 3 + 2];
		}
		
		for (int i = 0; i < N; i++) {
			vertices[i]->update_geo();
			if (vertices[i]->area <= 0) {
				accepted = false; // define H = infty, also need to backtrack
			}
		}
		if (!accepted) {
			backtrack = true;
			std::cout << "[LINE] [BACK] area" << std::endl;
		}

		// Renew energy
		H_new = 0;
		for (int i = 0; i < N; i++) {
			H_new += MS::h_all(vertices[i]);
		}

		if (USE_LINE_SEARCH && (H_new >= H_p)) { // Armijo condition not satisfied. simply taking c1=0
			std::cout << "[LINE] [BACK] energy" << std::endl;
			backtrack = true;
		} // Armijo condition satisfied.

		// Renew energy derivatives and m value
		m_new = 0;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < 3; j++) {
				d_H_new[i * 3 + j] = MS::d_h_all(vertices[i], j);
				m_new += p[i * 3 + j] * d_H_new[i * 3 + j];
			}
		}
		if (m_new > 0) {
			std::cout << "[LINE] [BACK] force" << std::endl;
			backtrack = true;
		}
		
		std::cout << "[LINE] H_new: " << H_new << "\tm_new: " << m_new << std::endl;

		if(backtrack){
			alpha -= d_alpha;

			if(false && m_new > 0){ // The force has changed sign
				d_alpha *= m_p / (m_p - m_new); // Linearized force profile
			}else{
				d_alpha *= tau; // Just to make it small
			}

			// Consider cases where moves are simply too small
			if (d_H_max * d_alpha <= MIN_D_ALPHA_FAC) {
				std::cout << "[LINE] d_alpha too small" << std::endl;
				return alpha;
			}

			continue;

			// H_p and m_p not updated, because we moved back.

		}

		// No backtracking case

		if (false && abs(H_new - 1.06455e-13) < 0.0001e-13) {
			//test derivative

			for (int ind = 0; ind < N; ind++) {
				
				vertices[ind]->point->x = vertices[ind]->point_last->x + (alpha+2e-2) * p[ind * 3];
				vertices[ind]->point->y = vertices[ind]->point_last->y + (alpha + 2e-2) * p[ind * 3 + 1];
				vertices[ind]->point->z = vertices[ind]->point_last->z + (alpha + 2e-2) * p[ind * 3 + 2];

				vertices[ind]->update_geo();

				for (int i = 0; i < vertices[ind]->neighbours; i++) {
					vertices[ind]->n[i]->update_geo();
				}

				// Renew energy
				double H_new_n = 0;
				for (int i = 0; i < N; i++) {
					H_new_n += MS::h_all(vertices[i]);
				}

				// Renew energy derivatives and m value
				double m_new_n = 0;
				for (int i = 0; i < N; i++) {
					for (int j = 0; j < 3; j++) {
						d_H_new[i * 3 + j] = MS::d_h_all(vertices[i], j);
						m_new_n += p[i * 3 + j] * d_H_new[i * 3 + j];
					}
				}
				double a = 2e-2 * (p[ind * 3] * d_H_new[ind * 3]+ p[ind * 3+1] * d_H_new[ind * 3+1]+ p[ind * 3+2] * d_H_new[ind * 3+2]);

				std::cout << ind << "\tDE: " << H_new_n - H_new << "\ta: " << a << std::endl;
				H_new = H_new_n;
			}
		}
		
		if (!USE_LINE_SEARCH || abs(m_new) <= -c2 * m) { // Curvature condition satisfied. Good.
			return alpha;
		} // Curvature condition not satisfied

		// Getting a new alpha using linearized force
		double boostFactor = 4;
		if (m_p < m_new) boostFactor = m_new / (m_p - m_new);
		if (boostFactor > 4) boostFactor = 4;
		d_alpha *= boostFactor;

		std::cout << "[LINE] alpha: " << alpha << "\td_alpha: " << d_alpha << std::endl;

		// start over
		m_p = m_new;
		H_p = H_new;

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
	/*
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
	*/

	// Test local free energies with neighbours
	
	int vind = 1078;
	double increment = 0.00001;

	double l1 = MS::h_curv_h(vertices[vind]), l2 = MS::h_pressure(vertices[vind]);
	for (int i = 0; i < vertices[vind]->neighbours; i++) {
		l1 += MS::h_curv_h(vertices[vind]->n[i]);
		l2 += MS::h_pressure(vertices[vind]->n[i]);
	}
	double d1 = MS::d_h_curv_h(vertices[vind], 1), d2 = MS::d_h_pressure(vertices[vind], 1);


	vertices[vind]->point->y += increment;
	vertices[vind]->update_geo();
	for (int i = 0; i < vertices[vind]->neighbours; i++) {
		vertices[vind]->n[i]->update_geo();
	}

	double n1 = MS::h_curv_h(vertices[vind]), n2 = MS::h_pressure(vertices[vind]);
	for (int i = 0; i < vertices[vind]->neighbours; i++) {
		n1 += MS::h_curv_h(vertices[vind]->n[i]);
		n2 += MS::h_pressure(vertices[vind]->n[i]);
	}
	double c1 = (n1 - l1) / increment,
		c2 = (n2 - l2) / increment;

	std::cout << "dy h:\tactual: " << c1 << "\tsupposed: " << d1 << std::endl;
	std::cout << "dy p:\tactual: " << c2 << "\tsupposed: " << d2 << std::endl;
	

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
