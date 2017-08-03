#define _USE_MATH_DEFINES

#include"simulation_process.h"

#include"common.h"
#include"math_public.h"
#include"surface_mesh.h"
#include"surface_mesh_tip.h"

#define USE_STEEPEST_DESCENT false
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

const double h_eps = 1e-12; // Maximum tolerance for forces
const double d_eps = 1e-8; // Maximum tolerance for coordinates
const double max_move = 5e-8; // Maximum displacement for each step in any direction


int minimization(MS::surface_mesh &sm);

double line_search(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets, int N, int N_f, double H, double &H_new, double* p, double d_H_max, double *d_H_new, double m, double &m_new, double alpha0);

void test_derivatives(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets);
void force_profile(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets);

void update_all_energy(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets, bool update_geometry=true, bool update_energy=true) {
	int len = vertices.size(), len_f = facets.size();
	if (update_geometry) {
		for (int i = 0; i < len_f; i++) {
			facets[i]->update_geo();
		}
		for (int i = 0; i < len; i++) {
			vertices[i]->update_geo();
		}
	}
	if (update_energy) {
		for (int i = 0; i < len; i++) {
			vertices[i]->update_energy();
		}
		for (int i = 0; i < len_f; i++) {
			facets[i]->update_energy(MS::po);
			facets[i]->inc_H_int(MS::po); // This could also change energy derivatives in vertices
		}
	}
}
int MS::simulation_start(MS::surface_mesh &sm, std::vector<MS::filament_tip*> &tips) {
	auto &vertices = sm.vertices;
	auto &facets = sm.facets;

	int N = vertices.size(),
		N_f = facets.size();

	for (int i = 0; i < N; i++) {
		vertices[i]->count_neighbors();
		if (i == 0) {
			int a = 1;
		}
		vertices[i]->update_geo();
		vertices[i]->make_initial();
	}

	{ // Doing some statistics
		double total_edge_length = 0, total_area = 0, total_edge_length_sq = 0, total_area_sq = 0;
		int num_edge2 = 0;
		for (int i = 0; i < N; i++) {
			total_area += vertices[i]->area;
			total_area_sq += vertices[i]->area*vertices[i]->area;
			num_edge2 += vertices[i]->neighbors;
			for (int j = 0; j < vertices[i]->neighbors; j++) {
				total_edge_length += vertices[i]->r_p_n[j];
				total_edge_length_sq += vertices[i]->r_p_n[j] * vertices[i]->r_p_n[j];
			}
		}
		double avg_edge_length = total_edge_length / num_edge2,
			avg_edge_length_sq = total_edge_length_sq / num_edge2;
		double stdev_edge_length = sqrt(avg_edge_length_sq - avg_edge_length*avg_edge_length);
		double avg_area = total_area / N,
			avg_area_sq = total_area_sq / N;
		double stdev_area = sqrt(avg_area_sq - avg_area*avg_area);

		std::stringstream ss;
		std::string big_divider(40, '='),
			small_divider(40, '-');
		ss << big_divider << std::endl
			<< "Number of vertices: " << N << std::endl
			<< "Number of facets: " << N_f << std::endl
			<< "Number of edges: " << num_edge2 / 2 << std::endl
			<< small_divider << std::endl
			<< "Total surface area: " << total_area << std::endl
			<< "Diameter if spherical: " << sqrt(total_area / M_PI) << std::endl
			<< "Average vertex area: " << avg_area << std::endl
			<< "Stdev vertex area: " << stdev_area << std::endl
			<< small_divider << std::endl
			<< "Average edge length: " << avg_edge_length << std::endl
			<< "Stdev edge length: " << stdev_edge_length << std::endl
			<< big_divider;
		LOG(INFO) << "Meshwork properties: " << std::endl << ss.str();
	} // End doing statistics


	std::ofstream p_out, f_out, a_out;
	p_out.open("F:\\p_out.txt");
	f_out.open("F:\\f_out.txt");
	a_out.open("F:\\a_out.txt");


	switch (RUN_MODE) {
	case 0:
		for (double a = 0.994e-6; a < 1.010e-6; a += 0.001e-6) {
			// Update filament tip position
			LOG(INFO) << "Polymer tip x position: " << update_len(a);

			minimization(sm);

			update_all_energy(vertices, facets);

			for (int i = 0; i < N; i++) {
				p_out << vertices[i]->point->x << '\t' << vertices[i]->point->y << '\t' << vertices[i]->point->z << '\t';
				math_public::Vec3 cur_d_h_all = vertices[i]->d_H;
				f_out << cur_d_h_all.x << '\t' << cur_d_h_all.y << '\t' << cur_d_h_all.z << '\t';
				a_out << vertices[i]->area << '\t' << vertices[i]->area0<<'\t';
			}
			p_out << std::endl;
			f_out << std::endl;
			a_out << std::endl;

		}

		minimization(sm);
		update_all_energy(vertices, facets);

		for (int i = 0; i < N; i++) {
			p_out << vertices[i]->point->x << '\t' << vertices[i]->point->y << '\t' << vertices[i]->point->z << '\t';
			math_public::Vec3 cur_d_h_all = vertices[i]->d_H;
			f_out << cur_d_h_all.x << '\t' << cur_d_h_all.y << '\t' << cur_d_h_all.z << '\t';
		}
		p_out << std::endl;
		f_out << std::endl;
		break;

	case 1:
		test_derivatives(vertices, facets);
		break;

	case 2:
		force_profile(vertices, facets);
		break;
	}
	

	p_out.close();
	f_out.close();
	a_out.close();

	return 0;
}

int minimization(MS::surface_mesh &sm) {
	/**************************************************************************
		This function uses the conjugate gradient method to do the energy
		minimization for vertices/facets system.
	**************************************************************************/
	auto &vertices = sm.vertices;
	auto &facets = sm.facets;

	bool finished = false;
	int N = vertices.size(); // Number of vertices
	int N_f = facets.size(); // Number of facets
	double H = 0, H_new = 0;
	double *d_H = new double[3 * N];
	double *d_H_new = new double[3 * N];
	double *p = new double[3 * N]; // Search direction
	double alpha0;
	double alpha; // alpha is the "portion" of distance that each vertex should go along the search vector.
	double beta;

	std::ofstream p_min_out, f_min_out, sd_min_out;
	p_min_out.open("F:\\p_min_out.txt");
	f_min_out.open("F:\\f_min_out.txt");
	sd_min_out.open("F:\\sd_min_out.txt");
	
	// Initializing
	update_all_energy(vertices, facets);
	for (int i = 0; i < N; i++) {
		// Get H and d_H
		H += vertices[i]->H; // This is only the vertices part of the energy
		math_public::Vec3 cur_d_h_all = vertices[i]->d_H;
		d_H[i * 3] = cur_d_h_all.x;
		d_H[i * 3 + 1] = cur_d_h_all.y;
		d_H[i * 3 + 2] = cur_d_h_all.z;
		for (int j = 0; j < 3; j++) {
			// Initialize search direction
			p[i * 3 + j] = -d_H[i * 3 + j];
		}

		// Store vertices location as the last location
		vertices[i]->make_last();

		//temp_out << MS::d_h_curv_g(&vertices[i], 0) << '\t' << MS::d_h_curv_g(&vertices[i], 1) << '\t' << MS::d_h_curv_g(&vertices[i], 2) << '\t';
	}
	for (int i = 0; i < N_f; i++) {
		H += facets[i]->H; // Now add the facets part of the energy
	}
	//temp_out.close();

	int k = 0; // Iteration counter. Only for debug use.

	while (true) {
		k++;
		LOG(INFO) << "Iteration " << k << " starting...";

		// Find alpha and update coordinates
		double d_H_max = 0.0;
		for (int i=0; i<3*N; i++){
			if(d_H_max < abs(d_H[i])) d_H_max = abs(d_H[i]);
		}
		if(d_H_max < h_eps) break; // Force is almost zero
		alpha0 = max_move / d_H_max; // This ensures that no vertex would have greater step than max_move
		LOG(INFO) << "Max gradient: " << d_H_max << " alpha0: " << alpha0;

		// m is the inner product of the gradient and the search direcion, and must be non-negative.
		double m = 0, m_new = 0;
		if (USE_STEEPEST_DESCENT) {
			for (int i = 0; i < 3 * N; i++) {
				p[i] = -d_H[i];
				m += p[i] * d_H[i];
			}
		}
		else { // Use conjugate gradient
			for (int i = 0; i < 3*N; i++) {
				m += p[i] * d_H[i];
			}
			if (m > 0) {
				LOG(WARNING) << "Warning: along search direction is increasing energy. Reassigning search direction.";
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

		LOG(INFO) << "Current H: " << H << " m: " << m;

		if (false) { // Data verification
			int vinds[] = { 0,100,200,1078 };
			for (int i = 0; i < 4; i++) {
				MS::vertex* v = vertices[vinds[i]];
				std::cout << vinds[i] << "\tarea: " << v->area<< std::endl;
				std::cout << std::endl;
			}
		}

		alpha = line_search(vertices, facets, N, H, N_f, H_new, p, d_H_max, d_H_new, m, m_new, alpha0);
		// So far, H_new and d_H_new have already been updated in line_search.

		// Temporary debugging output
		std::ofstream t1;
		t1.open("F:\\t1.txt");
		for (int i = 0; i < 3 * N; i++) {
			t1 << p[i] << '\t' << d_H[i] << '\t' << d_H_new[i] << std::endl;
		}
		LOG(DEBUG) << "t1 data dump complete.";
		t1.close();
		//std::cout << "New! Hn-H-c1*a*m=" << H_new - H - c1*alpha*m << "\t|mn|+c2*m=" << abs(m_new) + c2*m << std::endl;
		

		if (!USE_STEEPEST_DESCENT) { // Conjugate gradient method renewal of search direction.
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
			// Store coordinates as last-time coordinates
			vertices[i]->make_last();
		}

		// Renew H and d_H
		H = H_new;
		for (int i = 0; i < 3 * N; i++) {
			d_H[i] = d_H_new[i];
		}

		// Finish off and get ready for the next iteration.
		LOG(INFO) << "H_new: " << H_new << " m_new: " << m_new;
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
double line_search(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets, int N, int N_f, double H, double &H_new, double* p, double d_H_max, double *d_H_new, double m, double &m_new, double alpha0) {
	/**************************************************************************
	Purpose:
		This function does the line search for a given search direction.

	Parameters:
		alpha0: max value that alpha could take.
	**************************************************************************/

	double MIN_D_ALPHA_FAC = 1e-15; // Minimum delta

	double alpha_init = -0.1 * abs(H) / m; // In m^2/J
	double alpha = 0;
	double d_alpha = std::fmin(alpha_init, alpha0 * 0.5); // To ensure that 1st alpha is not larger than alpha0
	double H_p = H, m_p = m; // Previous H and m
	bool accepted = false;
	bool backtrack = true; // Whether alpha is too big for some reason so that we need to decrease alpha and do further iterations

	while (true) {
		accepted = true;
		backtrack = false;

		alpha += d_alpha;
		if(alpha > alpha0){
			alpha -= d_alpha;
			LOG(INFO) << "Returning alpha as " << alpha << " as it reaches maximum";
			return alpha; // Ensure this won't happen for the 1st iteration, because we cannot let alpha to be zero.
		}

		// Change the position and renew energy
		for (int i = 0; i < N; i++) {
			vertices[i]->point->x = vertices[i]->point_last->x + alpha * p[i * 3];
			vertices[i]->point->y = vertices[i]->point_last->y + alpha * p[i * 3 + 1];
			vertices[i]->point->z = vertices[i]->point_last->z + alpha * p[i * 3 + 2];
		}
		update_all_energy(vertices, facets);

		// Make sure that area is not negative
		for (int i = 0; i < N; i++) {
			if (vertices[i]->area <= 0) {
				accepted = false;
				backtrack = true; // Because H = infty, we also need to do backtracking
				LOG(INFO) << "[BACKTRACK] Area is negative.";
				break;
			}
		}

		// Renew the sum of energy
		H_new = 0;
		for (int i = 0; i < N; i++) {
			H_new += vertices[i]->H;
		}
		for (int i = 0; i < N_f; i++) {
			H_new += facets[i]->H;
		}

		if (USE_LINE_SEARCH && (H_new >= H_p)) { // Armijo condition not satisfied. simply taking c1=0
			LOG(INFO) << "[BACKTRACK] Energy is increasing.";
			backtrack = true;
		} // Otherwise, Armijo condition is satisfied.

		// Renew energy derivatives and m value
		m_new = 0;
		for (int i = 0; i < N; i++) {
			math_public::Vec3 cur_d_h_all = vertices[i]->d_H;
			d_H_new[i * 3] = cur_d_h_all.x;
			d_H_new[i * 3 + 1] = cur_d_h_all.y;
			d_H_new[i * 3 + 2] = cur_d_h_all.z;
			for (int j = 0; j < 3; j++) {
				m_new += p[i * 3 + j] * d_H_new[i * 3 + j];
			}
		}
		if (m_new > 0) {
			LOG(INFO) << "[BACKTRACK] New force along search direction.";
			backtrack = true;
		}
		
		LOG(DEBUG) << "H_new: " << H_new << " m_new: " << m_new;

		if(backtrack){
			alpha -= d_alpha; // Get back to last alpha

			if(false && m_new > 0){ // The force has changed sign
				d_alpha *= m_p / (m_p - m_new); // Linearized force profile. Currently we don't use that.
			}else{
				d_alpha *= tau; // Just to make it small
			}

			// Consider cases where moves are simply too small
			if (d_H_max * d_alpha <= MIN_D_ALPHA_FAC) {
				LOG(INFO) << "Returning alpha as " << alpha << " as d_alpha is too small";
				if (alpha == 0.0) LOG(WARNING) << "d_alpha is too small, and returned alpha is zero.";
				return alpha;
			}

			continue;

			// H_p and m_p not updated, because we moved back.

		}

		// No backtracking case

		if (false && abs(H_new - 1.06455e-13) < 0.0001e-13) {
			//test derivative
			// TODO: Consider facet interactions

			for (int ind = 0; ind < N; ind++) {
				
				vertices[ind]->point->x = vertices[ind]->point_last->x + (alpha+2e-2) * p[ind * 3];
				vertices[ind]->point->y = vertices[ind]->point_last->y + (alpha + 2e-2) * p[ind * 3 + 1];
				vertices[ind]->point->z = vertices[ind]->point_last->z + (alpha + 2e-2) * p[ind * 3 + 2];

				vertices[ind]->update_geo();
				for (int i = 0; i < vertices[ind]->neighbors; i++) {
					vertices[ind]->n[i]->update_geo();
				}
				vertices[ind]->update_energy();
				for (int i = 0; i < vertices[ind]->neighbors; i++) {
					vertices[ind]->n[i]->update_energy();
				}

				// Renew energy
				double H_new_n = 0;
				for (int i = 0; i < N; i++) {
					H_new_n += vertices[i]->H;
				}

				// Renew energy derivatives and m value
				double m_new_n = 0;
				for (int i = 0; i < N; i++) {
					math_public::Vec3 cur_d_h_all = vertices[i]->d_H;
					d_H_new[i * 3] = cur_d_h_all.x;
					d_H_new[i * 3 + 1] = cur_d_h_all.y;
					d_H_new[i * 3 + 2] = cur_d_h_all.z;
					for (int j = 0; j < 3; j++) {
						m_new_n += p[i * 3 + j] * d_H_new[i * 3 + j];
					}
				}
				double a = 2e-2 * (p[ind * 3] * d_H_new[ind * 3]+ p[ind * 3+1] * d_H_new[ind * 3+1]+ p[ind * 3+2] * d_H_new[ind * 3+2]);

				std::cout << ind << "\tDE: " << H_new_n - H_new << "\ta: " << a << std::endl;
				H_new = H_new_n;
			}
		} // End of derivative test
		
		if (!USE_LINE_SEARCH || abs(m_new) <= -c2 * m) { // Curvature condition satisfied. Good.
			LOG(INFO) << "Returning alpha as " << alpha << " as it fits search criteria.";
			return alpha;
		} // Curvature condition not satisfied

		// Getting a new alpha using linearized force
		double boostFactor = 4;
		if (m_p < m_new) boostFactor = m_new / (m_p - m_new);
		if (boostFactor > 4) boostFactor = 4;
		d_alpha *= boostFactor;

		LOG(DEBUG) << "alpha: " << alpha << " d_alpha: " << d_alpha<<" before next iteration.";

		// start over
		m_p = m_new;
		H_p = H_new;

	}
}

void test_derivatives(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets) {
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

	//// Test local free energies with neighbours
	//
	//int vind = 1078;
	//double increment = 0.00001;

	//double l1 = MS::h_curv_h(vertices[vind]), l2 = MS::h_pressure(vertices[vind]);
	//for (int i = 0; i < vertices[vind]->neighbors; i++) {
	//	l1 += MS::h_curv_h(vertices[vind]->n[i]);
	//	l2 += MS::h_pressure(vertices[vind]->n[i]);
	//}
	//double d1 = MS::d_h_curv_h(vertices[vind]).y, d2 = MS::d_h_pressure(vertices[vind]).y;


	//vertices[vind]->point->y += increment;
	//vertices[vind]->update_geo();
	//for (int i = 0; i < vertices[vind]->neighbors; i++) {
	//	vertices[vind]->n[i]->update_geo();
	//}

	//double n1 = MS::h_curv_h(vertices[vind]), n2 = MS::h_pressure(vertices[vind]);
	//for (int i = 0; i < vertices[vind]->neighbors; i++) {
	//	n1 += MS::h_curv_h(vertices[vind]->n[i]);
	//	n2 += MS::h_pressure(vertices[vind]->n[i]);
	//}
	//double c1 = (n1 - l1) / increment,
	//	c2 = (n2 - l2) / increment;

	//std::cout << "dy h:\tactual: " << c1 << "\tsupposed: " << d1 << std::endl;
	//std::cout << "dy p:\tactual: " << c2 << "\tsupposed: " << d2 << std::endl;
	//

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

void force_profile(std::vector<MS::vertex*> &vertices, std::vector<MS::facet*> &facets) {
	// TODO: Consider facet interactions
	int v_index = 10;
	int N = vertices.size();
	int len = vertices[v_index]->n.size();

	std::ofstream nfp, lfp1;
	nfp.open("F:\\nfp.txt");
	lfp1.open("F:\\lfp1.txt");

	vertices[v_index]->update_geo();
	double H = vertices[v_index]->H;
	for (int j = 0; j < len; j++) {
		vertices[v_index]->n[j]->update_geo();
		vertices[v_index]->n[j]->update_energy();
		H += vertices[v_index]->n[j]->H;
	}
	vertices[v_index]->update_energy();

	vertices[v_index]->make_last();
	double n_x = vertices[v_index]->n_vec.x;
	double n_y = vertices[v_index]->n_vec.y;
	double n_z = vertices[v_index]->n_vec.z;
	double l1_x = -n_y;
	double l1_y = n_x;
	double l1_z = 0;

	double H_new;

	for (double move = -0.02e-6; move < 0.02e-6; move += 0.001e-6) {
		vertices[v_index]->point->x = vertices[v_index]->point_last->x + n_x*move;
		vertices[v_index]->point->y = vertices[v_index]->point_last->y + n_y*move;
		vertices[v_index]->point->z = vertices[v_index]->point_last->z + n_z*move;
		vertices[v_index]->update_geo();
		vertices[v_index]->update_energy();
		H_new = 0;
		H_new = vertices[v_index]->H;
		for (int j = 0; j < len; j++) {
			vertices[v_index]->n[j]->update_geo();
			H_new += vertices[v_index]->n[j]->H;
		}
		//std::cout << "move: " << move << "\tdH: " << H_new - H << std::endl;
		nfp << move << '\t' << H_new - H << std::endl;
	}

	for (double move = -0.02e-6; move < 0.02e-6; move += 0.001e-6) {
		vertices[v_index]->point->x = vertices[v_index]->point_last->x + l1_x*move;
		vertices[v_index]->point->y = vertices[v_index]->point_last->y + l1_y*move;
		vertices[v_index]->point->z = vertices[v_index]->point_last->z + l1_z*move;
		vertices[v_index]->update_geo();
		vertices[v_index]->update_energy();
		H_new = 0;
		H_new = vertices[v_index]->H;
		for (int j = 0; j < len; j++) {
			vertices[v_index]->n[j]->update_geo();
			H_new += vertices[v_index]->n[j]->H;
		}
		//std::cout << "move: " << move << "\tdH: " << H_new - H << std::endl;
		lfp1 << move << '\t' << H_new - H << std::endl;
	}

	nfp.close();
	lfp1.close();
}
