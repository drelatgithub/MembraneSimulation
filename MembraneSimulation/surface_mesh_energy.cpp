#define _USE_MATH_DEFINES

#include<math.h>

#include"common.h"
#include"surface_mesh_tip.h"
#include"surface_mesh.h"

using namespace math_public;


const double k_c = 1e-19; // Bending modulus
const double k_g = -2 * k_c; // Saddle-splay modulus
const double c_0 = 0.0; // Spontaneous curvature
const double gamma = 0.4; // Surface tension

// Use power n=4 for surface repulsive potential
// @ h0, H = k * h^(-2)
const double surface_repulsion_en_0 = 6.9e-20; // k_B * 5000K
const double surface_repulsion_h0 = 5e-9;
const double surface_repulsion_k = surface_repulsion_en_0 * surface_repulsion_h0 * surface_repulsion_h0;


void MS::vertex::clear_energy() {
	H = 0;
	d_H.set(0, 0, 0);
}

void MS::vertex::calc_H_area() {
	H_area = gamma / 2 / area0 * (area - area0) * (area - area0);
	d_H_area = gamma / area0 * (area - area0) * d_area;
	for each(vertex* each_n in n) {
		d_H_area += gamma / each_n->area0 * (each_n->area - each_n->area0) * each_n->dn_area[each_n->neighbor_indices_map[this]];
	}
}
void MS::vertex::calc_H_curv_h() {
	H_curv_h = 2 * k_c*(curv_h - c_0)*(curv_h - c_0) * area;
	d_H_curv_h = 4 * k_c * (curv_h - c_0) * d_curv_h * area + 2 * k_c * (curv_h - c_0) * (curv_h - c_0) * d_area;
	for each(vertex* each_n in n) {
		int i = each_n->neighbor_indices_map[this];
		d_H_curv_h += 4 * k_c * (each_n->curv_h - c_0) * each_n->dn_curv_h[i] * each_n->area + 2 * k_c * (each_n->curv_h - c_0) * (each_n->curv_h - c_0) * each_n->dn_area[i];
	}
}
void MS::vertex::calc_H_curv_g() {
	H_curv_g = k_g * curv_g * area;
	d_H_curv_g = k_g * (d_curv_g * area + curv_g * d_area);
	for each(vertex* each_n in n) {
		int i = each_n->neighbor_indices_map[this];
		d_H_curv_g += k_g * (each_n->dn_curv_g[i] * each_n->area + each_n->curv_g * each_n->dn_area[i]);
	}
}
void MS::vertex::inc_d_H_int(const Vec3 &d) {
	d_H_int += d;
	d_H += d;
}


void MS::vertex::update_energy() {
	calc_H_area();
	calc_H_curv_h();
	//calc_H_curv_g();
	calc_H_int();
	sum_energy();
}


void MS::filament_tip::calc_repulsion_facet(MS::facet& f) {
	// Calculate interaction energy between the filament tip and a certain facet

	if (is_in_a_plane(*(f.v[0]->point), *(f.v[1]->point), *(f.v[2]->point), *point)) {
		*point -= f.n_vec*1e-10; // Move into the cell a little bit.
	}
	Vec3 r01 = *(f.v[1]->point) - *(f.v[0]->point), r12 = *(f.v[2]->point) - *(f.v[1]->point), rp0 = *(f.v[0]->point) - *point;

	// See notebook starting page 75
	double A = r01.get_norm2();
	double B = r12.get_norm2();
	double C = rp0.get_norm2();
	double D = 2 * dot(r01, rp0);
	double E = 2 * dot(r12, rp0);
	double F = 2 * dot(r01, r12);

	// Derivative subscript is 0,1,2,p
	Vec3 d_A[4] = { -2 * r01, 2 * r01, Vec3(), Vec3() };
	Vec3 d_B[4] = { Vec3(), -2 * r12, 2 * r12, Vec3() };
	Vec3 d_C[4] = { 2 * rp0, Vec3(), Vec3(), -2 * rp0 };
	Vec3 d_D[4] = { 2 * (-rp0 + r01), 2 * rp0, Vec3(), -2 * r01 };
	Vec3 d_E[4] = { 2 * r12, -2 * rp0, 2 * rp0, -2 * r12 };
	Vec3 d_F[4] = { -2 * r12, 2 * (r12 - r01), 2 * r01, Vec3() };

	double A1 = 2 * A*E - D*F;
	double A2 = 2 * B*D - 2 * A*E + (D - E)*F;
	double A3 = -4 * A*B - 2 * B*D + F*(E + F);

	Vec3 d_A1[4], d_A2[4], d_A3[4];
	for (int i = 0; i < 4; i++) {
		d_A1[i] = 2 * (d_A[i] * E + A*d_E[i]) - (d_D[i] * F + D*d_F[i]);
		d_A2[i] = 2 * (d_B[i] * D + B*d_D[i]) - 2 * (d_A[i] * E + A*d_E[i]) + ((d_D[i] - d_E[i])*F + (D - E)*d_F[i]);
		d_A3[i] = -4 * (d_A[i] * B + A*d_B[i]) - 2 * (d_B[i] * D + B*d_D[i]) + (d_F[i] * (E + F) + F*(d_E[i] + d_F[i]));
	}

	double B1 = 4 * A*C - D*D;
	double B2 = 4 * A*C - D*D + 4 * B*C - E*E + 4 * C*F - 2 * D*E;
	double B3 = 4 * B*A - F*F + 4 * B*C - E*E + 4 * B*D - 2 * E*F;
	double BB1 = sqrt(B1);
	double BB2 = sqrt(B2);
	double BB3 = sqrt(B3);

	Vec3 d_B1[4], d_B2[4], d_B3[4], d_BB1[4], d_BB2[4], d_BB3[4];
	for (int i = 0; i < 4; i++) {
		d_B1[i] = 4 * (d_A[i] * C + A*d_C[i]) - 2 * D*d_D[i];
		d_B2[i] = 4 * (d_A[i] * C + A*d_C[i]) - 2 * D*d_D[i] + 4 * (d_B[i] * C + B*d_C[i]) - 2 * E*d_E[i] + 4 * (d_C[i] * F + C*d_F[i]) - 2 * (d_D[i] * E + D*d_E[i]);
		d_B3[i] = 4 * (d_B[i] * A + B*d_A[i]) - 2 * F*d_F[i] + 4 * (d_B[i] * C + B*d_C[i]) - 2 * E*d_E[i] + 4 * (d_B[i] * D + B*d_D[i]) - 2 * (d_E[i] * F + E*d_F[i]);
		d_BB1[i] = d_B1[i] / 2 / BB1;
		d_BB2[i] = d_B2[i] / 2 / BB2;
		d_BB3[i] = d_B3[i] / 2 / BB3;
	}

	double C1 = 2 * A + D;
	double C2 = 2 * A + D + E + 2 * (B + F);
	double C3 = 2 * B + E + F;
	double D1 = D;
	double D2 = D + E;
	double D3 = E + F;
	Vec3 d_C1[4], d_C2[4], d_C3[4], d_D1[4], d_D2[4], d_D3[4];
	for (int i = 0; i < 4; i++) {
		d_C1[i] = 2 * d_A[i] + d_D[i];
		d_C2[i] = 2 * d_A[i] + d_D[i] + d_E[i] + 2 * (d_B[i] + d_F[i]);
		d_C3[i] = 2 * d_B[i] + d_E[i] + d_F[i];
		d_D1[i] = d_D[i];
		d_D2[i] = d_D[i] + d_E[i];
		d_D3[i] = d_E[i] + d_F[i];
	}

	double E1 = atan(C1 / BB1);
	double E2 = atan(C2 / BB2);
	double E3 = atan(C3 / BB3);
	double F1 = atan(D1 / BB1);
	double F2 = atan(D2 / BB2);
	double F3 = atan(D3 / BB3);
	Vec3 d_E1[4], d_E2[4], d_E3[4], d_F1[4], d_F2[4], d_F3[4];
	for (int i = 0; i < 4; i++) {
		d_E1[i] = 1 / (1 + (C1 / BB1)*(C1 / BB1)) * (BB1*d_C1[i] - C1*d_BB1[i]) / B1;
		d_E2[i] = 1 / (1 + (C2 / BB2)*(C2 / BB2)) * (BB2*d_C2[i] - C2*d_BB2[i]) / B2;
		d_E3[i] = 1 / (1 + (C3 / BB3)*(C3 / BB3)) * (BB3*d_C3[i] - C3*d_BB3[i]) / B3;
		d_F1[i] = 1 / (1 + (D1 / BB1)*(D1 / BB1)) * (BB1*d_D1[i] - D1*d_BB1[i]) / B1;
		d_F2[i] = 1 / (1 + (D2 / BB2)*(D2 / BB2)) * (BB2*d_D2[i] - D2*d_BB2[i]) / B2;
		d_F3[i] = 1 / (1 + (D3 / BB3)*(D3 / BB3)) * (BB3*d_D3[i] - D3*d_BB3[i]) / B3;
	}

	double G1 = A1 / BB1;
	double G2 = A2 / BB2;
	double G3 = A3 / BB3;
	Vec3 d_G1[4], d_G2[4], d_G3[4];
	for (int i = 0; i < 4; i++) {
		d_G1[i] = (BB1*d_A1[i] - A1*d_BB1[i]) / B1;
		d_G2[i] = (BB2*d_A2[i] - A2*d_BB2[i]) / B2;
		d_G3[i] = (BB3*d_A3[i] - A3*d_BB3[i]) / B3;
	}

	double I_numerator = G1*(E1 - F1) + G2*(E2 - F2) + G3*(E3 - F3);
	double I_denominator = B*D*D + A*(-4 * B*C + E*E) + F*(-D*E + C*F);
	Vec3 d_I_numerator[4], d_I_denominator[4];
	for (int i = 0; i < 4; i++) {
		d_I_numerator[i] = d_G1[i] * (E1 - F1) + G1*(d_E1[i] - d_F1[i]) + d_G2[i] * (E2 - F2) + G2*(d_E2[i] - d_F2[i]) + d_G3[i] * (E3 - F3) + G3*(d_E3[i] - d_F3[i]);
		d_I_denominator[i] = d_B[i] * D*D + B * 2 * D*d_D[i] + d_A[i] * (-4 * B*C + E*E) + A*(-4 * (d_B[i] * C + B*d_C[i]) + 2 * E*d_E[i]) + d_F[i] * (-D*E + C*F) + F*(-(d_D[i] * E + D*d_E[i]) + d_C[i] * F + C*d_F[i]);
	}

	double I = I_numerator / I_denominator;
	Vec3 d_I[4];
	for (int i = 0; i < 4; i++) {
		d_I[i] = (I_denominator*d_I_numerator[i] - I_numerator*d_I_denominator[i]) / (I_denominator*I_denominator);
	}

	// Calculate energy
	// H = k * |r01 x r12| * I
	double en = surface_repulsion_k * f.S * I;
	H += en;

	for (int i = 0; i < 3; i++) {
		f.v[i]->inc_d_H_int(surface_repulsion_k*(f.d_S[i] * I + d_I[i] * f.S));
	}


}
void MS::filament_tip::calc_repulsion(MS::surface_mesh& sm) {
	H = 0;
	int n_f = sm.facets.size();
	for (int i = 0; i < n_f; i++) {
		calc_repulsion_facet(*(sm.facets[i]));
	}
}


// The following function calculating point-facet interaction energy using facet integral is archived and currently not used.
//void MS::facet::inc_H_int(math_public::Vec3 *p) {
//	/*
//	Purpose:
//		This function calculates the interaction energy between this facet and a certain point (filament tip).
//		The energy would be added to both the interaction energy (H_int) and the total energy (H) of this facet.
//		The derivative of this energy should go to the point and each vertex of this facet.
//
//	Method:
//		We first calculate the integral of interaction energy over the whole plane,
//		and then use the cut-off method to get approximately how much portion of that integral
//		lies in the facet.
//
//	Limits:
//		In order that this function works, one has to make sure that the triangle has positive area,
//		i.e. the 3 vertices are not in a line (which also implies that no 2 vertices could be at the same position).
//	*/
//
//	// The total integral over the whole plane is
//	// H = d0^n * d^(2-n) * f(n) * k
//	// where k is the coefficient, d0 is the distance of force, f(n) is a constexpr polynomial of n.
//	static int po_pwr = 6; // negative power of distance in expression of potential energy
//	static double r_0_factor = 1.0 / sqrt(po_pwr-1);
//	static double d0 = 1e-9;
//	static double energy_coe = 1e-1; // in J/m^2
//	static double d0_pwr = pow(d0, po_pwr);
//	static double pre_calc = d0_pwr * energy_coe;
//
//	/**********************************
//	find the foot of perpendicular O on the triangle
//	rO = r0 + alpha v1 + beta v2
//	**********************************/
//	Vec3 v1 = *(v[1]->point) - *(v[0]->point), v2 = *(v[2]->point) - *(v[0]->point);
//	Vec3 r12 = v2 - v1;
//	Mat3 d0_v1 = -Eye3, d1_v1 = Eye3, d0_v2 = -Eye3, d2_v2 = Eye3, d1_r12 = -Eye3, d2_r12 = Eye3;
//
//	Vec3 r0p = *p - *(v[0]->point);
//	Mat3 d0_r0p = -Eye3;
//
//	// Calculate the double area of the triangle.
//	double S2 = cross(v1, v2).get_norm();
//	if (S2 <= 0) {
//		LOG(WARNING) << "Facet area is not positive. S2 = " << S2;
//	}
//
//	// alpha and beta need to satisfy the perpendicular condition
//	// A * (alpha, beta)' = B
//	// So (alpha, beta)' = A^(-1) * B
//	double dot12 = dot(v1, v2);
//	Vec3 d0_dot12 = -v2 - v1, d1_dot12 = v2, d2_dot12 = v1; // Already taken into account those "Eye"-derivatives.
//
//	v1.calc_norm();
//	Vec3 d0_norm2_v1 = -v1 * 2, d1_norm2_v1 = v1 * 2, d0_norm_v1 = -v1 / v1.norm, d1_norm_v1 = v1 / v1.norm;
//	v2.calc_norm();
//	Vec3 d0_norm2_v2 = -v2 * 2, d2_norm2_v2 = v2 * 2, d0_norm_v2 = -v2 / v2.norm, d2_norm_v2 = v2 / v2.norm;
//	r12.calc_norm();
//	Vec3 d1_norm2_r12 = -r12 * 2, d2_norm2_r12 = r12 * 2, d1_norm_r12 = -r12 / r12.norm, d2_norm_r12 = r12 / r12.norm;
//
//	Vec3 d0_S2 = (-v1.norm2*v2 - v2.norm2*v1 + dot12*(v1 + v2)) / S2,
//		d1_S2 = (v2.norm2*v1 - dot12*v2) / S2,
//		d2_S2 = (v1.norm2*v2 - dot12*v1) / S2;
//
//	// det(A) = |v1|^2 |v2|^2 - (v1 * v2)^2, but theoretically this is essentially S2^2
//	double det_A = S2*S2;
//	double det_A2 = det_A*det_A;
//	Vec3 d0_det_A = 2 * S2*d0_S2,
//		d1_det_A = 2 * S2*d1_S2,
//		d2_det_A = 2 * S2*d2_S2;
//	double AR11 = v2.norm2 / det_A,
//		AR12 = -dot12 / det_A, // AR21 = AR12
//		AR22 = v1.norm2 / det_A;
//	Vec3 d0_AR11 = (det_A*d0_norm2_v2 - v2.norm2*d0_det_A) / det_A2,
//		d1_AR11 = -v2.norm2*d1_det_A / det_A2,
//		d2_AR11 = (det_A*d2_norm2_v2 - v2.norm2*d2_det_A) / det_A2,
//		d0_AR12 = -(det_A*d0_dot12 - dot12*d0_det_A) / det_A2,
//		d1_AR12 = -(det_A*d1_dot12 - dot12*d1_det_A) / det_A2,
//		d2_AR12 = -(det_A*d2_dot12 - dot12*d2_det_A) / det_A2,
//		d0_AR22 = (det_A*d0_norm2_v1 - v1.norm2*d0_det_A) / det_A2,
//		d1_AR22 = (det_A*d1_norm2_v1 - v1.norm2*d1_det_A) / det_A2,
//		d2_AR22 = -v1.norm2*d2_det_A / det_A2;
//	double B1 = r0p.dot(v1), B2 = r0p.dot(v2);
//	Vec3 d0_B1 = -v1 - r0p, d1_B1 = r0p,
//		d0_B2 = -v2 - r0p, d2_B2 = r0p;
//
//	double alpha = AR11*B1 + AR12*B2,
//		beta = AR12*B1 + AR22*B2; // Actually it's AR21 * B1 + AR22 * B2
//	Vec3 d0_alpha = d0_AR11*B1 + AR11*d0_B1 + d0_AR12*B2 + AR12*d0_B2,
//		d1_alpha = d1_AR11*B1 + AR11*d1_B1 + d1_AR12*B2,
//		d2_alpha = d2_AR11*B1 + d2_AR12*B2 + AR12*d2_B2,
//		d0_beta = d0_AR12*B1 + AR12*d0_B1 + d0_AR22*B2 + AR22*d0_B2,
//		d1_beta = d1_AR12*B1 + AR12*d1_B1 + d1_AR22*B2,
//		d2_beta = d2_AR12*B1 + d2_AR22*B2 + AR22*d2_B2;
//
//	Vec3 r0O = alpha*v1 + beta*v2;
//	Mat3 d0_r0O = d0_alpha.tensor(v1) + alpha*d0_v1 + d0_beta.tensor(v2) + beta*d0_v2,
//		d1_r0O = d1_alpha.tensor(v1) + alpha*d1_v1 + d1_beta.tensor(v2),
//		d2_r0O = d2_alpha.tensor(v1) + d2_beta.tensor(v2) + beta*d2_v2;
//	Vec3 rO = r0O + *(v[0]->point),
//		r1O = r0O - v1,
//		r2O = r0O - v2;
//	Mat3 d0_rO = d0_r0O + Eye3, d1_rO = d1_r0O, d2_rO = d2_r0O,
//		d0_r1O = d0_r0O - d0_v1, d1_r1O = d1_r0O - d1_v1, d2_r1O = d2_r0O,
//		d0_r2O = d0_r0O - d0_v2, d1_r2O = d1_r0O, d2_r2O = d2_r0O - d2_v2;
//	r0O.calc_norm();
//	r1O.calc_norm();
//	r2O.calc_norm();
//	Vec3 rOp = *p - rO; // derivatives w.r.t. r0, r1, r2 are just negative of d0_rO, d1_rO and d2_rO
//	if (false) { // Check whether they are indeed perpendicular
//		LOG(DEBUG) << "rOp dot v1 = " << rOp.dot(v1);
//		LOG(DEBUG) << "rOp dot v2 = " << rOp.dot(v2);
//	}
//
//	/**********************************
//	distance from the point to the triangle
//	**********************************/
//	double d = rOp.get_norm();
//	Vec3 d0_d = (-d0_rO)*rOp / d,
//		d1_d = (-d1_rO)*rOp / d,
//		d2_d = (-d2_rO)*rOp / d;
//	if (false) {
//		LOG(DEBUG) << "Distance to the plane: " << d;
//	}
//
//	// If the point is too far away from the triangle then cut it off
//	double c_f = v1.norm + v2.norm + r12.norm; // circumference of the facet
//	if (d > 5 * d0 || (r0O.norm > c_f && r1O.norm > c_f && r2O.norm > c_f)) {
//		H_int = 0;
//		return; // No increase in energy derivative
//	}
//
//	/**********************************
//	distance from point O to all 3 edges
//	**********************************/
//	double a1 = beta*S2;
//	Vec3 d0_a1 = d0_beta*S2 + beta*d0_S2,
//		d1_a1 = d1_beta*S2 + beta*d1_S2,
//		d2_a1 = d2_beta*S2 + beta*d2_S2;
//	double d1 = a1 / v1.norm;
//	Vec3 d0_d1 = (v1.norm*d0_a1 - a1*d0_norm_v1) / v1.norm2,
//		d1_d1 = (v1.norm*d1_a1 - a1*d1_norm_v1) / v1.norm2,
//		d2_d1 = d2_a1 / v1.norm;
//
//	double a2 = alpha*S2;
//	Vec3 d0_a2 = d0_alpha*S2 + alpha*d0_S2,
//		d1_a2 = d1_alpha*S2 + alpha*d1_S2,
//		d2_a2 = d2_alpha*S2 + alpha*d2_S2;
//	double d2 = a2 / v2.norm;
//	Vec3 d0_d2 = (v2.norm*d0_a2 - a2*d0_norm_v2) / v2.norm2,
//		d1_d2 = d1_a2 / v2.norm,
//		d2_d2 = (v2.norm*d2_a2 - a2*d2_norm_v2) / v2.norm2;
//
//	double a3 = (1 - alpha - beta)*S2;
//	Vec3 d0_a3 = (-d0_alpha - d0_beta)*S2 + (1 - alpha - beta)*d0_S2,
//		d1_a3 = (-d1_alpha - d1_beta)*S2 + (1 - alpha - beta)*d1_S2,
//		d2_a3 = (-d2_alpha - d2_beta)*S2 + (1 - alpha - beta)*d2_S2;
//	double d3 = a3 / r12.norm;
//	Vec3 d0_d3 = d0_a3 / r12.norm,
//		d1_d3 = (r12.norm*d1_a3 - a3*d1_norm_r12) / r12.norm2,
//		d2_d3 = (r12.norm*d2_a3 - a3*d2_norm_r12) / r12.norm2;
//
//	if (false) {
//		LOG(DEBUG) << "Distance to the edges: " << d1 << ", " << d2 << ", " << d3;
//	}
//
//	/**********************************
//	find the affecting region and calculate energy
//	**********************************/
//	double sigma = d * r_0_factor;
//	Vec3 d0_sigma = d0_d*r_0_factor,
//		d1_sigma = d1_d*r_0_factor,
//		d2_sigma = d2_d*r_0_factor;
//	if (sigma >= 0.05*(v1.norm + v2.norm + r12.norm))
//		LOG(WARNING) << "Sigma is not significantly smaller than the scale of the facet.";
//
//	double I = pow(d, 2 - po_pwr) * pre_calc; // The energy integral
//	Vec3 d0_I = (2 - po_pwr)*pow(d, 1 - po_pwr)*pre_calc*d0_d,
//		d1_I = (2 - po_pwr)*pow(d, 1 - po_pwr)*pre_calc*d1_d,
//		d2_I = (2 - po_pwr)*pow(d, 1 - po_pwr)*pre_calc*d2_d;
//	double vcoe1 = exp(-r0O.norm2 / (sigma*sigma)),
//		vcoe2 = exp(-r1O.norm2 / (sigma*sigma)),
//		vcoe3 = exp(-r2O.norm2 / (sigma*sigma));
//	Vec3 d0_vcoe1 = vcoe1 * 2 * (-sigma*d0_r0O*r0O + r0O.norm2 * d0_sigma) / (sigma*sigma*sigma),
//		d1_vcoe1 = vcoe1 * 2 * (-sigma*d1_r0O*r0O + r0O.norm2*d1_sigma) / (sigma*sigma*sigma),
//		d2_vcoe1 = vcoe1 * 2 * (-sigma*d2_r0O*r0O + r0O.norm2*d2_sigma) / (sigma*sigma*sigma),
//		d0_vcoe2 = vcoe2 * 2 * (-sigma*d0_r1O*r1O + r1O.norm2*d0_sigma) / (sigma*sigma*sigma),
//		d1_vcoe2 = vcoe2 * 2 * (-sigma*d1_r1O*r1O + r1O.norm2*d1_sigma) / (sigma*sigma*sigma),
//		d2_vcoe2 = vcoe2 * 2 * (-sigma*d2_r1O*r1O + r1O.norm2*d2_sigma) / (sigma*sigma*sigma),
//		d0_vcoe3 = vcoe3 * 2 * (-sigma*d0_r2O*r2O + r2O.norm2*d0_sigma) / (sigma*sigma*sigma),
//		d1_vcoe3 = vcoe3 * 2 * (-sigma*d1_r2O*r2O + r2O.norm2*d1_sigma) / (sigma*sigma*sigma),
//		d2_vcoe3 = vcoe3 * 2 * (-sigma*d2_r2O*r2O + r2O.norm2*d2_sigma) / (sigma*sigma*sigma);
//	double ecoe1 = smooth_sgn(d1 / sigma),
//		ecoe2 = smooth_sgn(d2 / sigma),
//		ecoe3 = smooth_sgn(d3 / sigma);
//	Vec3 d0_ecoe1 = d_smooth_sgn(d1 / sigma)*(sigma*d0_d1 - d1*d0_sigma) / (sigma*sigma),
//		d1_ecoe1 = d_smooth_sgn(d1 / sigma)*(sigma*d1_d1 - d1*d1_sigma) / (sigma*sigma),
//		d2_ecoe1 = d_smooth_sgn(d1 / sigma)*(sigma*d2_d1 - d1*d2_sigma) / (sigma*sigma),
//		d0_ecoe2 = d_smooth_sgn(d2 / sigma)*(sigma*d0_d2 - d2*d0_sigma) / (sigma*sigma),
//		d1_ecoe2 = d_smooth_sgn(d2 / sigma)*(sigma*d1_d2 - d2*d1_sigma) / (sigma*sigma),
//		d2_ecoe2 = d_smooth_sgn(d2 / sigma)*(sigma*d2_d2 - d2*d2_sigma) / (sigma*sigma),
//		d0_ecoe3 = d_smooth_sgn(d3 / sigma)*(sigma*d0_d3 - d3*d0_sigma) / (sigma*sigma),
//		d1_ecoe3 = d_smooth_sgn(d3 / sigma)*(sigma*d1_d3 - d3*d1_sigma) / (sigma*sigma),
//		d2_ecoe3 = d_smooth_sgn(d3 / sigma)*(sigma*d2_d3 - d3*d2_sigma) / (sigma*sigma);
//	double not_at_corner = (1 - vcoe1)*(1 - vcoe2)*(1 - vcoe3),
//		in_triangle = ecoe1*ecoe2*ecoe3;
//	double en_fact = vcoe1*(v[0]->theta[ind[0]] / (2 * M_PI))
//		+ vcoe2*(v[1]->theta[ind[1]] / (2 * M_PI))
//		+ vcoe3*(v[2]->theta[ind[2]] / (2 * M_PI))
//		+ not_at_corner * in_triangle;
//	Vec3 d0_en_fact = d0_vcoe1*(v[0]->theta[ind[0]] / (2 * M_PI)) + vcoe1*(v[0]->d_theta[ind[0]] / (2 * M_PI))
//		+ d0_vcoe2*(v[1]->theta[ind[1]] / (2 * M_PI)) + vcoe2*(v[1]->dnn_theta[ind[1]] / (2 * M_PI))
//		+ d0_vcoe3*(v[2]->theta[ind[2]] / (2 * M_PI)) + vcoe3*(v[2]->dn_theta[ind[2]] / (2 * M_PI))
//		- d0_vcoe1*(1 - vcoe2)*(1 - vcoe3)*in_triangle
//		- d0_vcoe2*(1 - vcoe1)*(1 - vcoe3)*in_triangle
//		- d0_vcoe3*(1 - vcoe1)*(1 - vcoe2)*in_triangle
//		+ not_at_corner*d0_ecoe1*ecoe2*ecoe3
//		+ not_at_corner*ecoe1*d0_ecoe2*ecoe3
//		+ not_at_corner*ecoe1*ecoe2*d0_ecoe3,
//		d1_en_fact = d1_vcoe1*(v[0]->theta[ind[0]] / (2 * M_PI)) + vcoe1*(v[0]->dn_theta[ind[0]] / (2 * M_PI))
//		+ d1_vcoe2*(v[1]->theta[ind[1]] / (2 * M_PI)) + vcoe2*(v[1]->d_theta[ind[1]] / (2 * M_PI))
//		+ d1_vcoe3*(v[2]->theta[ind[2]] / (2 * M_PI)) + vcoe3*(v[2]->dnn_theta[ind[2]] / (2 * M_PI))
//		- d1_vcoe1*(1 - vcoe2)*(1 - vcoe3)*in_triangle
//		- d1_vcoe2*(1 - vcoe1)*(1 - vcoe3)*in_triangle
//		- d1_vcoe3*(1 - vcoe1)*(1 - vcoe2)*in_triangle
//		+ not_at_corner*d1_ecoe1*ecoe2*ecoe3
//		+ not_at_corner*ecoe1*d1_ecoe2*ecoe3
//		+ not_at_corner*ecoe1*ecoe2*d1_ecoe3,
//		d2_en_fact = d2_vcoe1*(v[0]->theta[ind[0]] / (2 * M_PI)) + vcoe1*(v[0]->dnn_theta[ind[0]] / (2 * M_PI))
//		+ d2_vcoe2*(v[1]->theta[ind[1]] / (2 * M_PI)) + vcoe2*(v[1]->dn_theta[ind[1]] / (2 * M_PI))
//		+ d2_vcoe3*(v[2]->theta[ind[2]] / (2 * M_PI)) + vcoe3*(v[2]->d_theta[ind[2]] / (2 * M_PI))
//		- d2_vcoe1*(1 - vcoe2)*(1 - vcoe3)*in_triangle
//		- d2_vcoe2*(1 - vcoe1)*(1 - vcoe3)*in_triangle
//		- d2_vcoe3*(1 - vcoe1)*(1 - vcoe2)*in_triangle
//		+ not_at_corner*d2_ecoe1*ecoe2*ecoe3
//		+ not_at_corner*ecoe1*d2_ecoe2*ecoe3
//		+ not_at_corner*ecoe1*ecoe2*d2_ecoe3;
//
//	double en = en_fact * I;
//	Vec3 d0_en = d0_en_fact * I + en_fact * d0_I,
//		d1_en = d1_en_fact * I + en_fact * d1_I,
//		d2_en = d2_en_fact * I + en_fact * d2_I;
//
//	H_int += en; H += en;
//	v[0]->inc_d_H_int(d0_en);
//	v[1]->inc_d_H_int(d1_en);
//	v[2]->inc_d_H_int(d2_en);
//	
//}

void MS::surface_mesh::update_energy() {
	int N;
	N = vertices.size();
	for (int i = 0; i < N; i++) {
		vertices[i]->update_energy();
	}
}
double MS::surface_mesh::get_sum_of_energy() {
	double res = 0;
	int N;
	N = vertices.size();
	for (int i = 0; i < N; i++) {
		res += vertices[i]->H;
	}
	return res;
}
