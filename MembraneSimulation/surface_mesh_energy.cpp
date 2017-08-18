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
	d_H += surface_repulsion_k*(d_I[3] * f.S);

	for (int i = 0; i < 3; i++) {
		f.v[i]->inc_d_H_int(surface_repulsion_k*(f.d_S[i] * I + d_I[i] * f.S));
	}


}
void MS::filament_tip::calc_repulsion(MS::surface_mesh& sm) {
	H = 0;
	d_H.set(0, 0, 0);
	int n_f = sm.facets.size();
	for (int i = 0; i < n_f; i++) {
		calc_repulsion_facet(*(sm.facets[i]));
	}
}


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
