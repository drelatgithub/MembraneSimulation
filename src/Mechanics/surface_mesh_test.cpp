#define _USE_MATH_DEFINES
#include"surface_mesh.h"
#include"surface_mesh_tip.h"

using namespace math_public;

test::TestCase MS::vertex::test_case("Vertex Test", []() {
	test_case.new_step("Initializing");
	LOG(TEST_DEBUG) << "Generating a hexagonal mesh which includes 7 vertices...";
	std::vector<vertex*> vertices;

	// We are only interested in the vertex in the center.
	vertices.push_back(new vertex(new Vec3(0, 0, 0)));

	for (int i = 1; i <= 6; i++) {
		vertices.push_back(new vertex(new Vec3(cos(M_PI / 3 * (i-1)), sin(M_PI / 3 * (i-1)), 0)));
		vertices[0]->n.push_back(vertices[i]);
		vertices[i]->n.push_back(vertices[0]);
		vertices[0]->dump_data_vectors();
		vertices[i]->dump_data_vectors();
		// Fixing neighbor parameters
		vertices[i]->area = vertices[i]->area0 = 1.0; // arbitrary number
		vertices[i]->dn_area[0] = Vec3();
		vertices[i]->dn_curv_h[0] = Vec3();
		vertices[i]->dn_curv_g[0] = Vec3();
		vertices[i]->gen_next_prev_n();
	}
	vertices[0]->gen_next_prev_n();

	for (int i = 0; i < 6; i++) {
		vertices[0]->f.push_back(new facet(vertices[0], vertices[i+1], vertices[loop_add(i, 1, 6)+1]));
		vertices[0]->f[i]->calc_vec();
		vertices[0]->f[i]->calc_normal();
	}

	int N = vertices.size();

	test_case.new_step("Check neighbor counts");
	test_case.assert_bool(N - 1 == vertices[0]->neighbors);

	vertices[0]->calc_angle();
	vertices[0]->calc_area();
	vertices[0]->calc_curv_h();
	vertices[0]->calc_curv_g();
	vertices[0]->calc_normal();
	test_case.new_step("Check angles, area, curvature and normal vector");
	for (int i = 0; i < 6; i++) {
		test_case.assert_bool(equal(vertices[0]->theta[i], M_PI / 3), "Theta is incorrect.");
		test_case.assert_bool(equal(vertices[0]->theta2[i], M_PI / 3), "Theta2 is incorrect.");
		test_case.assert_bool(equal(vertices[0]->theta3[i], M_PI / 3), "Theta3 is incorrect.");
	}
	test_case.assert_bool(equal(vertices[0]->area, sqrt(3) / 2.0), "Area is incorrect.");
	test_case.assert_bool(equal(vertices[0]->curv_h, 0), "Mean curvature is incorrect.");
	test_case.assert_bool(equal(vertices[0]->curv_g, 0), "Gaussian curvature is incorrect.");

	LOG(TEST_DEBUG) << vertices[0]->n_vec.str(1);
	test_case.assert_bool(vertices[0]->n_vec.equal_to(Vec3(0, 0, 1)), "Normal vector is incorrect.");

	test_case.new_step("Check self derivatives");
	double osm_p = 0.1;
	vertices[0]->point->set(0, 0, 0.5);
	for (int i = 0; i < 6; i++) {
		vertices[0]->f[i]->update_geo();
	}
	vertices[0]->calc_angle();
	vertices[0]->calc_area();
	vertices[0]->calc_curv_h();
	vertices[0]->calc_curv_g();
	vertices[0]->calc_normal();
	vertices[0]->calc_volume_op();
	vertices[0]->make_last();
	vertices[0]->make_initial();
	vertices[0]->update_energy(osm_p);
	Vec3 dx(-0.001, 0.00001, 0.0005);
	double cur_area = vertices[0]->area,
		cur_H_area = vertices[0]->H_area,
		cur_curv_h = vertices[0]->curv_h,
		cur_H_curv_h = vertices[0]->H_curv_h,
		cur_curv_g = vertices[0]->curv_g,
		cur_H_curv_g = vertices[0]->H_curv_g,
		cur_volume_op = vertices[0]->volume_op,
		cur_H_osm = vertices[0]->H_osm;
	Vec3 cur_n_vec = vertices[0]->n_vec;
	LOG(TEST_DEBUG) << "-------------------- Before change --------------------";
	LOG(TEST_DEBUG) << "area: " << cur_area << " curv_h: " << cur_curv_h << " curv_g: " << cur_curv_g;
	LOG(TEST_DEBUG) << "H_area: " << cur_H_area << " H_curv_h: " << cur_H_curv_h << " H_curv_g: " << cur_H_curv_g;
	LOG(TEST_DEBUG) << "volume integrand: " << cur_volume_op << " osmotic energy: " << cur_H_osm;
	double ex_del_area = vertices[0]->d_area.dot(dx),
		ex_del_H_area = vertices[0]->d_H_area.dot(dx),
		ex_del_curv_h = vertices[0]->d_curv_h.dot(dx),
		ex_del_H_curv_h = vertices[0]->d_H_curv_h.dot(dx),
		ex_del_curv_g = vertices[0]->d_curv_g.dot(dx),
		ex_del_H_curv_g = vertices[0]->d_H_curv_g.dot(dx),
		ex_del_volume_op = vertices[0]->d_volume_op.dot(dx),
		ex_del_H_osm = vertices[0]->d_H_osm.dot(dx);
	Vec3 ex_del_n_vec = vertices[0]->d_n_vec.transpose()*dx;

	*(vertices[0]->point) += dx; // Move the central vertex by a little bit.
	for (int i = 0; i < 6; i++) {
		vertices[0]->f[i]->update_geo();
	}
	vertices[0]->calc_angle();
	vertices[0]->calc_area();
	vertices[0]->calc_curv_h();
	vertices[0]->calc_curv_g();
	vertices[0]->calc_normal();
	vertices[0]->calc_volume_op();
	vertices[0]->update_energy(osm_p);

	LOG(TEST_DEBUG) << "-------------------- After change --------------------";
	double del_area = vertices[0]->area - cur_area,
		del_H_area = vertices[0]->H_area - cur_H_area,
		del_curv_h = vertices[0]->curv_h - cur_curv_h,
		del_H_curv_h = vertices[0]->H_curv_h - cur_H_curv_h,
		del_curv_g = vertices[0]->curv_g - cur_curv_g,
		del_H_curv_g = vertices[0]->H_curv_g - cur_H_curv_g,
		del_volume_op = vertices[0]->volume_op - cur_volume_op,
		del_H_osm = vertices[0]->H_osm - cur_H_osm;
	Vec3 del_n_vec = vertices[0]->n_vec - cur_n_vec;
	double log_diff_del_area = log10(fabs(ex_del_area - del_area)),
		log_diff_del_H_area = log10(fabs(ex_del_H_area - del_H_area)),
		log_diff_del_curv_h = log10(fabs(ex_del_curv_h - del_curv_h)),
		log_diff_del_H_curv_h = log10(fabs(ex_del_H_curv_h - del_H_curv_h)),
		log_diff_del_curv_g = log10(fabs(ex_del_curv_g - del_curv_g)),
		log_diff_del_H_curv_g = log10(fabs(ex_del_H_curv_g - del_H_curv_g)),
		log_diff_del_volume_op = log10(fabs(ex_del_volume_op - del_volume_op)),
		log_diff_del_H_osm = log10(fabs(ex_del_H_osm - del_H_osm));

	LOG(TEST_DEBUG) << "Area changed: " << del_area << " Expected: " << ex_del_area << " lg diff: " << log_diff_del_area;
	test_case.assert_bools_lv(log_diff_del_area <= -6, log_diff_del_area <= -5, "Area change inaccurate", "Area change incorrect");

	LOG(TEST_DEBUG) << "Area energy changed: " << del_H_area << " Expected: " << ex_del_H_area << " lg diff: " << log_diff_del_H_area;
	test_case.assert_bools_lv(log_diff_del_H_area <= -8, log_diff_del_H_area <= -6, "Area energy change inaccurate", "Area energy change incorrect");

	LOG(TEST_DEBUG) << "Mean curv changed: " << del_curv_h << " Expected: " << ex_del_curv_h << " lg diff: " << log_diff_del_curv_h;
	test_case.assert_bools_lv(log_diff_del_curv_h <= -6, log_diff_del_curv_h <= -5, "Mean Curv change inaccurate", "Mean Curv change incorrect");

	LOG(TEST_DEBUG) << "Mean curv energy changed: " << del_H_curv_h << " Expected: " << ex_del_H_curv_h << " lg diff: " << log_diff_del_H_curv_h;
	test_case.assert_bools_lv(log_diff_del_H_curv_h <= -25, log_diff_del_H_curv_h <= -24, "Mean Curv energy change inaccurate", "Mean Curv energy change incorrect");

	LOG(TEST_DEBUG) << "Gaussian curv changed: " << del_curv_g << " Expected: " << ex_del_curv_g << " lg diff: " << log_diff_del_curv_g;
	test_case.assert_bools_lv(log_diff_del_curv_g <= -6, log_diff_del_curv_g <= -5, "Gaussian Curv change inaccurate", "Gaussian Curv change incorrect");

	LOG(TEST_DEBUG) << "Gaussian curv energy changed: " << del_H_curv_g << " Expected: " << ex_del_H_curv_g << " lg diff: " << log_diff_del_H_curv_g;
	test_case.assert_bools_lv(log_diff_del_H_curv_g <= -25, log_diff_del_H_curv_g <= -24, "Gaussian Curv energy change inaccurate", "Gaussian Curv energy change incorrect");
	
	LOG(TEST_DEBUG) << "Volume integrand changed: " << del_volume_op << " Expected: " << ex_del_volume_op << " lg diff: " << log_diff_del_volume_op;
	test_case.assert_bools_lv(log_diff_del_volume_op <= -7, log_diff_del_volume_op <= -6, "Volume integrand change inaccurate", "Volume integrand change incorrect");

	LOG(TEST_DEBUG) << "Osmotic energy changed: " << del_H_osm << " Expected: " << ex_del_H_osm << " lg diff: " << log_diff_del_H_osm;
	test_case.assert_bools_lv(log_diff_del_H_osm <= -7, log_diff_del_H_osm <= -6, "Osmotic energy change inaccurate", "Osmotic energy change incorrect");

	LOG(TEST_DEBUG) << "Normal vector changed: " << del_n_vec.str(1) << " Expected: " << ex_del_n_vec.str(1);
	test_case.assert_bool(del_n_vec.equal_to(ex_del_n_vec, 1e-3), "Normal vector change incorrect.");

	test_case.new_step("Cleaning");
	for (int i = 0; i < 6; i++) {
		delete vertices[0]->f[i];
	}
	for (int i = 0; i < N; i++) {
		vertices[i]->release_point();
		delete vertices[i];
	}
});

test::TestCase MS::facet::test_case("Facet Test", []() {
	test_case.new_step("Initializing");

	LOG(TEST_DEBUG) << "Generating a triangle mesh which includes 3 vertices...";
	std::vector<vertex*> vertices;
	int N = 3;

	vertices.push_back(new vertex(new Vec3(1e-7, 0, 0)));
	vertices.push_back(new vertex(new Vec3(-0.5e-7, sqrt(3) /2 *1e-7, 0)));
	vertices.push_back(new vertex(new Vec3(-0.5e-7, -sqrt(3) /2 *1e-7, 0)));

	facet f(vertices[0], vertices[1], vertices[2]);

	f.update_geo();

	Vec3 move[3] = {
		Vec3(0.005e-7,0.001e-7,-0.001e-7),
		Vec3(-0.001e-7,0.005e-7,0.001e-7),
		Vec3(0.001e-7,-0.001e-7,0.005e-7)
	};

	test_case.new_step("Check normal vector");
	LOG(TEST_DEBUG) << "Normal vector: " << f.n_vec.str(1);
	test_case.assert_bool(f.n_vec.equal_to(Vec3(0, 0, 1)), "Normal vector incorrect.");

	test_case.new_step("Check area");
	double ex_S = 3 * sqrt(3) / 4 * 1e-14;
	LOG(TEST_DEBUG) << "Area: " << f.S << " Expected: " << ex_S;
	test_case.assert_bool(equal(f.S, ex_S), "Area incorrect.");

	test_case.new_step("Check derivatives");
	Vec3 old_n_vec = f.n_vec;
	double old_S = f.S;
	Vec3 diff_n_vec_ex;
	double diff_S_ex = 0;
	for (int i = 0; i < 3; i++) {
		diff_n_vec_ex += f.d_n_vec[i].transpose() * move[i];
		diff_S_ex += f.d_S[i] * move[i]; // dot product
	}

	for (int i = 0; i < 3; i++) {
		*(f.v[i]->point) += move[i];
	}
	f.update_geo();

	Vec3 diff_n_vec = f.n_vec - old_n_vec;
	double diff_S = f.S - old_S;
	LOG(TEST_DEBUG) << "Normal vector difference: " << diff_n_vec.str(1) << " Expected: " << diff_n_vec_ex.str(1);
	test_case.assert_bool(diff_n_vec.equal_to(diff_n_vec_ex,1e-2), "Normal vector derivative incorrect.");
	test_case.assert_bool(equal(diff_S, diff_S_ex), "Area derivative incorrect.");

	test_case.new_step("Cleaning");
	for (int i = 0; i < N; i++) {
		vertices[i]->release_point();
		delete vertices[i];
	}

});

test::TestCase MS::filament_tip::test_case("Filament Tip Test", []() {
	test_case.new_step("Initializing...");

	LOG(TEST_DEBUG) << "Generating a triangle mesh which includes 3 vertices...";
	std::vector<vertex*> vertices;
	int N = 3;

	vertices.push_back(new vertex(new Vec3(1e-7, 0, 0)));
	vertices.push_back(new vertex(new Vec3(-0.5e-7, sqrt(3) / 2 * 1e-7, 0)));
	vertices.push_back(new vertex(new Vec3(-0.5e-7, -sqrt(3) / 2 * 1e-7, 0)));

	facet f(vertices[0], vertices[1], vertices[2]);

	f.update_geo();

	surface_mesh sm;
	sm.facets.push_back(&f);

	LOG(TEST_DEBUG) << "Generating a filament tip...";
	filament_tip ft(new Vec3(0, 0, 0.5e-8));

	Vec3 move[3] = {
		Vec3(0.005e-7,0.001e-7,-0.001e-7),
		Vec3(-0.001e-7,0.005e-7,0.001e-7),
		Vec3(0.001e-7,-0.001e-7,0.005e-7)
	};
	Vec3 movep(0.001e-7, 0.001e-7, -0.0001e-7);

	test_case.new_step("Calculate energy");
	ft.calc_repulsion(sm);
	LOG(TEST_DEBUG) << "Interaction energy: " << ft.H;
	LOG(TEST_DEBUG) << "Energy derivative on tip: " << ft.d_H.str(1);
	LOG(TEST_DEBUG) << "Energy derivative on vertices: "
		<< f.v[0]->d_H.str(1) << " "
		<< f.v[1]->d_H.str(1) << " "
		<< f.v[2]->d_H.str(1);

	test_case.new_step("Check derivatives");
	double old_H = ft.H;
	double diff_H_ex = 0;
	for (int i = 0; i < 3; i++) {
		diff_H_ex += f.v[i]->d_H * move[i];
	}
	diff_H_ex += ft.d_H * movep;

	for (int i = 0; i < 3; i++) {
		*(f.v[i]->point) += move[i];
	}
	*(ft.point) += movep;
	f.update_geo();
	ft.calc_repulsion(sm);

	double diff_H = ft.H - old_H;
	LOG(TEST_DEBUG) << "Energy difference: " << diff_H << " Expected: " << diff_H_ex;
	test_case.assert_bool(equal(diff_H, diff_H_ex, 1e-20), "Energy derivative incorrect.");

	test_case.new_step("Cleaning");
	for (int i = 0; i < N; i++) {
		vertices[i]->release_point();
		delete vertices[i];
	}

});
