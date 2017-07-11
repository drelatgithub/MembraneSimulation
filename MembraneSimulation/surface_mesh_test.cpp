#define _USE_MATH_DEFINES
#include"surface_mesh.h"

using namespace math_public;

test::TestCase MS::vertex::test_case("Vertex Test", []() {
	test_case.new_step("Initializing");
	LOG(TEST_DEBUG) << "Generating a hexagonal mesh which includes 7 vertices...";
	std::vector<vertex*> vertices;

	// We are only interested in the vertex in the center.
	vertices.push_back(new vertex(new Vec3(0, 0, 0)));

	for (int i = 0; i < 6; i++) {
		vertices.push_back(new vertex(new Vec3(cos(M_PI / 3 * i), sin(M_PI / 3 * i), 0)));
		vertices[0]->n.push_back(vertices[i + 1]);
		vertices[i + 1]->n.push_back(vertices[0]);
		vertices[0]->dump_data_vectors();
		vertices[i + 1]->dump_data_vectors();
	}
	vertices[0]->gen_next_prev_n();

	int N = vertices.size();

	test_case.new_step("Check neighbor counts");
	test_case.assert_bool(N - 1 == vertices[0]->neighbors);

	vertices[0]->update_geo();
	test_case.new_step("Check angles, area and curvature");
	for (int i = 0; i < 6; i++) {
		test_case.assert_bool(equal(vertices[0]->theta[i], M_PI / 3), "Theta is incorrect.");
		test_case.assert_bool(equal(vertices[0]->theta2[i], M_PI / 3), "Theta2 is incorrect.");
		test_case.assert_bool(equal(vertices[0]->theta3[i], M_PI / 3), "Theta3 is incorrect.");
	}
	test_case.assert_bool(equal(vertices[0]->area, sqrt(3) / 2.0), "Area is incorrect.");
	test_case.assert_bool(equal(vertices[0]->curv_h, 0), "Mean curvature is incorrect.");
	test_case.assert_bool(equal(vertices[0]->curv_g, 0), "Gaussian curvature is incorrect.");

	test_case.new_step("Check self derivatives");
	vertices[0]->point->set(0, 0, 0.5);
	vertices[0]->update_geo();
	vertices[0]->make_last();
	Vec3 dx(-0.001, 0.00001, 0.0005);
	double cur_area = vertices[0]->area,
		cur_curv_h = vertices[0]->curv_h,
		cur_curv_g = vertices[0]->curv_g;
	double ex_del_area = vertices[0]->d_area.dot(dx),
		ex_del_curv_h = vertices[0]->d_curv_h.dot(dx),
		ex_del_curv_g = vertices[0]->d_curv_g.dot(dx);

	*(vertices[0]->point) += dx; // Move the central vertex by a little bit.
	vertices[0]->update_geo();

	double del_area = vertices[0]->area - cur_area,
		del_curv_h = vertices[0]->curv_h - cur_curv_h,
		del_curv_g = vertices[0]->curv_g - cur_curv_g;
	double log_diff_del_area = log10(fabs(ex_del_area - del_area)),
		log_diff_del_curv_h = log10(fabs(ex_del_curv_h - del_curv_h)),
		log_diff_del_curv_g = log10(fabs(ex_del_curv_g - del_curv_g));

	LOG(TEST_DEBUG) << "Area changed: " << del_area << " Expected: " << ex_del_area << " lg diff: " << log_diff_del_area;
	test_case.assert_bools_lv(log_diff_del_area <= -6, log_diff_del_area <= -5, "Area change inaccurate", "Area change incorrect");

	LOG(TEST_DEBUG) << "Mean curv changed: " << del_curv_h << " Expected: " << ex_del_curv_h << " lg diff: " << log_diff_del_curv_h;
	test_case.assert_bools_lv(log_diff_del_curv_h <= -6, log_diff_del_curv_h <= -5, "Mean Curv change inaccurate", "Mean Curv change incorrect");

	LOG(TEST_DEBUG) << "Gaussian curv changed: " << del_curv_g << " Expected: " << ex_del_curv_g << " lg diff: " << log_diff_del_curv_g;
	test_case.assert_bools_lv(log_diff_del_curv_g <= -6, log_diff_del_curv_g <= -5, "Gaussian Curv change inaccurate", "Gaussian Curv change incorrect");

	test_case.new_step("Cleaning");
	for (int i = 0; i < N; i++) {
		vertices[i]->release_point();
		delete vertices[i];
	}
});
