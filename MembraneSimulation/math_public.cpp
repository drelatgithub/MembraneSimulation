#include"common.h"
#include"math_public.h"

using namespace math_public;

test::TestCase math_public::test_case_arithmetics("Math arithmetics", []() {
	auto& tc = test_case_arithmetics;
	{
		tc.new_step("Check smooth sign function");
		double x[3] = { -40,0,40 },
			dx = 1e-4;
		tc.assert_bool(equal(smooth_sgn(x[0]), 0), "Function incorrect at x -> -infty.");
		tc.assert_bool(equal(smooth_sgn(x[1]), 0.5), "Function incorrect at x = 0.");
		tc.assert_bool(equal(smooth_sgn(x[2]), 1), "Function incorrect at x -> +infty.");
		for (int i = 0; i < 3; i++) {
			double x0 = smooth_sgn(x[i]);
			double ex_del = d_smooth_sgn(x[i]);
			double del = (smooth_sgn(x[i] + dx) - x0) / dx;
			tc.assert_bool(equal(del, ex_del), std::string("Derivative is incorrect at x = ").append(std::to_string(x[i])));
		}
	}
});

Mat3 math_public::Vec3::tensor(const Vec3& operand)const {
	Mat3 ans(operator*(operand.x), operator*(operand.y), operator*(operand.z));
	return ans;
}
Mat3 math_public::Vec3::to_skew_cross()const {
	return Mat3(
		0, -z, y,
		z, 0, -x,
		-y, x, 0
	);
}
std::string math_public::Vec3::str(int mode)const {
	std::stringstream ss;
	switch (mode) {
	case 0: // plain numbers separated by '\t'
		ss << x << '\t' << y << '\t' << z;
		break;
	case 1: // like (x, y, z)
		ss << "(" << x << ", " << y << ", " << z << ")";
		break;
	}
	return ss.str();
}

test::TestCase math_public::Vec3::test_case("Vec3", []() {
	{
		test_case.new_step("Check norm");
		Vec3 op(3, 4, 12);
		double ex = 13;
		test_case.assert_bool(equal(op.get_norm(), ex));
		test_case.assert_bool(equal(op.norm, ex));
	}
	{
		test_case.new_step("Check plus");
		Vec3 op1(1, 2, 3), op2(-3, -2, -1);
		Vec3 ex(-2, 0, 2);
		test_case.assert_bool((op1 + op2).equal_to(ex));
	}
	{
		test_case.new_step("Check minus");
		Vec3 op1(1, 2, 3), op2(-3, -2, -1);
		Vec3 ex(4, 4, 4);
		test_case.assert_bool((op1 - op2).equal_to(ex));
	}
	{
		test_case.new_step("Check multiplication and division");
		Vec3 op1(1, 2, 3); double op2 = -0.5;
		Vec3 ex1(-0.5, -1, -1.5), ex2(-2, -4, -6);
		test_case.assert_bool((op1 * op2).equal_to(ex1));
		test_case.assert_bool((op2 * op1).equal_to(ex1));
		test_case.assert_bool((op1 / op2).equal_to(ex2));
	}
	{
		test_case.new_step("Check compound operators");
		Vec3 op1(1, 2, 3), op2(-3, -2, -1); double op3 = -2;
		Vec3 ex1(-2, 0, 2); double ex2 = sqrt(8);
		Vec3 ex3(1, 2, 3); double ex4 = sqrt(14);
		Vec3 ex5(-2, -4, -6); double ex6 = sqrt(56);
		Vec3 ex7(1, 2, 3); double ex8 = sqrt(14);
		op1 += op2;
		test_case.assert_bool(op1.equal_to(ex1));
		test_case.assert_bool(equal(op1.norm, ex2));
		op1 -= op2;
		test_case.assert_bool(op1.equal_to(ex3));
		test_case.assert_bool(equal(op1.norm, ex4));
		op1 *= op3;
		test_case.assert_bool(op1.equal_to(ex5));
		test_case.assert_bool(equal(op1.norm, ex6));
		op1 /= op3;
		test_case.assert_bool(op1.equal_to(ex7));
		test_case.assert_bool(equal(op1.norm, ex8));
	}
	{
		test_case.new_step("Check dot product");
		Vec3 op1(1, 2, 3), op2(-3, -2, -1);
		double ex = -10;
		test_case.assert_bool(equal(op1*op2, ex));
	}
	{
		test_case.new_step("Check cross product");
		Vec3 op1(1, 2, 3), op2(-3, -2, -1);
		Vec3 ex(4, -8, 4);
		test_case.assert_bool(op1.cross(op2).equal_to(ex));
		test_case.assert_bool(equal(op1*ex, 0));
		test_case.assert_bool(equal(op2*ex, 0));
	}
	{
		test_case.new_step("Check cross product via skew-sym matrix");
		Vec3 op1(1, 2, 3), op2(-3, -2, -1);
		test_case.assert_bool(op1.cross(op2).equal_to(op1.to_skew_cross()*op2));
	}
	{
		test_case.new_step("Check tensor product");
		Vec3 op1(1, 2, 3), op2(-3, -2, -1);
		Mat3 ex(-3, -2, -1, -6, -4, -2, -9, -6, -3);
		test_case.assert_bool(op1.tensor(op2).equal_to(ex));
	}
});


test::TestCase math_public::Mat3::test_case("Mat3", []() {
	{
		test_case.new_step("Check plus and minus");
		Mat3 op1(1, 2, 3, 4, 5, 6, 7, 8, 9);
		Mat3 ex1(2, 2, 3, 4, 6, 6, 7, 8, 10), ex2(0, 2, 3, 4, 4, 6, 7, 8, 8);
		test_case.assert_bool((op1 + Eye3).equal_to(ex1));
		test_case.assert_bool((op1 - Eye3).equal_to(ex2));
	}
	{
		test_case.new_step("Check multiplication and division");
		Mat3 op1(1, 2, 3, 4, 5, 6, 7, 8, 9); double op2 = -0.5;
		Mat3 ex1(-0.5, -1, -1.5, -2, -2.5, -3, -3.5, -4, -4.5), ex2(-2, -4, -6, -8, -10, -12, -14, -16, -18);
		test_case.assert_bool((op1 * op2).equal_to(ex1));
		test_case.assert_bool((op2 * op1).equal_to(ex1));
		test_case.assert_bool((op1 / op2).equal_to(ex2));
	}
	{
		test_case.new_step("Check compound operators");
		Mat3 op1(1, 2, 3, 4, 5, 6, 7, 8, 9), op2(-9, -8, -7, -6, -5, -4, -3, -2, -1); double op3 = -2;
		Mat3 ex1(-8, -6, -4, -2, 0, 2, 4, 6, 8);
		Mat3 ex2(op1);
		Mat3 ex3(-2, -4, -6, -8, -10, -12, -14, -16, -18);
		Mat3 ex4(op1);
		op1 += op2;
		test_case.assert_bool(op1.equal_to(ex1));
		op1 -= op2;
		test_case.assert_bool(op1.equal_to(ex2));
		op1 *= op3;
		test_case.assert_bool(op1.equal_to(ex3));
		op1 /= op3;
		test_case.assert_bool(op1.equal_to(ex4));
	}
	{
		test_case.new_step("Check linear transformation");
		Mat3 op1(1, 2, 3, 4, 5, 6, 7, 8, 9); Vec3 op2(10, 11, 12);
		Vec3 ex(68, 167, 266);
		test_case.assert_bool((op1*op2).equal_to(ex));
	}
	{
		test_case.new_step("Check matrix multiplication");
		Mat3 op1(1, 2, 3, 4, 5, 6, 7, 8, 9), op2(9, 8, 7, 6, 5, 4, 3, 2, 1);
		Mat3 ex(30, 24, 18, 84, 69, 54, 138, 114, 90);
		test_case.assert_bool((op1*op2).equal_to(ex));
	}
});