#pragma once

#include<math.h>

#include"Test/test.h"

namespace math_public {

	// Declarations first due to circular references
	class Vec3;
	class Mat3;

	// Useful arithmetics
	inline int loop_add(const int lhs, const int rhs, const int loop_size) {
		// The answer is in {0, 1, ..., loop_size - 1}.
		// Either operand (lhs or rhs) could be negative.
		// loop_size must be positive.
		int raw = (lhs + rhs) % loop_size;
		if (raw < 0)raw += loop_size;
		return raw;
	}
	bool is_in_a_plane(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3, double eps = 1e-25);
	inline double smooth_sgn(double x) {
		return 0.5 + 0.5*tanh(x);
	}
	inline double d_smooth_sgn(double x) {
		double tanh_x = tanh(x);
		return 0.5*(1 - tanh_x*tanh_x);
	}
	extern test::TestCase test_case_arithmetics;

	inline bool equal(double op1, double op2, double eps = 1e-6) {
		return fabs(op1 - op2) < eps;
	}

	// define class for 3d vector
	class Vec3 {
	public:
		double x, y, z;
		Vec3() :x(0), y(0), z(0) {}
		Vec3(double nx, double ny, double nz) :x(nx), y(ny), z(nz) {}

		inline Vec3& set(double nx, double ny, double nz) { x = nx; y = ny; z = nz; return *this; }

		// comparisons
		inline bool equal_to(const Vec3& operand, double eps_ratio=1e-10)const {
			double eps = eps_ratio * get_norm();
			return equal(x, operand.x, eps) && equal(y, operand.y, eps) && equal(z, operand.z, eps);
		}
		inline bool equal_to_in_norm(const Vec3& operand, double eps_ratio=1e-10)const {
			double eps = eps_ratio * get_norm();
			return equal(norm, operand.norm, eps);
		}

		// vector plus, minus, multiplication and division
		inline Vec3 operator+(const Vec3 &operand)const {
			return Vec3(x + operand.x, y + operand.y, z + operand.z);
		}
		inline Vec3& operator+=(const Vec3 &operand) {
			x += operand.x; y += operand.y; z += operand.z;
			update();
			return *this;
		}
		inline Vec3 operator-()const {
			return Vec3(-x, -y, -z);
		}
		inline Vec3 operator-(const Vec3 &operand)const {
			return Vec3(x - operand.x, y - operand.y, z - operand.z);
		}
		inline Vec3& operator-=(const Vec3 &operand) {
			x -= operand.x; y -= operand.y; z -= operand.z;
			update();
			return *this;
		}
		inline Vec3 operator*(const double operand)const {
			return Vec3(x*operand, y*operand, z*operand);
		}
		inline Vec3& operator*=(const double operand) {
			x *= operand; y *= operand; z *= operand;
			update();
			return *this;
		}
		friend inline Vec3 operator*(const double op1, const Vec3 &op2);
		inline Vec3 operator/(const double operand)const {
			return Vec3(x / operand, y / operand, z / operand);
		}
		inline Vec3& operator/=(const double operand) {
			x /= operand; y /= operand; z /= operand;
			update();
			return *this;
		}

		// dot product and cross product
		inline double dot(const Vec3 &operand)const {
			return x*operand.x + y*operand.y + z*operand.z;
		}
		inline double operator*(const Vec3 &operand)const { // "*" as dot product
			return this->dot(operand);
		}
		inline Vec3 cross(const Vec3 &operand)const {
			return Vec3(y*operand.z - z*operand.y, z*operand.x - x*operand.z, x*operand.y - y*operand.x);
		}
		// tensor product
		Mat3 tensor(const Vec3& operand)const;

		// norms
		mutable double norm, norm2;
		// norms are not initially calculated.
		// To use norm, one has to first calc_norm()
		// For temporary use just use the return value of get_norm() or get_norm2()
		inline void calc_norm()const {
			norm2 = x*x + y*y + z*z;
			norm = sqrt(norm2);
		}
		inline double get_norm()const {
			calc_norm();
			return norm;
		}
		inline double get_norm2()const {
			calc_norm();
			return norm2;
		}

		// Calculate skew-symmetric matrix to be used in a cross product
		// (Vec3) a => (Mat3) [a]_x,
		// where a x b = [a]_x * b
		Mat3 to_skew_cross()const;

		// update properties
		inline void update() {
			calc_norm();
		}

		// string display
		std::string str(int mode = 0)const; // mode 0: numbers separated by '\t', 1: like (x, y, z)

		// test
		static test::TestCase test_case;

	};
	inline Vec3 operator*(const double op1, const Vec3 &op2) {
		return op2 * op1;
	}
	inline double dot(const Vec3 &op1, const Vec3 &op2) {
		return op1.dot(op2);
	}
	inline Vec3 cross(const Vec3 &op1, const Vec3 &op2) {
		return op1.cross(op2);
	}
	inline double triple_product(const Vec3& op1, const Vec3& op2, const Vec3& op3) {
		return op1.cross(op2).dot(op3);
	}

	inline double dist2(const Vec3 &op1, const Vec3 &op2) {
		return (op1 - op2).get_norm2();
	}
	inline double dist(const Vec3 &op1, const Vec3 &op2) {
		return (op1 - op2).get_norm();
	}


	// define class for 3x3 matrix
	class Mat3 {
	public:
		Vec3 x, y, z;
		Mat3() :x(), y(), z() { update(); }
		Mat3(const Vec3& nx, const Vec3& ny, const Vec3& nz) :x(nx), y(ny), z(nz) { update(); }
		Mat3(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33) :
			x(a11, a21, a31), y(a12, a22, a32), z(a13, a23, a33) {
			update();
		}

		// comparisons
		inline bool equal_to(const Mat3& operand, double eps_ratio = 1e-10)const {
			return x.equal_to(operand.x, eps_ratio) && y.equal_to(operand.y, eps_ratio) && z.equal_to(operand.z, eps_ratio);
		}

		// plus, minus, multiplication and division
		inline Mat3 operator+(const Mat3& operand)const {
			return Mat3(x + operand.x, y + operand.y, z + operand.z);
		}
		inline Mat3& operator+=(const Mat3 &operand) {
			x += operand.x; y += operand.y; z += operand.z;
			update();
			return *this;
		}
		inline Mat3 operator-()const {
			return Mat3(-x, -y, -z);
		}
		inline Mat3 operator-(const Mat3& operand)const {
			return Mat3(x - operand.x, y - operand.y, z - operand.z);
		}
		inline Mat3& operator-=(const Mat3 &operand) {
			x -= operand.x; y -= operand.y; z -= operand.z;
			update();
			return *this;
		}
		inline Mat3 operator*(double operand)const {
			return Mat3(x*operand, y*operand, z*operand);
		}
		friend inline Mat3 operator*(double op1, const Mat3& op2);
		inline Mat3& operator*=(const double operand) {
			x *= operand; y *= operand; z *= operand;
			update();
			return *this;
		}
		inline Mat3 operator/(double operand)const {
			return Mat3(x / operand, y / operand, z / operand);
		}
		inline Mat3& operator/=(const double operand) {
			x /= operand; y /= operand; z /= operand;
			update();
			return *this;
		}

		// vector and matrix multiplication
		inline Vec3 operator*(const Vec3& operand)const {
			return Vec3(x_row.dot(operand), y_row.dot(operand), z_row.dot(operand));
		}
		inline Mat3 operator*(const Mat3& operand)const {
			return Mat3(operator*(operand.x), operator*(operand.y), operator*(operand.z));
		}

		Vec3 x_row, y_row, z_row;
		inline void get_row_vec() {
			x_row.set(x.x, y.x, z.x);
			y_row.set(x.y, y.y, z.y);
			z_row.set(x.z, y.z, z.z);
		}

		inline void update() {
			get_row_vec();
		}

		inline Mat3 transpose()const {
			// This function only works when no changes are made after update() is called.
			return Mat3(x_row, y_row, z_row);
		}

		// test
		static test::TestCase test_case;
	};
	inline Mat3 operator*(double op1, const Mat3& op2) {
		return op2*op1;
	}

	const Mat3 Eye3(1, 0, 0, 0, 1, 0, 0, 0, 1);
}