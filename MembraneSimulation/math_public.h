#pragma once

#include<math.h>

namespace math_public {

	inline int loop_add(const int lhs, const int rhs, const int loop_size) {
		// The answer is in {0, 1, ..., loop_size - 1}.
		// Either operand (lhs or rhs) could be negative.
		// loop_size must be positive.
		int raw = (lhs + rhs) % loop_size;
		if (raw < 0)raw += loop_size;
		return raw;
	}

	// define class for 3d vector
	class Vec3 {
	public:
		double x, y, z;
		Vec3() :x(0), y(0), z(0) { update(); }
		Vec3(double nx, double ny, double nz) :x(nx), y(ny), z(nz) { update(); }
		Vec3(const Vec3 &another) :x(another.x), y(another.y), z(another.z) { update(); }

		// vector plus, minus and multiplication
		inline Vec3 operator+(const Vec3 &operand)const {
			return Vec3(x + operand.x, y + operand.y, z + operand.z);
		}
		inline const Vec3& operator+=(const Vec3 &operand) {
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
		inline const Vec3& operator-=(const Vec3 &operand) {
			x -= operand.x; y -= operand.y; z -= operand.z;
			update();
			return *this;
		}
		inline Vec3 operator*(const double operand)const {
			return Vec3(x*operand, y*operand, z*operand);
		}
		inline const Vec3& operator*=(const double operand) {
			x *= operand; y *= operand; z *= operand;
			update();
			return *this;
		}
		friend inline Vec3 operator*(const double op1, const Vec3 &op2);
		inline Vec3 operator/(const double operand)const {
			return Vec3(x / operand, y / operand, z / operand);
		}
		inline const Vec3& operator/=(const double operand) {
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

		// norms
		double norm, norm2;
		inline double get_norm() {
			norm2 = x*x + y*y + z*z;
			norm = sqrt(norm2);
			return norm;
		}

		// update properties
		inline void update() {
			get_norm();
		}

	};
	inline Vec3 operator*(const double op1, const Vec3 &op2) {
		return op2 * op1;
	}

	inline double dist2(const Vec3 &op1, const Vec3 &op2) {
		return (op1 - op2).norm2;
	}
	inline double dist(const Vec3 &op1, const Vec3 &op2) {
		return (op1 - op2).norm;
	}
}