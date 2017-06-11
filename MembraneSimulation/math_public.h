#pragma once

namespace math_public {

	inline int loop_add(const int lhs, const int rhs, const int loop_size) {
		// The answer is in {0, 1, ..., loop_size - 1}.
		// Either operand (lhs or rhs) could be negative.
		// loop_size must be positive.
		int raw = (lhs + rhs) % loop_size;
		if (raw < 0)raw += loop_size;
		return raw;
	}
}