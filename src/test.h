#pragma once

#include<string>
#include<vector>
#include<map>

namespace test {

	class TestCase {
	public:
		TestCase(const std::string& n_name, void(*n_func)());
		std::string name;
		bool run_test();

		// Test flow control
		void new_step(const std::string& step_name);

		// Assertions
		inline bool assert_bool(bool x, const std::string& failure_info = std::string()) {
			assertion_result(x, failure_info);
			return x;
		}
		inline bool assert_bools_lv(bool good, bool pass, const std::string& warning_info = std::string(), const std::string& failure_info = std::string()) {
			if (pass)
				assertion_result(good, warning_info, false);
			else
				assertion_result(pass, failure_info);
			return pass;
		}

	private:
		void(*test_content)();

		std::string name_block;
		bool test_success;
		int test_step;
		std::map<int, bool> step_success;

		// Test flow control
		void test_start();
		void test_end();
		
		// Pass and fail
		void step_result(int which_step);
		void assertion_result(bool res, const std::string& failure_info = std::string(), bool make_error = true);
	};

	// Test cases registering
	std::vector<TestCase*>& get_test_cases();
	void run_all_tests();

}