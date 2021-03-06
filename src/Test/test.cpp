#include"test.h"
#include"log.h"

using namespace test;

test::TestCase::TestCase(const std::string& n_name, void(*n_func)()) :name(n_name), test_content(n_func) {
	get_test_cases().push_back(this);
	name_block.append("[").append(name).append("] ");
}
bool test::TestCase::run_test() {
	test_start();
	test_content();
	test_end();
	return test_success;
}
void test::TestCase::test_start() {
	test_success = true;
	test_step = 0;
	LOG(TEST_INFO) << name_block << "======== Test case \"" << name << "\" started. ========";
}
void test::TestCase::test_end() {
	if (test_step > 0) {
		step_result(test_step);
	}
	LOG((test_success ? TEST_INFO : TEST_ERROR)) << name_block << "======== Test case " << (test_success ? "passed" : "failed") << ". ========";
}
void test::TestCase::new_step(const std::string& step_name) {
	if (test_step > 0) {
		step_result(test_step);
	}
	++test_step;
	LOG(TEST_STEP) << name_block << "Step " << test_step << ": " << step_name;
	step_success[test_step] = true;
}

void test::TestCase::step_result(int which_step) {
	LOG((step_success[which_step] ? TEST_INFO : TEST_ERROR)) << name_block << "Step " << which_step << " " << (step_success[which_step] ? "passed" : "failed") << ".";
}
void test::TestCase::assertion_result(bool res, const std::string& failure_info, bool make_error) {
	if (!res) {
		if (make_error) {
			if (test_step > 0 && step_success[test_step])step_success[test_step] = false;
			if (test_success)test_success = false;

			if (test_step > 0)
				LOG(TEST_ERROR) << name_block << "Assertion failed in step " << test_step << ". " << failure_info;
			else
				LOG(TEST_ERROR) << name_block << "Assertion failed. " << failure_info;
		}
		else {
			LOG(TEST_WARNING) << name_block << failure_info;
		}

	}
}


std::vector<TestCase*>& test::get_test_cases() {
	static std::vector<TestCase*>* ptr = new std::vector<TestCase*>;
	return *ptr;
}

void test::run_all_tests() {
	int num_test_cases = 0, num_passed_test_cases = 0;
	int N = get_test_cases().size();
	for (int i = 0; i < N; i++) {
		if (get_test_cases()[i]->run_test())num_passed_test_cases++;
		num_test_cases++;
	}
	if (num_test_cases == num_passed_test_cases) {
		LOG(TEST_INFO) << "All tests passed.";
	}
	else {
		LOG(TEST_ERROR) << num_test_cases - num_passed_test_cases << " test(s) out of " << num_test_cases << " failed.";
	}
}