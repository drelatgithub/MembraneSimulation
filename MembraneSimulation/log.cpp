#include<chrono>
#include<ctime>
#include<iomanip>

#include"log.h"
using namespace logger;


Writer::Writer(const char* run_file, const int run_line, const char* run_func, Level new_level) :
	log_buffer(), level(new_level),
	run_file_name(run_file), run_line_name(run_line), run_func_name(run_func) {

	Logger::build_writer(*this, new_level);
}
Writer::~Writer() {
	log_dispatch();
}
void Writer::log_dispatch() {
	if (FILE_OUT) {
		(*file_o) << file_prefix << log_buffer.str() << file_suffix;
		file_o->flush();
	}
	if (SCN_OUT) {
		change_color(level);
		std::cout << scn_prefix << log_buffer.str() << scn_suffix;
		restore_color();
	}
}


std::map<Level, std::string> Logger::level_literal;
std::ofstream Logger::file_o;
std::map<Disp_lv, int> Logger::file_lv, Logger::scn_lv;

void Logger::uni_init() {
	level_literal[Debug] = std::string("DEBUG");
	level_literal[Info] = std::string("INFO");
	level_literal[Warning] = std::string("WARNING");
	level_literal[Error] = std::string("ERROR");
	level_literal[TestDebug] = std::string("TEST DEBUG");
	level_literal[TestInfo] = std::string("TEST INFO");
	level_literal[TestStep] = std::string("TEST STEP");
	level_literal[TestWarning] = std::string("TEST WARNING");
	level_literal[TestError] = std::string("TEST ERROR");
}

void Logger::default_init(const char *file_name) {
	file_o.open(file_name);

	// The default settings are as follows
	int	all = Debug | Info | Warning | Error,
		test_all = TestDebug | TestInfo | TestStep | TestWarning | TestError;
	int min_warning = Warning | Error,
		min_info = Info | Warning | Error,
		test_min_info = TestInfo | TestStep | TestWarning | TestError,
		not_info = Info ^ all;

	file_lv[On_lv] = all | test_all;
	scn_lv[On_lv] = min_info | test_min_info;

	file_lv[Date_lv] = all | test_all;
	file_lv[File_lv] = all | test_all;
	file_lv[Line_lv] = all | test_all;
	file_lv[Function_lv] = all;
	scn_lv[Date_lv] = all | test_all;
	scn_lv[File_lv] = not_info;
	scn_lv[Line_lv] = min_warning;
	scn_lv[Function_lv] = all;

	uni_init();

}

void Logger::build_writer(Writer& writer, Level level) {
	writer.FILE_OUT = (level & file_lv[On_lv]);
	writer.SCN_OUT = (level & scn_lv[On_lv]);
	if (writer.FILE_OUT) {
		writer.file_o = &file_o;
	}
	
	info_gen(writer, writer.file_prefix, writer.file_suffix, file_lv, level);
	info_gen(writer, writer.scn_prefix, writer.scn_suffix, scn_lv, level);

}

void Logger::info_gen(const Writer& writer, std::string &prefix, std::string &suffix, std::map<Disp_lv, int> &disp_setting, Level level) {
	std::stringstream prefix_buffer, suffix_buffer;

	if (disp_setting[Date_lv] & level)
		prefix_buffer << time_gen() << " ";
	prefix_buffer << "[" << level_literal[level] << "] ";
	if (disp_setting[File_lv] & level)
		prefix_buffer << "[File " << writer.run_file_name << "] ";
	if (disp_setting[Line_lv] & level)
		prefix_buffer << "[Line " << writer.run_line_name << "] ";
	if (disp_setting[Function_lv] & level)
		prefix_buffer << "[Function " << writer.run_func_name << "] ";
	prefix = prefix_buffer.str();

	suffix_buffer << std::endl;
	suffix = suffix_buffer.str();
}
std::string Logger::time_gen() {
	using namespace std::chrono;

	system_clock::time_point p = system_clock::now();
	milliseconds ms = duration_cast<milliseconds>(p.time_since_epoch());
	seconds s = duration_cast<seconds>(ms);

	std::time_t time_to_sec = s.count();
	tm timeinfo_to_sec;
	localtime_s(&timeinfo_to_sec, &time_to_sec);
	std::size_t ms_remain = ms.count() % 1000;

	std::stringstream ss;
	ss << timeinfo_to_sec.tm_year + 1900 << '-'
		<< std::setfill('0') << std::setw(2) << timeinfo_to_sec.tm_mon + 1 << '-'
		<< std::setw(2) << timeinfo_to_sec.tm_mday << ' '
		<< std::setw(2) << timeinfo_to_sec.tm_hour << ':'
		<< std::setw(2) << timeinfo_to_sec.tm_min << ':'
		<< std::setw(2) << timeinfo_to_sec.tm_sec << '.'
		<< std::setw(3) << ms_remain;

	return ss.str();
	
}