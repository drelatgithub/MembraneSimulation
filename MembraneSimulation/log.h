#pragma once

#include<iostream>
#include<fstream>
#include<sstream>
#include<map>

#if defined(_MSC_VER)
#  define LOG_COMPILER_MSVC 1
#else
#  define LOG_COMPILER_MSVC 0
#endif

#if LOG_COMPILER_MSVC
#  define LOG_FUNC __FUNCTION__
#elif defined(__func__)
#  define LOG_FUNC __func__
#else
#  define LOG_FUNC ""
#endif

namespace logger {
	enum Level {
		Debug = 1 << 0,
		Info = 1 << 1,
		Warning = 1 << 2,
		Error = 1 << 3
	};

	enum Disp_lv {
		On_lv, Date_lv, File_lv, Line_lv, Function_lv
	};

	class Writer {
		friend class Logger;
	public:
		Writer(const char* run_file, const int run_line, const char* run_func, Level new_level);
		Writer(const Writer& origin) = delete;
		~Writer();

		template <typename T>
		Writer& operator<<(const T& msg) {
			log_buffer << msg;
			return *this;
		}
		Writer& operator<<(std::ostream& (*os)(std::ostream&)) {
			log_buffer << os;
			return *this;
		}
	private:
		const char *run_file_name, *run_func_name;
		const int run_line_name;

		std::stringstream log_buffer;
		std::string file_prefix, scn_prefix, file_suffix, scn_suffix;
		std::ofstream *file_o;

		Level level;
		bool FILE_OUT, SCN_OUT;

		void log_dispatch();
	};

	class Logger {
	public:
		
		static void default_init(const char* file_name);

		static void build_writer(Writer &writer, Level level);
	private:
		static std::map<Level, std::string> level_literal;

		static std::ofstream file_o;
		static std::map<Disp_lv, int> file_lv, scn_lv;

		static void uni_init();

		static void info_gen(const Writer& writer, std::string &prefix, std::string &suffix, std::map<Disp_lv, int> &disp_setting, Level level);
		static std::string time_gen();
	};

}


// User interface
#define DEBUG logger::Debug
#define INFO logger::Info
#define WARNING logger::Warning
#define ERROR logger::Error

#define LOG(LEVEL) logger::Writer(__FILE__, __LINE__, LOG_FUNC, LEVEL)