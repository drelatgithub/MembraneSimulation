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

#ifdef _WIN32
#  include<Windows.h>
#endif

namespace logger {
	enum Level {
		Debug = 1 << 0,
		Info = 1 << 1,
		Warning = 1 << 2,
		Error = 1 << 3,
		TestDebug = 1 << 4,
		TestInfo = 1 << 5,
		TestStep = 1 << 6,
		TestWarning = 1 << 7,
		TestError = 1 << 8
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

	inline void change_color(Level level, bool restore=false) {
#ifdef _WIN32
		static HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);

		if (restore) {
			SetConsoleTextAttribute(hConsole, 7); // White on Black
		}
		else {
			switch (level) {
			case Debug: case TestDebug: SetConsoleTextAttribute(hConsole, 8); break; // Gray on Black
			case Info: SetConsoleTextAttribute(hConsole, 7); break; // White on Black
			case TestInfo: SetConsoleTextAttribute(hConsole, 10); break; // Light Green on Black
			case Warning: case TestWarning: SetConsoleTextAttribute(hConsole, 14); break; // Light Yellow on Black
			case Error: case TestError: SetConsoleTextAttribute(hConsole, 12); break; // Light Red on Black
			case TestStep: SetConsoleTextAttribute(hConsole, 11); break; // Light Aqua on Black
			}
		}
#else
		// do nothing
#endif
	}
	inline void restore_color() {
		change_color(Debug, true);
	}

}


// User interface
#define DEBUG logger::Debug
#define INFO logger::Info
#define WARNING logger::Warning
#define ERROR logger::Error
#define TEST_DEBUG logger::TestDebug
#define TEST_INFO logger::TestInfo
#define TEST_STEP logger::TestStep
#define TEST_WARNING logger::TestWarning
#define TEST_ERROR logger::TestError

#define LOG(LEVEL) logger::Writer(__FILE__, __LINE__, LOG_FUNC, LEVEL)