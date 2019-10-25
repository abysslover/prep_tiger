/*
 * OptionParser.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "OptionParser.hpp"

namespace castle {

void OptionParser::show_help() {
	cout << "A preprocessor for TIGER (v.1.0)\n"
			"==============================================================================\n";
	cout << "Syntax:    prep_tiger -f [STR] -v [STR] -s [STR] -h\n";
	cout << "Options:   -f STR  a path of a VCF file to be extracted.\n";
	cout << "           -v STR  a path of file that contains known small variant positions.\n";
	cout << "           -s STR  the output suffix. e.g.) .selected.txt\n";
	cout << "           -h      show this help message.\n\n";
	cout << "Note:      the variant file should be stored in the following manner:\n";
	cout << "           [CHROMOSME]\t[POSITION]\tetc.\n";
	cout << "Contact:   Euncheon Lim <abysslover@gmail.com>\n";
}

OptionParser::OptionParser(int argc, char **argv) {
	string line;
	vector<string> cols;
	// the first line

	for (int i = 1; i < argc; ++i) {
		string argument = argv[i];
		if ("-f" == argument) {
			file_list_path = string(argv[i+1]);
			expand_path(file_list_path);
			working_path = get_working_path(file_list_path) + "/";
			if(!boost::filesystem::exists(file_list_path)) {
				cout << (boost::format("[OptionParser.OptionParser] %s does not exists\n") % file_list_path).str();
				exit(1);
			}
		} else if ("-v" == argument) {
			variant_path = string(argv[i+1]);
			expand_path(variant_path);
			if(!boost::filesystem::exists(variant_path)) {
				cout << (boost::format("[OptionParser.OptionParser] %s does not exists\n") % variant_path).str();
				exit(1);
			}
		} else if ("-s" == argument) {
			output_suffix = string(argv[i+1]);
		} else if ("-h" == argument) {
			show_help();
			break;
		}
	}
}

OptionParser::~OptionParser() {
}

bool OptionParser::is_complete() {
	return !file_list_path.empty() && !variant_path.empty() && !output_suffix.empty();
}
void OptionParser::expand_home_path(string& a_path) {
   if (a_path.empty() || '~' != a_path[0]) {
		   return;
   }
   char const* home = getenv("HOME");
	if (home || ((home = getenv("USERPROFILE")))) {
		a_path.replace(0, 1, home);
	} else {
		char const *hdrive = getenv("HOMEDRIVE"), *hpath = getenv("HOMEPATH");
		a_path.replace(0, 1, string(hdrive) + hpath);
	}
}

string OptionParser::get_current_directory(void) {

	char buff[FILENAME_MAX];
	GetCurrentDir( buff, FILENAME_MAX );
	std::string current_working_dir(buff);
	return current_working_dir;
}

void OptionParser::expand_path(string& a_path) {
	if ('/' != a_path[0] && '~' != a_path[0]) {
		if ('.' == a_path[0] && '/' == a_path[1]) {
			a_path = a_path.substr(2);
		}
		string cur_dir = get_current_directory();
		a_path = cur_dir + "/" + a_path;
	} else if('~' == a_path[0]) {
		expand_home_path(a_path);
	}
}
string OptionParser::get_working_path(const string& a_path) {
	string copied_path(a_path);
	boost::filesystem::path current_path(get_current_directory());
	string file_prefix;
	if ('/' != copied_path[0] && '~' != copied_path[0]) {
		file_prefix = current_path.string();
	} else if('~' == copied_path[0]) {
		expand_home_path(copied_path);
	}
	file_prefix = copied_path.substr(0, copied_path.rfind('/'));
	return file_prefix;
}

string OptionParser::get_working_folder_path() const {
	return working_path;
}

} /* namespace castle */
