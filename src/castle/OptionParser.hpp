/*
 * OptionParser.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef OPTIONPARSER_HPP_
#define OPTIONPARSER_HPP_
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/cstdint.hpp>
#include "StringUtils.hpp"
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif
#include<iostream>

namespace castle {
using namespace std;
class OptionParser {
public:
	string file_list_path;
	string variant_path;
	string working_path;
	string output_suffix;
public:
	void show_help();
	OptionParser(int argc, char **argv);
	~OptionParser();
	bool is_complete();
	bool is_memory_efficient();
	void expand_path(string& a_path);
	void expand_home_path(string& a_path);
	string get_current_directory();
	uint64_t get_k() const;
	string get_working_path(const string& a_path);
	string get_working_folder_path() const;

};

} /* namespace castle */
#endif /* OPTIONPARSER_HPP_ */
