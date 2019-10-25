/*
 * clean_vcf.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "dependency.hpp"
using namespace std;
using namespace pgr;

void inline read_by_line(int64_t block_id, boost::mutex& print_mutex, const std::set<string>& known_positions, const string& line, const char* delims, vector<string>& fields, ofstream& out) {
	if ('#' == line[0]) {
		return;
	}
	const char* info_delims = ";";
	const char* dp4_delims = ",";
	vector<string> info_fields;
	vector<string> dp4_fields;

	StringUtils::c_string_multi_split(line, delims, fields);
	StringUtils::c_string_multi_split(fields[7], info_delims, info_fields);
	boost::replace_all(fields[0], "Chr", "");
	if ("C" == fields[0] || "M" == fields[0] || "Pt" == fields[0] || "Mt" == fields[0]) {
		return;
	}
	string combined_key = (boost::format("%s_%s") % fields[0] % fields[1]).str();
	if (known_positions.end() == known_positions.find(combined_key)) {
		return;
	}
//	if ("1" != fields[0]) {
//		return;
//	}
	string found_DP = "";
	string found_DP4 = "";
	string found_indel = ".";
	string formatted_line;

	int64_t sum_dp4_first = 0;
	int64_t sum_dp4_second = 0;
//	bool debug = "1" == fields[0] && "9573760" == fields[1];
//	if (debug) {
//		boost::lock_guard<boost::mutex> a_lock(print_mutex);
//		cout << line << "\n";
//	}
	for (auto& a_field : info_fields) {
		if("INDEL" == a_field) {
			found_indel = "1";
		} else if(string::npos != a_field.find("DP=")) {
			boost::replace_all(a_field, "DP=", "");
			found_DP = a_field;
//			if (debug) {
//				boost::lock_guard<boost::mutex> a_lock(print_mutex);
//				cout << found_DP << "\n";
//			}
		} else if(string::npos != a_field.find("DP4=")) {
			boost::replace_all(a_field, "DP4=", "");
			found_DP4 = a_field;
			StringUtils::c_string_multi_split(found_DP4, dp4_delims, dp4_fields);
//			if (debug) {
//				boost::lock_guard<boost::mutex> a_lock(print_mutex);
//				cout << "(" << dp4_fields[0] << "," << dp4_fields[1] << "," << dp4_fields[2] << "," << dp4_fields[3] << ")\n";
//			}
			sum_dp4_first = boost::lexical_cast<int64_t>(dp4_fields[0]) + boost::lexical_cast<int64_t>(dp4_fields[1]);
			sum_dp4_second = boost::lexical_cast<int64_t>(dp4_fields[2]) + boost::lexical_cast<int64_t>(dp4_fields[3]);
//				boost::replace_all(found_DP4, ",", "\t");
		}
	}
//		formatted_line = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") %
//				fields[0] % fields[1] % fields[3] % fields[4] % fields[5] %
//				found_indel % found_DP % found_DP4).str();
	formatted_line = (boost::format("%s\t%s\t%s\t%s\t%s\t%s\n") %
					fields[0] % fields[1] % fields[3] % sum_dp4_first % fields[4] % sum_dp4_second).str();
	out << formatted_line;
}

int main(int argc, char **argv) {
	setbuf(stdout, NULL);
	castle::IOUtils::liftrlimit();
	castle::OptionParser opt(argc, argv);
	if(!opt.is_complete()) {
		opt.show_help();
		return 0;
	}
	castle::TimeChecker checker;
	checker.setTarget("main");
	checker.start();
	string input_file_path = opt.file_list_path;
	BlockReader br;
	br.set_n_cores(checker.get_number_of_cores());
	br.set_known_position_path(opt.variant_path);
	br.set_output_suffix(opt.output_suffix);
	br.set_path(input_file_path);
	br.set_callback(read_by_line);
	br.set_delimters("\t");
//	br.read_serial_mode();
	br.read_lines();
	br.combine_results();
	cout << checker;
}
