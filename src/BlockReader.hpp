/*
 * BlockReader.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#ifndef BLOCKREADER_HPP_
#define BLOCKREADER_HPP_

#include "dependency.hpp"
using namespace castle;
// function(block_id, current_line)
typedef boost::function<void(int64_t, boost::mutex&, const std::set<string>&, const string&, const char*, vector<string>&, ofstream& out) > callback_type;
namespace pgr {
	class BlockReader {
	private:
		std::string file_path;
		std::set<string> known_positions;
		std::string output_prefix;
		std::string output_file_path;
//		std::string working_dir;
		std::vector<uint64_t> skip_points;
		std::vector<uint64_t> chunk_sizes;
		int64_t total_file_size;
		int64_t n_cores;
		callback_type read_callback;
		string delims;
	public:
		BlockReader();
		~BlockReader();
		void set_n_cores(int64_t n_cores);
		void set_output_suffix(const std::string& output_prefix);
		void set_known_position_path(const std::string& a_known_position_path);
		void set_path(const std::string& a_path);
		void set_delimters(const std::string& a_delim);
		void calculate_skip_points();
		void set_callback(callback_type read_callback);
		void read_lines();
		void read_serial_mode();
		void combine_results();
		void remove_temporary_files(const string& a_path_prefix);
	};
}
#endif /* BLOCKREADER_HPP_ */
