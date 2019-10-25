/*
 * BlockReader.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomic Recombination Lab.
 *      Email : abysslover@gmail.com
 *    License : GPLv3
 */

#include "BlockReader.hpp"

namespace pgr {
BlockReader::BlockReader() : total_file_size(-1), n_cores(1) {
}

BlockReader::~BlockReader() {
}

void BlockReader::set_n_cores(int64_t n_cores) {
	this->n_cores = n_cores;
}
void BlockReader::set_output_suffix(const std::string& output_prefix) {
	this->output_prefix = output_prefix;
}
void BlockReader::set_known_position_path(const std::string& a_known_position_path) {
	TimeChecker checker;
	checker.setTarget("BlockReader.set_known_position_path");
	checker.start();
	string line;
	ifstream in(a_known_position_path, ios::binary);
	vector<string> fields;
	const char* local_delims = "\t ";
	while(getline(in, line, '\n')) {
		StringUtils::c_string_multi_split(line, local_delims, fields);
		string a_combined_key = (boost::format("%s_%s") % fields[0] % fields[1]).str();
		known_positions.insert(a_combined_key);
	}
	cout << checker;
}
void BlockReader::set_path(const std::string& a_path) {
	this->file_path = a_path;
	string copied_path(a_path);
	boost::replace_all(copied_path, ".vcf", output_prefix);
	this->output_file_path = copied_path;
	cout << (boost::format("[BlockReader.set_path] INPUT: %s\n") % file_path).str();
	cout << (boost::format("[BlockReader.set_path] OUTPUT: %s\n") % output_file_path).str();
//	boost::filesystem::path p(a_path);
//	this->working_dir = p.parent_path().string();
	calculate_skip_points();
}
void BlockReader::set_delimters(const std::string& a_delim) {
	this->delims = a_delim;
}

void BlockReader::calculate_skip_points() {
	TimeChecker checker;
	checker.setTarget("BlockReader.calculate_skip_points");
	checker.start();
	this->total_file_size = IOUtils::get_file_size(file_path);
	double size_of_block = IOUtils::MIDDLE_BLOCK_SIZE;
	int64_t n_blocks  = total_file_size / (double) size_of_block;
	if(n_blocks < n_cores) {
		size_of_block = total_file_size / (double) n_cores;
		n_blocks = total_file_size / (double) size_of_block;
	}
	skip_points.resize(n_blocks + 1);
	chunk_sizes.resize(n_blocks);
	IOUtils::find_skip_points(skip_points, file_path,
				size_of_block, total_file_size, n_blocks, n_cores);
	for (uint64_t skip_id = 0; skip_id < skip_points.size() - 1; ++skip_id) {
		chunk_sizes[skip_id] = skip_points[skip_id + 1] - skip_points[skip_id];
//		cout << (boost::format("[BlockReader.calculate_skip_points] %d: %d~%d\n") % skip_id % skip_points[skip_id] % skip_points[skip_id + 1]).str();
	}
	cout << "[BlockReader.calculate_skip_points] # skip points: " << skip_points.size() << "\n";
	cout << checker;
}

void BlockReader::set_callback(callback_type read_callback) {
	this->read_callback = read_callback;
}

void BlockReader::read_lines() {
	TimeChecker checker;
	checker.setTarget("BlockReader.read_lines");
	checker.start();
	vector<boost::function<void()> > tasks;
	boost::mutex print_mutex;
	uint64_t max_block_index = chunk_sizes.size();
//	vector<int64_t> progress_vector;
//	progress_vector.resize(max_block_index);
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
			tasks.push_back([&, block_id] {
				const uint64_t skip_bytes = skip_points[block_id];
				const uint64_t chunk_bytes = chunk_sizes[block_id];
				string local_output_file_path = (boost::format("%s.%s") % output_file_path % block_id).str();
				ofstream out(local_output_file_path, ios::binary);
				ifstream input_read_stream(file_path, ios::binary);
				input_read_stream.seekg(skip_bytes, ios::beg);
				string line;
				uint64_t local_processed_bytes = 0;
				vector<string> fields;
				const char * local_delims = delims.c_str();
				while(getline(input_read_stream, line, '\n')) {
					local_processed_bytes += line.size() + 1;
					read_callback(block_id, print_mutex, known_positions, line, local_delims, fields, out);
					if(local_processed_bytes >= chunk_bytes) {
						break;
					}
				}
//				progress_vector[block_id] = local_processed_bytes;
//				{
//					boost::lock_guard<boost::mutex> a_lock(print_mutex);
//					int64_t current_processed = 0;
//					for(auto processed_bytes: progress_vector) {
//						current_processed += processed_bytes;
//					}
//					double processed_percent = current_processed / (double)total_file_size * 100.0;
//					cout << (boost::format("%.2f %% processed.\n") % processed_percent).str();
//				}
			});
	}
	ParallelRunner::run_unbalanced_load(n_cores, tasks);
	cout << checker;
}

void BlockReader::read_serial_mode() {
	TimeChecker checker;
	checker.setTarget("BlockReader.read_serial_mode");
	checker.start();
	vector<boost::function<void()> > tasks;
	boost::mutex print_mutex;
	uint64_t max_block_index = chunk_sizes.size();
//	vector<int64_t> progress_vector;
//	progress_vector.resize(max_block_index);
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		cout << (boost::format("[BlockReader.read_serial_mode] block: %s\n") % block_id).str();
		const uint64_t skip_bytes = skip_points[block_id];
		const uint64_t chunk_bytes = chunk_sizes[block_id];
		string local_output_file_path = (boost::format("%s.%s") % output_file_path % block_id).str();
		ofstream out(local_output_file_path, ios::binary);
		ifstream input_read_stream(file_path, ios::binary);
		input_read_stream.seekg(skip_bytes, ios::beg);
		string line;
		uint64_t local_processed_bytes = 0;
		vector<string> fields;
		const char * local_delims = delims.c_str();
		while(getline(input_read_stream, line, '\n')) {
			local_processed_bytes += line.size() + 1;
			read_callback(block_id, print_mutex, known_positions, line, local_delims, fields, out);
			if(local_processed_bytes >= chunk_bytes) {
				break;
			}
		}
	}
	cout << checker;
}
void BlockReader::combine_results() {
	uint64_t max_block_index = chunk_sizes.size();
	ofstream out(output_file_path, ios::binary);
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		string local_output_file_path = (boost::format("%s.%s") % output_file_path % block_id).str();
		out << IOUtils::read_fully(local_output_file_path);
	}
	remove_temporary_files(output_file_path);
}
void BlockReader::remove_temporary_files(const string& a_path_prefix) {
	vector<boost::function<void()> > tasks;
	uint64_t max_block_index = chunk_sizes.size();
	for (uint64_t block_id = 0; block_id < max_block_index; ++block_id) {
		tasks.push_back([&, block_id] {
			string local_output_file_path = (boost::format("%s.%s") % output_file_path % block_id).str();
			boost::filesystem::remove(local_output_file_path);
		});
	}
	ParallelRunner::run_unbalanced_load(n_cores, tasks);
}

}
