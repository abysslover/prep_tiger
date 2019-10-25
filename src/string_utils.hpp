/*
 * string_utils.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: Eun-Cheon Lim @ Postech Plant Genomics Lab.
 */

#ifndef STRING_UTILS_HPP_
#define STRING_UTILS_HPP_

#include <string>
using namespace std;
#include <bitset>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

template<typename C>
void split_string(string const& s, char const* d, C& ret) {
	C output;

	bitset<255> delims;
	while (*d) {
		unsigned char code = *d++;
		delims[code] = true;
	}
	typedef string::const_iterator iter;
	iter beg;
	bool in_token = false;
	for (string::const_iterator it = s.begin(), end = s.end(); it != end;
			++it) {
		if (delims[*it]) {
			if (in_token) {
				output.push_back(typename C::value_type(beg, it));
				in_token = false;
			}
		} else if (!in_token) {
			beg = it;
			in_token = true;
		}
	}
	if (in_token)
		output.push_back(typename C::value_type(beg, s.end()));
	output.swap(ret);
}
typedef string::const_iterator iter;
typedef boost::iterator_range<iter> string_view;
#endif /* STRING_UTILS_HPP_ */
