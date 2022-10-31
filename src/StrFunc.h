
#ifndef STRFUNC_H_
#define STRFUNC_H_

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

extern const int na_col;

namespace StrFunc
{
	int find_token_in_header(std::string token_name, const std::vector<std::string>& header, std::string file);
	// split string 
	std::vector<std::string> split(const std::string &s, char delim);
	std::vector<std::string> get_token_names(std::string file);
}

#endif

