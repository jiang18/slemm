
#include "StrFunc.h"

const int na_col = -1;
const char delimiter = ',';

int StrFunc::find_token_in_header(std::string token_name, const std::vector<std::string>& header, std::string file)
{
	int token_col = na_col;
	if(!token_name.empty()) {
		for (unsigned i=1; i<header.size(); ++i) {
			if(header[i].compare(token_name) == 0) {
				token_col = i;
				break;
			}
		}
		if(token_col == na_col) {
			std::cout<<token_name<<" not found in the header of "<<file<<std::endl;
			exit(1);
		}
	}
	return token_col;
}

std::vector<std::string> StrFunc::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
	std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> StrFunc::get_token_names(std::string file) {
	std::ifstream ifs;
	ifs.open(file);
	if(!ifs.is_open()) {
		std::cout<<"Cannot open "<<file<<std::endl;
		exit(1);
	}
	
	std::string line;
	std::getline(ifs, line);
	ifs.close();
	if(line[line.size()-1] == '\r') line.erase(line.size()-1);
	std::vector<std::string> tokens = StrFunc::split(line, delimiter);
	tokens.erase(tokens.begin());
	return tokens;
}

