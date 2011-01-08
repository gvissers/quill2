#include <locale>
#include <algorithm>
#include <tr1/functional>
#include "support.hh"

using namespace std::tr1::placeholders; 

std::string ucFirst(const std::string& str)
{
	if (str.empty()) return std::string();

	std::locale loc;
	std::string res(1, std::toupper<char>(str[0], loc));
	std::transform(str.begin()+1, str.end(), std::back_inserter(res),
		std::tr1::bind(std::tolower<char>, _1, loc));

	return res;
}
