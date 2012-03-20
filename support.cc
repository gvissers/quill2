#include <locale>
#include <algorithm>
#include "support.hh"

std::string lower(const std::string& str)
{
	std::locale loc;
	std::string res;
	std::transform(str.begin(), str.end(), std::back_inserter(res),
		std::bind(std::tolower<char>, _1, loc));

	return res;
}

std::string ucFirst(const std::string& str)
{
	if (str.empty()) return std::string();

	std::locale loc;
	std::string res(1, std::toupper<char>(str[0], loc));
	std::transform(str.begin()+1, str.end(), std::back_inserter(res),
		std::bind(std::tolower<char>, _1, loc));

	return res;
}

std::string& toUCFirst(std::string& str)
{
	if (str.empty())
		return str;

	std::locale loc;
	str[0] = std::toupper<char>(str[0], loc);
	std::transform(str.begin()+1, str.end(), str.begin()+1,
		std::bind(std::tolower<char>, _1, loc));

	return str;
}
