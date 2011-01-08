#include <cctype>
#include "support.hh"

std::string ucFirst(const std::string& str)
{
	if (str.empty()) return std::string();

	std::string res(1, std::toupper(str[0]));
	std::transform(str.begin()+1, str.end(), std::back_inserter(res),
		std::tolower);

	return res;
}
