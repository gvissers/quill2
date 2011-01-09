#include <fstream>
#include "BasisSet.hh"

int main()
{
	std::ifstream is("basis/STO-3G.dalton");
	if (!is.good())
	{
		std::cerr << "Failed to open file\n";
		return 1;
	}

	try
	{
		BasisSet set;
		set.read<BasisSet::Dalton>(is);
		std::cout << set << "\n";
	}
	catch (const BasisSet::ParseError& e)
	{
		std::cerr << "Parse error on line " << e.line << ": "
			<< e.what() << "\nBacktrace:\n" << e.backtrace()
			<< "\n";
		return 1;
	}

	return 0;
}
