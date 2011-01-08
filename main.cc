#include <fstream>
#include "BasisSet.hh"

int main()
{
	std::ifstream is("basis/STO-3G.molpro");
	if (!is.good())
	{
		std::cerr << "Failed to open file\n";
		return 1;
	}

	BasisSet set;
	set.read<BasisSet::Molpro>(is);
	std::cout << set << "\n";

	return 0;
}
