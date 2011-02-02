#include <fstream>
#include "PeriodicTable.hh"
#include "Geometry.hh"
#include "IndentingOStream.hh"
#include "BasisSet.hh"

int main()
{
	try
	{
		PeriodicTable table("data/elements.dat");
/*
		std::string job("C 0 0 0\nH 1 1 1\nH -1 -1 1\nH -1 1 -1\nH 1 -1 -1");
		std::istringstream iss(job);
		JobIStream jis(iss);

		Geometry geom;
		jis >> geom;

		IndentingOStream os(std::cout);
		os << geom << "\n";
*/
		BasisSet set;
		std::ifstream is("basis_sets/STO-3G.molpro");
		if (!is.good())
			throw Li::Exception("Failed to open basis set");
		set.read<BasisSet::Molpro>(is);
		IndentingOStream os(std::cout);
		os << set << "\n";
	}
	catch (const Li::Exception& e)
	{
		std::cerr << "Exception: " << e.what()
			<< "\nBacktrace:\n" << e.backtrace()
			<< "\n";
		return 1;
	}
	catch (const std::exception& e)
	{
		std::cerr << "Exception: " << e.what() << "\n";
		return 1;
	}

	return 0;
}
