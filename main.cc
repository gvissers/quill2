#include <fstream>
#include "PeriodicTable.hh"
#include "Geometry.hh"
#include "IndentingOStream.hh"
#include "BasisSet.hh"
#include "XYZMatrix.hh"
#include "ZMatrix.hh"

int main()
{
	try
	{
		PeriodicTable table("data/elements.dat");

		std::string job("C 0 0 0\nH 1 1 1\nH -1 -1 1\nH -1 1 -1\nH 1 -1 -1");
		std::istringstream iss(job);
		JobIStream jis(iss);

		IndentingOStream os(std::cout);

		Geometry geom;
		jis >> geom;

		os << geom << "\n";
/*
		BasisSet set;
		std::ifstream is("basis_sets/STO-3G.molpro");
		if (!is.good())
			throw Li::Exception("Failed to open basis set");
		set.read<BasisSet::Molpro>(is);
		IndentingOStream os(std::cout);
		os << set << "\n";
*/
/*
		std::string mat = "C\n"
			"H   1 1.089000\n"
			"H   1 1.089000  2  109.4710\n"
			"H   1 1.089000  2  109.4710  3  120.0000\n"
			"H   1 1.089000  2  109.4710  3 -120.0000\n";
*/
/*
		std::string mat = "C\n"
			"O     1  1.4\n"
			"H     2  0.95      1  109.471\n"
			"H     1  1.089      2  109.471     3   0\n"
			"H     1  1.089      2  109.471     3   120\n"
			"H     1  1.089      2  109.471     3  -120";
		std::istringstream iss(mat);
		JobIStream jis(iss);
		Geometry geom;
		jis >> geom;

		IndentingOStream os(std::cout);
		os << geom << "\n";
*/
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
