#include <fstream>
#include "PeriodicTable.hh"
#include "CGTOShellList.hh"
#include "Geometry.hh"
#include "io/manipulators.hh"
#include "BasisSet.hh"
#include "XYZMatrix.hh"
#include "ZMatrix.hh"
#include "CGTO.hh"
#include "Basis.hh"
#include "Dispatcher.hh"
#include "boys.hh"
#include "HartreeFock.hh"

int main()
{
	try
	{
		PeriodicTable table("data/elements.dat");
		CGTOShellList shells;
		IndentingOStream os(std::cout);

		std::string job("C 0 0 0\nH 1 1 1\nH -1 -1 1\nH -1 1 -1\nH 1 -1 -1");
		//std::string job("H 0.7 0 0\nH -0.7 0 0");
		std::istringstream iss(job);
		JobIStream jis(iss);

		Geometry geom;
		jis >> geom;
		os << geom << "\n";
		//geom.toPrincipalAxes();
		//os << geom << "\n";

		BasisSet set;
		//std::ifstream is("basis_sets/STO-3G.molcas");
		std::ifstream is("basis_sets/6-31G**.turbomole");
		if (!is.good())
			throw Li::Exception("Failed to open basis set");
		//set.scan<BasisSet::Molcas>(is);
		set.scan<BasisSet::Turbomole>(is);
		os << set << "\n";

		Basis basis;
		set.expand(geom, &basis);
		os << basis << "\n";

		std::cout << CGTOShellList::singleton().nrShells() << " shells\n";
		std::cout << CGTOShellList::singleton().nrPairs() << " shell pairs\n";
		std::cout << CGTOShellList::singleton().nrQuads() << " shell quads\n";
		std::cout << basis.size() << " basis functions\n";
		std::cout << Dispatcher::singleton().nrPairs() << " bf pairs\n";
		std::cout << Dispatcher::singleton().nrQuads() << " bf quartets\n";
		os << basis.overlap() << "\n\n" << basis.kineticEnergy() << "\n\n"
			<< basis.nuclearAttraction(geom.positions(), geom.charges())
			<< "\n";

		HartreeFock hf(basis, geom, 1, false);
		std::cout << "Total energy: " << hf.energy() << "\n";
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
/*
		Eigen::Vector3d weights, widths;
		weights << 0.15432897, 0.53532814, 0.44463454;
		widths << 6.36242139, 1.15892300, 0.31364979;

		CGTO<0,0,0> bf(weights, widths, Eigen::Vector3d::Zero());
		IndentingOStream os(std::cout);
		os << bf << "\n";
		os << bf.eval(Eigen::Vector3d::Ones());
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
