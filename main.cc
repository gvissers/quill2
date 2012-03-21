#include <fstream>
#include <iomanip>
#include "PeriodicTable.hh"
#include "CGTOShellList.hh"
#include "Job.hh"
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

void printUsage(const std::string& name)
{
	std::cerr << "Usage: " << name << " <job-file>\n";
}

int main(int argc, const char* argv[])
{
	if (argc != 2)
	{
		printUsage(argv[0]);
		return 1;
	}
	
	try
	{
		PeriodicTable table("data/elements.dat");
		CGTOShellList shells;
		IndentingOStream os(std::cout);

		std::ifstream fin(argv[1]);
		if (!fin.good())
		{
			std::cerr << "Failed to open job file \"" << argv[1] << "\"\n";
			return 1;
		}
		JobIStream jis(fin);
		Job job(jis);
		fin.close();

 		Geometry geom = job.get<Geometry>("geometry");
 		os << geom << "\n";
 		geom.toPrincipalAxes();
 		os << geom << "\n";

		std::string setname = job.get("basis");
		BasisSet set;
		if (!set.findAndScan(setname, "basis_sets"))
		{
			std::cerr << "Failed to find basis set \"" + setname + "\"\n";
			return 1;
		}

		Basis basis;
		set.expand(geom, &basis);
		os << basis << "\n";

		std::cout << CGTOShellList::singleton().nrShells() << " shells\n";
		std::cout << CGTOShellList::singleton().nrPairs() << " shell pairs\n";
		std::cout << CGTOShellList::singleton().nrQuads() << " shell quads\n";
		std::cout << basis.size() << " basis functions\n";
		std::cout << Dispatcher::singleton().nrPairs() << " bf pairs\n";
		std::cout << Dispatcher::singleton().nrQuads() << " bf quartets\n";
		os << "diagonal of S: " << basis.overlap().diagonal() << "\n\n";
		os << basis.overlap() << "\n\n" << basis.kineticEnergy() << "\n\n"
			<< basis.nuclearAttraction(geom.positions(), geom.charges())
			<< "\n";

		HartreeFock hf(basis, geom, 1);
		std::cout << "Total energy: " << std::setprecision(15) << std::fixed
			<< hf.energy() << " Eh\n";
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
