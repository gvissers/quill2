#include <fstream>
#include "PeriodicTable.hh"
#include "IndentingOStream.hh"
#include "FilteringIStream.hh"
#include "JobFilter.hh"

int main()
{
	int i, j, k;
	std::string line("1 , 2,3,,4");
	std::cout << "line = \"" << line << "\"\n";

	std::istringstream is;
	is.clear();
	is.str(line);
	i = j = k = -1;
	is >> i >> j >> k;
	std::cout << "i = " << i << ", j = " << j << ", k = " << k << "\n";

	is.clear();
	is.str(line);
	FilteringIStream<JobFilter> fis(is);
	i = j = k = -1;
	fis >> i >> j >> k;
	std::cout << "i = " << i << ", j = " << j << ", k = " << k << "\n";

	try
	{
		PeriodicTable table("data/elements.dat");
		IndentingOStream os(std::cout);

		os << PeriodicTable::getSingleton() << "\n\n"
			<< table["Cl"] << "\n"
			<< table.findByNumber(12) << "\n"
			<< table.findByName("soDIUM") << "\n";
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
