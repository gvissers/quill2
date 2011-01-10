#include <fstream>
#include "PeriodicTable.hh"
#include "IndentingStreambuf.hh"

int main()
{
	try
	{
		PeriodicTable table("data/elements.dat");
		IndentingStreambuf ind_sb(std::cout.rdbuf());
		std::streambuf *old_sb = std::cout.rdbuf(&ind_sb);

		std::cout << PeriodicTable::getSingleton() << "\n\n"
			<< table["Cl"] << "\n";

		std::cout.rdbuf(old_sb);
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
