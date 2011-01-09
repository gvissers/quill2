#include <fstream>
#include "PeriodicTable.hh"

int main()
{
	try
	{
		PeriodicTable table("data/elements.dat");
		std::cout << PeriodicTable::getSingleton() << "\n";
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
