#include <algorithm>
#include <functional>
#include "BasisSet.hh"

void BasisSet::read(std::istream& is, Format format)
{
	switch (format)
	{
		case Auto:
			readAuto(is);
			break;
		case Turbomole:
			readTurbomole(is);
			break;
		default:
			throw UnknownFormat(format);
	}
}

void BasisSet::readAuto(std::istream& is)
{
}

void BasisSet::readTurbomole(std::istream& is)
{
	char buf[256];
	int line = 1;

	while (true)
	{
		is.getline(buf, sizeof(buf));
		if (!is.good()) break;
		if (buf[0] && buf[0] != '#')
	}
}

void BasisSet::clear()
{
	for (BFMap::iterator it = _functions.begin();
		it != _functions.end(); ++it)
	{
		delete it->second;
	}
}


