#include <fstream>
#include <sstream>
#include "PeriodicTable.hh"
#include "constants.hh"
#include "exceptions.hh"

SINGLETON_OBJECT(PeriodicTable);

PeriodicTable::PeriodicTable(const std::string& fname):
	Li::Singleton<PeriodicTable>(), _elements()
{
	std::ifstream is(fname.c_str());
	if (!is.good())
		throw NoFile(fname);
	read(is);
	is.close();
}

void PeriodicTable::read(std::istream& is)
{
	char line[1024];

	while (true)
	{
		is.getline(line, sizeof(line));
		if (is.eof()) break;
		if (line[0] == '#') continue;

		std::istringstream ss(line);
		std::string symbol, name;
		char comma;
		unsigned int nr;
		double mass, covrad, vdwrad = 0.0;
		bool H_bonding;

		ss >> symbol >> name >> nr >> comma >> mass >> comma
			>> H_bonding >> comma >> covrad;
		if (!ss.good()) // bad line
			continue;

		if (symbol[symbol.size()-1] == ',')
			symbol.erase(symbol.size()-1);
		if (name[name.size()-1] == ',')
			name.erase(name.size()-1);

		comma = 0;
		ss >> comma;
		if (ss.good() && comma == ',')
			ss >> vdwrad;

		insert(Element(symbol, name, nr, mass*Constants::amu,
			H_bonding, covrad*Constants::Ang,
			vdwrad*Constants::Ang));
	}
}

const Element& PeriodicTable::operator[](const std::string& elem) const
{
	ElementMap::const_iterator it = _elements.find(elem);
	if (it == _elements.end())
		throw UnknownElement(elem);
	return it->second;
}

Element& PeriodicTable::operator[](const std::string& elem)
{
	ElementMap::iterator it = _elements.find(elem);
	if (it == _elements.end())
		throw UnknownElement(elem);
	return it->second;
}

std::ostream& PeriodicTable::print(std::ostream& os) const
{
	os << "PeriodicTable (\n";
	for (ElementMap::const_iterator it = _elements.begin();
		it != _elements.end(); ++it)
	{
		os << it->second << "\n";
	}
	return os << ")";
}
