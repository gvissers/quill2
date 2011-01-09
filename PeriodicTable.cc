#include <fstream>
#include <sstream>
#include "PeriodicTable.hh"
#include "constants.hh"
#include "exceptions.hh"
#include "IndentingStreambuf.hh"

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
	SymbolElementMap::const_iterator it = _by_symbol.find(elem);
	if (it == _by_symbol.end())
		throw UnknownElement(elem);
	return *it->second;
}

std::ostream& PeriodicTable::print(std::ostream& os) const
{
	os << "PeriodicTable (\n" << indent;
	for (ElementList::const_iterator it = _elements.begin();
		it != _elements.end(); ++it)
	{
		os << *it << "\n";
	}
	return os << dedent << ")";
}

void PeriodicTable::insert(const Element& elem)
{
	SymbolElementMap::iterator it = _by_symbol.find(elem.symbol());
	Element *elem_ptr;
	if (it != _by_symbol.end())
	{
		elem_ptr = it->second;
		if (elem_ptr->number() != elem.number())
			_by_number.erase(elem_ptr->number());
		if (elem_ptr->name() != elem.name())
			_by_name.erase(elem_ptr->name());
		*elem_ptr = elem;
	}
	else
	{
		_elements.push_back(elem);
		elem_ptr = &_elements.back();
		_by_symbol.insert(std::make_pair(elem.symbol(), elem_ptr));
	}
	_by_number.insert(std::make_pair(elem.number(), elem_ptr));
	_by_name.insert(std::make_pair(elem.name(), elem_ptr));
}
