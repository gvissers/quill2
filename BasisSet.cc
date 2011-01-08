#include <sstream>
#include "BasisSet.hh"
#include "Deleter.hh"
#include "CGTO.hh"

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


std::ostream& BasisSet::print(std::ostream& os) const
{
	os << "BasisSet (\n";
	for (BFMap::const_iterator elem_it = _elements.begin();
		elem_it != _elements.end(); ++elem_it)
	{
		os << elem_it->first << " (\n";
		for (BFList::const_iterator fun_it = elem_it->second.begin();
			fun_it != elem_it->second.end(); ++fun_it)
		{
			os << **fun_it << "\n";
		}
		os << ")\n";
	}
	return os << ")";
}

void BasisSet::readAuto(std::istream& is)
{
	readTurbomole(is);
}

void BasisSet::readTurbomole(std::istream& is)
{
	const char* line;
	LineGetter getter(is, '#');

	try
	{
		line = getter.next();
		if (std::strcmp(line, "$basis") != 0)
			throw ParseError(getter.lineNumber(), Turbomole,
				"Expected \"$basis\"");

		line = getter.next();
		if (std::strcmp(line, "*") != 0)
			throw ParseError(getter.lineNumber(), Turbomole,
				"Expected \"*\"");
		while (true)
		{
			line = getter.next();
			if (std::strcmp(line, "$end") == 0)
				break;

			readTurbomoleElement(getter);
		}

		if (std::strcmp(getter.line(), "$end") != 0)
			throw ParseError(getter.lineNumber(), Turbomole,
				"Expected \"$end\"");
	}
	catch (const LineGetter::UnexpectedEOF&)
	{
		throw ParseError(getter.lineNumber(), Turbomole,
			"Unexpected end of file");
	}
}

void BasisSet::readTurbomoleElement(LineGetter& getter)
{
	std::string elem, name;
	std::istringstream is;

	const char* line = getter.line();
	is.str(line);
	is >> elem >> name;
	if (is.fail())
		throw ParseError(getter.lineNumber(), Turbomole,
			"Element name and basis name expected");
	elem = ucFirst(elem);

	line = getter.next();
	if (std::strcmp(line, "*") != 0)
		throw ParseError(getter.lineNumber(), Turbomole,
			"Expected \"*\"");
	BFList& elem_funs = _elements[elem];
	while (true)
	{
		line = getter.next();
		if (std::strcmp(line, "*") == 0)
			break;

		AbstractBF* f = readTurbomoleBF(getter);
		elem_funs.push_back(f);
	}
}

AbstractBF* BasisSet::readTurbomoleBF(LineGetter& getter)
{
	int nr_prim;
	char shell;

	const char *line = getter.line();
	std::istringstream is(line);
	is >> nr_prim >> shell;
	if (is.fail())
		throw ParseError(getter.lineNumber(), Turbomole,
			"Number of primitives and shell expected");

	std::vector< std::pair<double, double> > ww;
	for (int i = 0; i < nr_prim; i++)
	{
		double width, weight;

		line = getter.next();
		is.clear();
		is.str(line);
		is >> width >> weight;
		if (is.fail())
			throw ParseError(getter.lineNumber(), Turbomole,
				"Width and weight of primitive expected");

		ww.push_back(std::make_pair(weight, width));
	}

	AbstractBF *bf;
	switch (shell)
	{
		case 's': bf = new CGTO<0>(ww); break;
		case 'p': bf = new CGTO<1>(ww); break;
		case 'd': bf = new CGTO<2>(ww); break;
		case 'f': bf = new CGTO<3>(ww); break;
		case 'g': bf = new CGTO<4>(ww); break;
		case 'h': bf = new CGTO<5>(ww); break;
		case 'i': bf = new CGTO<6>(ww); break;
		default:
			throw ParseError(getter.lineNumber(), Turbomole,
				std::string("Unknown shell \'") + shell + "'");
	}
	return bf;
}

void BasisSet::clear()
{
	for (BFMap::iterator it = _elements.begin();
		it != _elements.end(); ++it)
	{
		std::for_each(it->second.begin(), it->second.end(),
			Deleter<AbstractBF>());
	}
}


