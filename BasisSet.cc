#include <sstream>
#include <cstdio>
#include <util.hh>
#include "BasisSet.hh"
#include "Deleter.hh"
#include "CGTO.hh"

template<>
void BasisSet::readElement<BasisSet::Turbomole>(LineGetter& getter)
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

template<>
void BasisSet::readElement<BasisSet::Molpro>(LineGetter& getter)
{
	const char* line = getter.line();
	Li::Vector<std::string> fields = Li::explode(",", line);
	if (fields.size() < 3)
		throw ParseError(getter.lineNumber(), Molpro,
			"Expected shell, element name, and primitive widths");

	std::string shell = Li::strip(fields[0]);
	std::string elem = ucFirst(Li::strip(fields[1]));
	std::vector<double> widths;
	std::transform(fields.begin()+2, fields.end(),
		std::back_inserter(widths), Li::fromString<double>);

	std::vector<double> weights;
	std::vector< std::pair<double, double> > ww;
	BFList& elem_funs = _elements[elem];
	while (true)
	{
		line = getter.next();
		fields = Li::explode(",", line);
		if (fields.size() < 3 || Li::strip(fields[0]) != "c")
			break;

		int istart, iend;
		if (std::sscanf(fields[1].c_str(), "%d.%d", &istart, &iend) != 2)
			throw ParseError(getter.lineNumber(), Molpro,
				"Expected start end end indices of primitive widths");
		if (istart < 1 || istart > int(widths.size())
			|| iend < istart || iend > int(widths.size()))
			throw ParseError(getter.lineNumber(), Molpro,
				"Invalid primitive index");

		weights.clear();
		std::transform(fields.begin()+2, fields.end(),
			std::back_inserter(weights), Li::fromString<double>);

		ww.clear();
		std::transform(weights.begin(), weights.end(),
			widths.begin()+(istart-1), std::back_inserter(ww),
			std::make_pair<double, double>);

		AbstractBF *bf = contractedGaussian(shell[0], ww);
		elem_funs.push_back(bf);
	}
}

template<>
void BasisSet::readElement<BasisSet::Dalton>(LineGetter& getter)
{
	static std::string last_elem;
	static char last_shell;

	std::string elem;
	char shell;
	int nprim, nbf;

	const char* line = getter.line();
	std::istringstream is(line);
	is >> elem >> nprim >> nbf;
	if (is.fail())
		throw ParseError(getter.lineNumber(), Dalton,
			"Expected element name, number of primitives, and number of contractions");
	elem = ucFirst(elem);

	if (elem != last_elem)
	{
		shell = last_shell = 's';
		last_elem = elem;
	}
	else
	{
		switch (last_shell)
		{
			case 's': shell = last_shell = 'p'; break;
			case 'p': shell = last_shell = 'd'; break;
			case 'd': shell = last_shell = 'f'; break;
			case 'f': shell = last_shell = 'g'; break;
			case 'g': shell = last_shell = 'h'; break;
			case 'h': shell = last_shell = 'i'; break;
			default:
				throw ParseError(getter.lineNumber(), Dalton,
					"Orbital shell is higher than supported");
		}
	}

	std::vector< std::vector< std::pair<double, double> > > wws(nbf);
	for (int iprim = 0; iprim < nprim; iprim++)
	{
		double width;

		line = getter.next();
		is.clear();
		is.str(line);

		is >> width;
		if (!is.good())
			throw ParseError(getter.lineNumber(), Dalton,
				"Expected width of primitive");
		for (int ibf = 0; ibf < nbf; ibf++)
		{
			double weight;
			is >> weight;
			if (!is.good())
				throw ParseError(getter.lineNumber(), Dalton,
					"Expected contraction coefficient");
			if (weight != 0.0)
				wws[ibf].push_back(std::make_pair(weight, width));
		}
	}

	BFList& elem_funs = _elements[elem];
	for (int ibf = 0; ibf < nbf; ibf++)
	{
		AbstractBF *bf = contractedGaussian(shell, wws[ibf]);
		elem_funs.push_back(bf);
	}
}

template<>
void BasisSet::read<BasisSet::Auto>(std::istream& is)
{
}

template<>
void BasisSet::read<BasisSet::Turbomole>(std::istream& is)
{
	const char* line;
	LineGetter getter(is, "#");

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

			readElement<Turbomole>(getter);
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

template<>
void BasisSet::read<BasisSet::Molpro>(std::istream& is)
{
	const char* line;
	LineGetter getter(is, "!");

	try
	{
		line = getter.next();
		if (strcmp(line, "basis={") != 0)
			throw ParseError(getter.lineNumber(), Molpro,
				"Expected \"basis={\"");

		line = getter.next();
		while (true)
		{
			if (strcmp(getter.line(), "}") == 0)
				break;

			readElement<Molpro>(getter);
		}

		if (std::strcmp(getter.line(), "}") != 0)
			throw ParseError(getter.lineNumber(), Turbomole,
				"Expected \"}\"");
	}
	catch (const LineGetter::UnexpectedEOF&)
	{
		throw ParseError(getter.lineNumber(), Molpro,
			"Unexpected end of file");
	}
}

template<>
void BasisSet::read<BasisSet::Dalton>(std::istream& is)
{
	LineGetter getter(is, "!");
	while (true)
	{
		const char* line = getter.next(false);
		if (!line) break;

		readElement<Dalton>(getter);
	}
}

void BasisSet::read(std::istream& is)
{
	read<Auto>(is);
}

std::ostream& BasisSet::print(std::ostream& os) const
{
	os << "BasisSet (\n" << indent;
	for (BFMap::const_iterator elem_it = _elements.begin();
		elem_it != _elements.end(); ++elem_it)
	{
		os << elem_it->first << " (\n" << indent;
		for (BFList::const_iterator fun_it = elem_it->second.begin();
			fun_it != elem_it->second.end(); ++fun_it)
		{
			os << **fun_it << "\n";
		}
		os << dedent << ")\n";
	}
	return os << dedent << ")";
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

	return contractedGaussian(shell, ww);
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

AbstractBF* BasisSet::contractedGaussian(char shell,
	const std::vector< std::pair<double, double> >& ww)
{
	switch (shell)
	{
		case 's': return new CGTO<0>(ww);
		case 'p': return new CGTO<1>(ww);
		case 'd': return new CGTO<2>(ww);
		case 'f': return new CGTO<3>(ww);
		case 'g': return new CGTO<4>(ww);
		case 'h': return new CGTO<5>(ww);
		case 'i': return new CGTO<6>(ww);
		default: throw UnknownShell(shell);
	}
}

