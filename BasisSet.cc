#include <sstream>
#include <cstdio>
#include <util.hh>
#include <tr1/functional>
#include "BasisSet.hh"
#include "Deleter.hh"
#include "CGTODef.hh"
#include "exceptions.hh"
#include "support.hh"
#include "manipulators.hh"

using namespace std::tr1::placeholders;

template<>
void BasisSet::readElement<BasisSet::Turbomole>(FilteringLineIStream<CommentFilter>& fis)
{
	std::string elem, name, sep;

	fis >> expectline >> element(elem, false) >> name;
	if (fis.fail())
		throw ParseError(fis.lineNumber(), Turbomole,
			"Element name and basis name expected");

	fis >> expectline >> sep;
	if (sep != "*")
		throw ParseError(fis.lineNumber(), Turbomole,
			"Expected \"*\"");
	BFList& elem_funs = _elements[elem];
	while (true)
	{
		fis >> expectline >> sep;
		if (sep == "*")
			break;
		fis.ungetLastLine();

		AbstractBFDef* f = readTurbomoleBF(fis);
		elem_funs.push_back(f);
	}
}

template<>
void BasisSet::readElement<BasisSet::Molpro>(FilteringLineIStream<CommentFilter>& fis)
{
	std::string shell, elem;
	std::vector<double> widths;

	fis >> expectline >> shell >> element(elem, false);
	if (fis.fail())
		throw ParseError(fis.lineNumber(), Molpro,
			"Expected shell, element name, and primitive widths");
	while (!fis.eof())
	{
		double width;
		fis >> width;
		if (!fis.fail())
			widths.push_back(width);
	}

	std::vector<double> weights;
	std::vector< std::pair<double, double> > ww;
	char ch;
	BFList& elem_funs = _elements[elem];
	while (true)
	{
		fis >> expectline >> ch;
		if (ch != 'c')
		{
			fis.ungetLastLine();
			break;
		}

		int istart, iend;
		fis >> istart >> ch >> iend;
		if (fis.fail() || ch !='.')
			throw ParseError(fis.lineNumber(), Molpro,
				"Expected start end end indices of primitive widths");
		if (istart < 1 || istart > int(widths.size())
			|| iend < istart || iend > int(widths.size()))
			throw ParseError(fis.lineNumber(), Molpro,
				"Invalid primitive index");

		weights.clear();
		for (int i = istart; i <= iend; i++)
		{
			double weight;
			fis >> weight;
			if (fis.fail())
				throw ParseError(fis.lineNumber(), Molpro,
					"Expected primitive weight");
			weights.push_back(weight);
		}

		ww.clear();
		std::transform(weights.begin(), weights.end(),
			widths.begin()+(istart-1), std::back_inserter(ww),
			std::make_pair<double, double>);

		AbstractBFDef *bf = contractedGaussian(shell[0], ww);
		elem_funs.push_back(bf);
	}
}

template<>
void BasisSet::readElement<BasisSet::Dalton>(FilteringLineIStream<CommentFilter>& fis)
{
	static std::string last_elem;
	static char last_shell;

	std::string elem;
	char shell;
	int nprim, nbf;

	fis >> expectline >> element(elem, false) >> nprim >> nbf;
	if (fis.fail())
		throw ParseError(fis.lineNumber(), Dalton,
			"Expected element name, number of primitives, and number of contractions");

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
				throw ParseError(fis.lineNumber(), Dalton,
					"Orbital shell is higher than supported");
		}
	}

	std::vector< std::vector< std::pair<double, double> > > wws(nbf);
	for (int iprim = 0; iprim < nprim; iprim++)
	{
		double width;

		fis >> expectline >> width;
		if (!fis.good())
			throw ParseError(fis.lineNumber(), Dalton,
				"Expected width of primitive");
		for (int ibf = 0; ibf < nbf; ibf++)
		{
			double weight;
			fis >> weight;
			if (!fis.good())
				throw ParseError(fis.lineNumber(), Dalton,
					"Expected contraction coefficient");
			if (weight != 0.0)
				wws[ibf].push_back(std::make_pair(weight, width));
		}
	}

	BFList& elem_funs = _elements[elem];
	for (int ibf = 0; ibf < nbf; ibf++)
	{
		AbstractBFDef *bf = contractedGaussian(shell, wws[ibf]);
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
	std::string str;
	CommentFilter filter('#');
	FilteringLineIStream<CommentFilter> fis(is, &filter);

	try
	{
		fis >> expectline >> str;
		if (str != "$basis")
			throw ParseError(fis.lineNumber(), Turbomole,
				"Expected \"$basis\"");
		fis >> expectline >> str;
		if (str != "*")
			throw ParseError(fis.lineNumber(), Turbomole,
				"Expected \"*\"");
		while (true)
		{
			fis >> expectline >> str;
			if (str == "$end")
				break;

			fis.ungetLastLine();
			readElement<Turbomole>(fis);
		}
	}
	catch (const UnexpectedEOF&)
	{
		throw ParseError(fis.lineNumber(), Turbomole,
			"Unexpected end of file");
	}
}

template<>
void BasisSet::read<BasisSet::Molpro>(std::istream& is)
{
	std::string str;
	CommentFilter filter('!');
	FilteringLineIStream<CommentFilter> fis(is, &filter);

	try
	{
		fis >> expectline;
		str = fis.line();
		remove(str, std::tr1::bind(std::isspace<char>, _1, std::locale()));
		if (str != "basis={")
			throw ParseError(fis.lineNumber(), Molpro,
				"Expected \"basis={\"");

		while (true)
		{
			fis >> expectline >> str;
			if (str == "}")
				break;

			fis.ungetLastLine();
			readElement<Molpro>(fis);
		}
	}
	catch (const UnexpectedEOF&)
	{
		throw ParseError(fis.lineNumber(), Molpro,
			"Unexpected end of file");
	}
}

template<>
void BasisSet::read<BasisSet::Dalton>(std::istream& is)
{
	CommentFilter filter('!');
	FilteringLineIStream<CommentFilter> fis(is, &filter);
	while (true)
	{
		fis >> getline;
		if (fis.eof()) break;

		fis.ungetLastLine();
		readElement<Dalton>(fis);
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

AbstractBFDef* BasisSet::readTurbomoleBF(FilteringLineIStream<CommentFilter>& fis)
{
	int nr_prim;
	char shell;

	fis >> expectline >> nr_prim >> shell;
	if (fis.fail())
		throw ParseError(fis.lineNumber(), Turbomole,
			"Number of primitives and shell expected");

	std::vector< std::pair<double, double> > ww;
	for (int i = 0; i < nr_prim; i++)
	{
		double width, weight;

		fis >> expectline >> width >> weight;
		if (fis.fail())
			throw ParseError(fis.lineNumber(), Turbomole,
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
			Deleter<AbstractBFDef>());
	}
}

AbstractBFDef* BasisSet::contractedGaussian(char shell,
	const std::vector< std::pair<double, double> >& ww)
{
	switch (shell)
	{
		case 's': return new CGTODef<0>(ww);
		case 'p': return new CGTODef<1>(ww);
		case 'd': return new CGTODef<2>(ww);
		case 'f': return new CGTODef<3>(ww);
		case 'g': return new CGTODef<4>(ww);
		case 'h': return new CGTODef<5>(ww);
		case 'i': return new CGTODef<6>(ww);
		default: throw UnknownShell(shell);
	}
}

