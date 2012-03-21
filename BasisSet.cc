#include <sstream>
#include <fstream>
#include <util.hh>
#include <filesystem.hh>
#include "BasisSet.hh"
#include "Deleter.hh"
#include "CGTODef.hh"
#include "PeriodicTable.hh"
#include "exceptions.hh"
#include "support.hh"
#include "io/manipulators.hh"

void BasisSet::expand(const Geometry& geom, Basis *basis) const
{
	for (int i = 0; i < geom.size(); i++)
	{
		BFMap::const_iterator eit = _elements.find(geom.symbol(i));
		if (eit == _elements.end())
			throw NoFunctions(geom.symbol(i));

		Eigen::Vector3d pos = geom.position(i);
		for (BFList::const_iterator fit = eit->second.begin();
			fit != eit->second.end(); ++fit)
		{
			(*fit)->expand(i, pos, basis);
		}
	}
}

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
			[](double w, double a) { return std::make_pair(w,a); } );

		AbstractBFDef *bf = contractedGaussian(shell[0], ww);
		elem_funs.push_back(bf);
	}
}

template<>
void BasisSet::readElement<BasisSet::Dalton>(FilteringLineIStream<CommentFilter>& fis)
{
	static std::string last_elem;
	static int last_l;

	std::string elem;
	int l, nprim, nbf;

	fis >> expectline >> element(elem, false) >> nprim >> nbf;
	if (fis.fail())
		throw ParseError(fis.lineNumber(), Dalton,
			"Expected element name, number of primitives, and number of contractions");

	if (elem != last_elem)
	{
		last_elem = elem;
		l = last_l = 0;
	}
	else
	{
		l = ++last_l;
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
			if (weight != 0)
				wws[ibf].push_back(std::make_pair(weight, width));
		}
	}

	BFList& elem_funs = _elements[elem];
	for (int ibf = 0; ibf < nbf; ibf++)
	{
		AbstractBFDef *bf = contractedGaussian(l, wws[ibf]);
		elem_funs.push_back(bf);
	}
}

template<>
void BasisSet::readElement<BasisSet::Molcas>(FilteringLineIStream<CommentFilter>& fis)
{
	std::string kw1, kw2, lbl;

	fis >> expectline >> kw1 >> kw2;
	if (kw1 != "Basis" || kw2 != "set")
		throw ParseError(fis.lineNumber(), Molcas,
			"Expected \"Basis set\"");

	char c;
	fis >> expectline >> lbl >> c >> kw1;
	if (fis.fail() || c != '/' || kw1 != "inline")
		throw ParseError(fis.lineNumber(), Molcas,
			"Expected atom label and inline keyword");

	double charge;
	int lmax;
	fis >> expectline >> charge >> lmax;
	if (fis.fail())
		throw ParseError(fis.lineNumber(), Molcas,
			"Expected nuclear charge and maximum angular momentum");

	int elem_nr = charge;
	const std::string& elem = PeriodicTable::singleton().findByNumber(elem_nr).symbol();
	BFList& elem_funs = _elements[elem];

	std::vector<double> widths;
	std::vector< std::vector< std::pair<double, double> > > wws;
	for (int l = 0; l <= lmax; l++)
	{
		int nprim, nbf;
		fis >> expectline >> nprim >> nbf;
		if (fis.fail())
			throw ParseError(fis.lineNumber(), Molcas,
				"Expected number of primitives and number of functions");

		widths.clear();
		for (int iprim = 0; iprim < nprim; iprim++)
		{
			double width;
			fis >> expectline >> width;
			if (fis.fail())
				throw ParseError(fis.lineNumber(), Molcas,
					"Expected width of primitive");

			widths.push_back(width);
		}

		wws.clear();
		wws.resize(nbf);
		for (int iprim = 0; iprim < nprim; iprim++)
		{
			double width = widths[iprim];
			fis >> expectline;
			for (int ibf = 0; ibf < nbf; ibf++)
			{
				double weight;
				fis >> weight;
				if (fis.fail())
					throw ParseError(fis.lineNumber(), Molcas,
						"Expected contraction coefficient");
				if (weight != 0)
					wws[ibf].push_back(std::make_pair(weight, width));
			}
		}

		for (int ibf = 0; ibf < nbf; ibf++)
		{
			AbstractBFDef *bf = contractedGaussian(l, wws[ibf]);
			elem_funs.push_back(bf);
		}
	}
}

void BasisSet::scan(std::istream& is, Format hint)
{
	switch (hint)
	{
		case Dalton:
			if (tryScan<Dalton>(is))
				return;
			break;
		case Molcas:
			if (tryScan<Molcas>(is))
				return;
			break;
		case Molpro:
			if (tryScan<Molpro>(is))
				return;
			break;
		case Turbomole:
			if (tryScan<Turbomole>(is))
				return;
			break;
		default:
			/* nothing */ ;
	}

	if (tryScan<Dalton>(is))
		return;
	if (tryScan<Molcas>(is))
		return;
	if (tryScan<Molpro>(is))
		return;
	if (tryScan<Turbomole>(is))
		return;

	throw ParseError(0, Auto, "Failed to parse basis set in any format");
}

template<>
void BasisSet::scan<BasisSet::Turbomole>(std::istream& is)
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
void BasisSet::scan<BasisSet::Molpro>(std::istream& is)
{
	std::string str;
	CommentFilter filter('!');
	FilteringLineIStream<CommentFilter> fis(is, &filter);

	try
	{
		fis >> expectline;
		str = fis.line();
		remove(str, std::bind(std::isspace<char>, _1, std::locale()));
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
void BasisSet::scan<BasisSet::Dalton>(std::istream& is)
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

template <>
void BasisSet::scan<BasisSet::Molcas>(std::istream& is)
{
	CommentFilter filter("!*");
	FilteringLineIStream<CommentFilter> fis(is, &filter);
	while (true)
	{
		fis >> getline;
		if (fis.eof()) break;

		fis.ungetLastLine();
		readElement<Molcas>(fis);
	}
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
	static const std::string shells("spdfghi");

	std::string::size_type pos = shells.find(shell);
	if (pos == std::string::npos)
		throw UnknownShell(shell);

	return contractedGaussian(int(pos), ww);
}

AbstractBFDef* BasisSet::contractedGaussian(int l,
	const std::vector< std::pair<double, double> >& ww)
{
	switch (l)
	{
		case 0: return new CGTODef<0>(ww);
		case 1: return new CGTODef<1>(ww);
		case 2: return new CGTODef<2>(ww);
		case 3: return new CGTODef<3>(ww);
		case 4: return new CGTODef<4>(ww);
		case 5: return new CGTODef<5>(ww);
		case 6: return new CGTODef<6>(ww);
		default: throw ShellTooHigh(l);
	}
}

bool BasisSet::findAndScan(const std::string& name, const std::string& dir)
{
	std::string lcname = lower(name);
	for (Li::DirectoryIterator it(dir); !it.end(); ++it)
	{
		std::string fname = *it;
		std::string fullname = Li::joinPath(dir, fname);
		if (!Li::isRegularFile(fullname))
			continue;

		std::string ext = lower(Li::getExtension(fname));
		std::string base = lower(fname.substr(0, fname.size()-ext.size()));
		if (base == lcname)
		{
			Format hint = Auto;
			if (ext == ".dalton")
				hint = Dalton;
			else if (ext == ".molcas")
				hint = Molcas;
			else if (ext == ".molpro")
				hint = Molpro;
			else if (ext == ".turbomole")
				hint = Turbomole;

			std::ifstream is(fullname);
			try
			{
				scan(is, hint);
				is.close();
				std::cout << "Read basis set from file " << fullname  << "\n";
				return true;
			}
			catch (const ParseError&)
			{
				/* ignore and continue with the next file name */
			}
			is.close();
		}
	}

	return false;
}