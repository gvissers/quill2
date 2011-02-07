#include "Geometry.hh"
#include "PeriodicTable.hh"
#include "IndentingOStream.hh"
#include "XYZMatrix.hh"
#include "ZMatrix.hh"
#include "manipulators.hh"
#include "exceptions.hh"

void Geometry::setAtom(int idx, const std::string& symbol,
	double x, double y, double z)
{
	checkIndex(idx);

	const Element& elem = PeriodicTable::getSingleton().findBySymbol(symbol);

	_positions.col(idx) << x, y, z;
	_masses(idx) = elem.mass();
	_charges(idx) = elem.number();
	_symbols[idx] = symbol;
}

std::ostream& Geometry::print(std::ostream& os) const
{
	os << "Geometry (\n" << indent;
	for (int i = 0; i < size(); i++)
		os << _charges(i) << "\t" << _masses(i) << "\t"
			<< _symbols[i] << "\t"
			<< _positions.col(i).transpose() << "\n";
	os << dedent << ")";
	return os;
}

JobIStream& Geometry::scan(JobIStream& is)
{
	is >> getline;
	if (is.eof())
		throw UnexpectedEOF();

	std::string elem;
	double x, y, z;
	is >> element(elem, false) >> x >> y >> z;
	if (!is.fail())
	{
		is.ungetLastLine();

		XYZMatrix mat;
		is >> mat;
		mat.fillGeometry(this);
	}
	else
	{
		is.ungetLastLine();

		ZMatrix mat;
		is >> mat;
		mat.fillGeometry(this);
	}

	return is;
}
