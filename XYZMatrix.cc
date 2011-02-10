#include "XYZMatrix.hh"
#include "io/manipulators.hh"

std::ostream& XYZMatrix::print(std::ostream& os) const
{
	os << "XYZMatrix (\n" << indent;
	for (AtomList::const_iterator nit = _positions.begin();
		nit != _positions.end(); ++nit)
	{
		os << nit->symbol << "\t" << nit->x
			<< "\t" << nit->y << "\t" << nit->z << "\n";
	}
	os << dedent << ")";
	return os;
}

JobIStream& XYZMatrix::scan(JobIStream& is)
{
	std::string elem;
	double x, y, z;

	while (true)
	{
		is >> getline;
		if (is.eof())
			break;

		is >> element(elem, false) >> x >> y >> z;
		if (is.fail())
		{
			is.ungetLastLine();
			break;
		}

		_positions.push_back(AtomPos(elem, x, y, z));
	}

	if (size() == 0)
		is.setstate(std::ios::failbit);
	return is;
}

void XYZMatrix::fillGeometry(Geometry *geom)
{
	geom->resize(size());
	for (int i = 0; i < size(); i++)
	{
		const AtomPos& ap = _positions[i];
		geom->setAtom(i, ap.symbol, ap.x, ap.y, ap.z);
	}
}
