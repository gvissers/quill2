#include "Geometry.hh"
#include "PeriodicTable.hh"
#include "exceptions.hh"
#include "support.hh"

JobIStream& Geometry::scan(JobIStream& is)
{
	std::string elem;
	double x, y, z;

	is >> getline;
	if (is.eof())
		throw UnexpectedEOF();

	is >> elem >> x >> y >> z;
	bool failed = is.fail();
	is.ungetLastLine();
	return failed ? scanZMatrix(is) : scanXYZ(is);
}

JobIStream& Geometry::scanXYZ(JobIStream& is)
{
	bool ok;

	scanXYZLine(is, true);
	do
	{
		ok = scanXYZLine(is);
	} while (ok);

	return is;
}

bool Geometry::scanXYZLine(JobIStream& is, bool except)
{
	std::string elem;
	double x, y, z;

	is >> getline;
	if (is.eof())
	{
		if (except) throw UnexpectedEOF();
		return false;
	}

	is >> elem >> x >> y >> z;
	if (is.fail())
	{
		is.ungetLastLine();
		if (except) throw ParseError();
		return false;
	}

	elem = ucFirst(elem);
	if (!PeriodicTable::getSingleton().exists(elem))
	{
		is.ungetLastLine();
		if (except) throw UnknownElement(elem);
		return false;
	}

	_nuclei.push_back(Nucleus(elem, x, y, z));
	return true;
}

JobIStream& Geometry::scanZMatrix(JobIStream& is)
{
	return is;
}
