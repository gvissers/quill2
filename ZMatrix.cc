#include <Eigen/Dense>
#include "ZMatrix.hh"
#include "IndentingOStream.hh"
#include "manipulators.hh"
#include "support.hh"

namespace {

struct ElementGetter
{
	ZMatrix *zmat;
	std::string elem;

	ElementGetter(ZMatrix* zmat): zmat(zmat), elem() {}
};

struct IndexGetter
{
	const ZMatrix *zmat;
	int idx;
	std::string lbl;

	IndexGetter(ZMatrix* zmat): zmat(zmat), idx(-1), lbl() {}
};

std::istream& operator>>(std::istream& is, ElementGetter& eg)
{
	is >> element(eg.elem, false);
	if (!is.fail())
	{
		eg.zmat->addLabel(eg.elem, eg.zmat->size());
		rtrim(eg.elem, std::tr1::bind(std::isdigit<char>, _1, std::locale()));
	}
	return is;
}

std::istream& operator>>(std::istream& is, IndexGetter& ig)
{
	ig.lbl.clear();
	is >> ig.idx;
	if (!is.fail())
	{
		ig.idx--;
	}
	else
	{
		is.clear();
		is >> element(ig.lbl, false);
		if (!is.fail())
			ig.idx = ig.zmat->findLabel(ig.lbl);
		else
			ig.idx = -1;
	}
	return is;
}

} // namespace

std::ostream& ZMatrix::print(std::ostream& os) const
{
	os << "ZMatrix (\n" << indent;

	AtomList::const_iterator it = _positions.begin();
	if (it != _positions.end())
	{
		os << it->symbol << "\n";
		if (++it != _positions.end())
		{
			os << it->symbol << "\t"
				<< (it->idx1+1) << "\t" << it->r
				<< "\n";
			if (++it != _positions.end())
			{
				os << it->symbol
					<< "\t" << (it->idx1+1) << "\t" << it->r
					<< "\t" << (it->idx2+1) << "\t" << radToDeg(it->theta)
					<< "\n";
				for (++it; it != _positions.end(); ++it)
				{
					os << it->symbol
						<< "\t" << (it->idx1+1) << "\t" << it->r
						<< "\t" << (it->idx2+1) << "\t" << radToDeg(it->theta)
						<< "\t" << (it->idx3+1) << "\t" << radToDeg(it->phi)
						<< "\n";
				}
			}
		}
	}
	os << dedent << ")";

	return os;
}

JobIStream& ZMatrix::scan(JobIStream& is)
{
	is >> getline;
	if (is.eof())
		return is;

	ElementGetter eg(this);
	is >> eg;

	if (is.fail())
	{
		is.ungetLastLine();
		is.setstate(std::ios::failbit);
		return is;
	}
	_positions.push_back(AtomPos(eg.elem, -1, 0, -1, 0, -1, 0));

	is >> getline;
	if (is.eof())
		return is;

	IndexGetter ig1(this);
	double r;
	is >> eg >> ig1 >> r;
	if (is.fail())
	{
		is.ungetLastLine();
		return is;
	}
	checkIndex(ig1.idx, ig1.lbl);

	_positions.push_back(AtomPos(eg.elem, ig1.idx, r, 0, 0, 0, 0));

	is >> getline;
	if (is.eof())
		return is;

	IndexGetter ig2(this);
	double theta;
	is >> eg >> ig1 >> r >> ig2 >> theta;
	if (is.fail())
	{
		is.ungetLastLine();
		return is;
	}
	checkIndex(ig1.idx, ig1.lbl);
	checkIndex(ig2.idx, ig2.lbl);

	_positions.push_back(AtomPos(eg.elem, ig1.idx, r,
		ig2.idx, degToRad(theta), 0, 0));

	IndexGetter ig3(this);
	double phi;
	while (true)
	{
		is >> getline;
		if (is.eof())
			return is;
		is >> eg >> ig1 >> r >> ig2 >> theta >> ig3 >> phi;
		if (is.fail())
		{
			is.ungetLastLine();
			return is;
		}
		checkIndex(ig1.idx, ig1.lbl);
		checkIndex(ig2.idx, ig2.lbl);
		checkIndex(ig3.idx, ig3.lbl);

		_positions.push_back(AtomPos(eg.elem, ig1.idx, r,
			ig2.idx, degToRad(theta), ig3.idx, degToRad(phi)));
	}
}

void ZMatrix::checkIndex(int idx, const std::string& lbl) const
{
	if (idx < 0 || idx >= size())
	{
		if (!lbl.empty())
			throw InvalidLabel(lbl);
		else
			throw InvalidLabel(idx+1);
	}
}

void ZMatrix::fillGeometry(Geometry *geom)
{
	geom->resize(size());

	for (int i = 0; i < size(); i++)
	{
		const AtomPos& ap = _positions[i];
		double ct = cos(ap.theta), st = sin(ap.theta);
		double cp = cos(ap.phi), sp = sin(ap.phi);
		geom->setAtom(i, ap.symbol, ap.r*st*cp, ap.r*st*sp, ap.r*ct);
	}

	if (size() > 2)
	{
		const AtomPos& ap = _positions[2];
		if (ap.idx1 == 1)
			geom->setPosition(2, geom->position(1)-geom->position(2));

		Eigen::Vector3d posA, posB, posC;
		Eigen::Vector3d cb, n;
		Eigen::Matrix3d M;
		for (int i = 3; i < size(); i++)
		{
			const AtomPos& ap = _positions[i];

			posC = geom->position(ap.idx1);
			posB = geom->position(ap.idx2);
			posA = geom->position(ap.idx3);

			cb = (posB - posC).normalized();
			n = (posA-posB).cross(cb).normalized();

			M.col(0) = cb.cross(n);
			M.col(1) = n;
			M.col(2) = cb;

			geom->setPosition(i, M*geom->position(i) + posC);
		}
	}
}
