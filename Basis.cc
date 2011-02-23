#include "Basis.hh"
#include "io/manipulators.hh"

std::ostream& Basis::print(std::ostream& os) const
{
	os << "Basis (\n" << indent;
	for (BasisFunList::const_iterator it = _funs.begin();
		it != _funs.end(); ++it)
	{
		os << **it << "\n";
	}
	os << dedent << ")";
	return os;
}

void Basis::calcOverlap() const
{
	int n = size();
	_overlap.resize(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			_overlap(i, j) = _overlap(j, i)
				= ::overlap(*_funs[i], *_funs[j]);
		}
		_overlap(i, i) = ::overlap(*_funs[i], *_funs[i]);
	}
	_status.set(OVERLAP_CURRENT);
}

void Basis::calcKinetic() const
{
	int n = size();
	_kinetic.resize(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			_kinetic(i, j) = _kinetic(j, i)
				= ::kineticEnergy(*_funs[i], *_funs[j]);
		}
		_kinetic(i, i) = ::kineticEnergy(*_funs[i], *_funs[i]);
	}
	_status.set(KINETIC_CURRENT);
}

void Basis::calcOneElectron() const
{
	int n = size();
	_overlap.resize(n, n);
	_kinetic.resize(n, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			::oneElectron(*_funs[i], *_funs[j],
				&_overlap(i,j), &_kinetic(i,j));
			_overlap(j,i) = _overlap(i,j);
			_kinetic(j,i) = _kinetic(i,j);
		}
		::oneElectron(*_funs[i], *_funs[i],
			&_overlap(i,i), &_kinetic(i,i));
	}
	_status.set(OVERLAP_CURRENT);
	_status.set(KINETIC_CURRENT);
}
