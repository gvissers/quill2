#include "Basis.hh"
#include "Dispatcher.hh"
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

void Basis::setPairs() const
{
	const Dispatcher& dispatcher = Dispatcher::singleton();

	_pairs.clear();
	for (BasisFunList::const_iterator iit = _funs.begin(); iit != _funs.end(); ++iit)
	{
		for (BasisFunList::const_iterator jit = _funs.begin(); jit <= iit; ++jit)
			_pairs.push_back(dispatcher.pair(**iit, **jit));
	}

	_status.set(PAIRS_CURRENT);
}

void Basis::calcOverlap() const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	PairList::const_iterator pit = _pairs.begin();
	int n = size();
	_overlap.resize(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j, ++pit)
			_overlap(i, j) = _overlap(j, i) = (*pit)->overlap();
		_overlap(i, i) = (*pit)->overlap();
		++pit;
	}

	_status.set(OVERLAP_CURRENT);
}

void Basis::calcKinetic() const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	PairList::const_iterator pit = _pairs.begin();
	int n = size();
	_kinetic.resize(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j, ++pit)
			_kinetic(i, j) = _kinetic(j, i) = (*pit)->kineticEnergy(); 
		_kinetic(i, i) = (*pit)->kineticEnergy();
		++pit;
	}

	_status.set(KINETIC_CURRENT);
}

void Basis::calcOneElectron() const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	PairList::const_iterator pit = _pairs.begin();
	int n = size();
	_overlap.resize(n, n);
	_kinetic.resize(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j, ++pit)
		{
			(*pit)->oneElectron(_overlap(i,j), _kinetic(i,j));
			_overlap(j,i) = _overlap(i,j);
			_kinetic(j,i) = _kinetic(i,j);
		}
		(*pit)->oneElectron(_overlap(i,i), _kinetic(i,i));
		++pit;
	}

	_status.set(OVERLAP_CURRENT);
	_status.set(KINETIC_CURRENT);
}
