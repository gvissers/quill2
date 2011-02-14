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

const Eigen::MatrixXd& Basis::overlap() const
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
	return _overlap;
}
