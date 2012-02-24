#include "FmCoefs.hh"

void FmCoefs::multiply(const EriCoefs& coefs, int l1, int l2)
{
	int lsum = l1+l2;
	for (int m = lsum+_mmax; m >= _mmin; --m)
	{
		// i >= 0
		// i <= lsum
		// m-i >= _mmin => i <= m-_mmin
		// m-i <= _mmax => i >= m-_mmax
		int i = std::max(0, m-_mmax);
		block(m) = coefs(l1, l2, i) * block(m-i);
		for (++i; i <= std::min(lsum, m-_mmin); ++i)
			block(m) += coefs(l1, l2, i) * block(m-i);
	}
	_mmax += lsum;
}
