#include <Exception.hh>
#include "CGTOShellQuad.hh"
#include "boys.hh"

CGTOShellQuad::CGTOShellQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD):
	_pAB(pAB), _pCD(pCD),
	_inv_widths_sum((widthsAB().replicate(1, pCD.size()).rowwise()
		+ widthsCD()).inverse()),
	_dPQ(_pAB.size(), 3*_pCD.size()), _m(-1), _Fm()
{
	if (_pAB.samePositionId())
	{
		if (_pCD.samePositionId())
		{
			_pos_sym = _pAB.positionIdA() == _pCD.positionIdA()
				? POS_SYM_AAAA : POS_SYM_AACC;
		}
		else
		{
			_pos_sym = POS_SYM_AACD;
		}
	}
	else if (_pCD.samePositionId())
	{
		_pos_sym = POS_SYM_ABCC;
	}
	else
	{
		_pos_sym = POS_SYM_ABCD;
	}

	_dPQ.block(0, 0, _pAB.size(), _pCD.size())
		= Q(0).replicate(_pAB.size(), 1).colwise() - P(0);
	_dPQ.block(0, _pCD.size(), _pAB.size(), _pCD.size())
		= Q(1).replicate(_pAB.size(), 1).colwise() - P(1);
	_dPQ.block(0, 2*_pCD.size(), _pAB.size(), _pCD.size())
		= Q(2).replicate(_pAB.size(), 1).colwise() - P(2);
	
	_T = widthsReduced() * (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	_expmT = (-_T).exp();
}

Eigen::ArrayXXd CGTOShellQuad::Fm(int m) const
{
	if (m == _m)
		return _Fm;
	if (m > _m)
	{
		_Fm = _T.boys(m, _expmT);
		_m = m;
		return _Fm;
	}
	throw Li::Exception("m less than stored m");
}

