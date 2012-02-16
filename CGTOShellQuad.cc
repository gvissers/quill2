#include "CGTOShellQuad.hh"

void CGTOShellQuad::setDPQ() const
{
	_dPQ.reset(new Eigen::ArrayXXd(_pAB.size(), 3*_pCD.size()));
	_dPQ->block(0, 0, _pAB.size(), _pCD.size())
		= Q(0).replicate(_pAB.size(), 1).colwise() - P(0);
	_dPQ->block(0, _pCD.size(), _pAB.size(), _pCD.size())
		= Q(1).replicate(_pAB.size(), 1).colwise() - P(1);
	_dPQ->block(0, 2*_pCD.size(), _pAB.size(), _pCD.size())
		= Q(2).replicate(_pAB.size(), 1).colwise() - P(2);
}

