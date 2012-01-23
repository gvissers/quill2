#ifndef ABSTRACTBFQUAD_HH
#define ABSTRACTBFQUAD_HH

#include "AbstractBFPair.hh"

class AbstractBFQuad
{
public:
	AbstractBFQuad(const AbstractBFPair& p,
		const AbstractBFPair& q): _p(p), _q(q) {}

	const AbstractBFPair& p() const { return _p; }
	const AbstractBFPair& q() const { return _q; }
	
private:
	const AbstractBFPair& _p;
	const AbstractBFPair& _q;
};

#endif // ABSTRACTBFQUAD_HH