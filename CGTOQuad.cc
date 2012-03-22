#include "CGTOQuad.hh"
#include "BFQuadPool.hh"

static const CGTOPair& firstPair(const CGTOPair& pp, const CGTOPair& qq)
{
	return pp.ishellPair() < qq.ishellPair() ? qq : pp;
}

static const CGTOPair& secondPair(const CGTOPair& pp, const CGTOPair& qq)
{
	return pp.ishellPair() < qq.ishellPair() ? pp : qq;
}

CGTOQuad::CGTOQuad(const CGTOPair& pp, const CGTOPair& qq):
	AbstractBFQuad(firstPair(pp, qq), secondPair(pp, qq)),
	_ishell_quad(CGTOShellList::pairIndex(p().ishellPair(), q().ishellPair())),
	_shell_quad(CGTOShellList::singleton().quad(_ishell_quad)) {}

AbstractBFQuad* CGTOQuad::create(const AbstractBFPair& p,
	const AbstractBFPair& q, BFQuadPool& pool)
{
	try
	{
#ifdef DEBUG
		const CGTOPair& pp = dynamic_cast< const CGTOPair& >(p);
		const CGTOPair& qq = dynamic_cast< const CGTOPair& >(q);
#else
		const CGTOPair& pp = static_cast< const CGTOPair& >(p);
		const CGTOPair& qq = static_cast< const CGTOPair& >(q);
#endif
		return new(pool) CGTOQuad(pp, qq);
	}
	catch (const std::bad_cast&)
	{
		throw Li::Exception("Invalid basis function pair type");
	}
}
