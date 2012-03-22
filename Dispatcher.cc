#include "Dispatcher.hh"
#include "AbstractBF.hh"
#include "limits.hh"
#include "exceptions.hh"

STATIC_SINGLETON_OBJECT(Dispatcher);

Dispatcher::Dispatcher(): Li::Singleton<Dispatcher, true>(),
	_hasher(), _pair_funs()
{
	size_t id = classID<CGTO>();
	setPairCreator(id, id, CGTOPair::create);

	id = classID<CGTOPair>();
	setQuadCreator(id, id, CGTOQuad::create);
}

AbstractBFPair* Dispatcher::pair(const AbstractBF& f, const AbstractBF& g) const
{
	PairMap::const_iterator pit = _pair_funs.find(std::make_pair(f.cid, g.cid));
#ifdef DEBUG
	if (pit == _pair_funs.end())
		throw InvalidIndex(f.cid, g.cid);
#endif
	return pit->second(f, g);
}

AbstractBFQuad* Dispatcher::quad(const AbstractBFPair& p,
	const AbstractBFPair& q, BFQuadPool& pool) const
{
	QuadMap::const_iterator pit = _quad_funs.find(std::make_pair(p.cid, q.cid));
#ifdef DEBUG
	if (pit == _quad_funs.end())
		throw InvalidIndex(p.cid, q.cid);
#endif
	return pit->second(p, q, pool);
}

