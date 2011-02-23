#include "Dispatcher.hh"
#include "DispatcherFiller.hh"
#include "AbstractBF.hh"
#include "limits.hh"
#include "exceptions.hh"

STATIC_SINGLETON_OBJECT(Dispatcher);

Dispatcher::Dispatcher(): Li::Singleton<Dispatcher, true>(),
	_hasher(), _pair_funs()
{
	PairMapFiller<Limits::lmax, Limits::lmax>::fill();
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
