#include "Dispatcher.hh"
#include "DispatcherFiller.hh"
#include "AbstractBF.hh"
#include "limits.hh"
#include "exceptions.hh"

STATIC_SINGLETON_OBJECT(Dispatcher);

Dispatcher::Dispatcher(): Li::Singleton<Dispatcher, true>(),
	_hasher(), _S_funs()
{
	PairMapFiller<Limits::lmax, Limits::lmax>::fill();
}

double Dispatcher::overlap(const AbstractBF& f, const AbstractBF& g) const
{
	OverlapMap::const_iterator fit = _S_funs.find(std::make_pair(f.cid, g.cid));
#ifdef DEBUG
	if (fit == _S_funs.end())
		throw InvalidIndex(f.cid, g.cid);
#endif
	return fit->second(f, g);
}
