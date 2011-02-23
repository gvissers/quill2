#include "Dispatcher.hh"
#include "DispatcherFiller.hh"
#include "AbstractBF.hh"
#include "limits.hh"
#include "exceptions.hh"

STATIC_SINGLETON_OBJECT(Dispatcher);

Dispatcher::Dispatcher(): Li::Singleton<Dispatcher, true>(),
	_hasher(), _S_funs(), _T_funs(), _one_elec_funs()
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

double Dispatcher::kineticEnergy(const AbstractBF& f, const AbstractBF& g) const
{
	KineticMap::const_iterator fit = _T_funs.find(std::make_pair(f.cid, g.cid));
#ifdef DEBUG
	if (fit == _T_funs.end())
		throw InvalidIndex(f.cid, g.cid);
#endif
	return fit->second(f, g);
}

void Dispatcher::oneElectron(const AbstractBF& f, const AbstractBF& g,
	double *S, double *T) const
{
	OneElecMap::const_iterator fit = _one_elec_funs.find(std::make_pair(f.cid, g.cid));
#ifdef DEBUG
	if (fit == _one_elec_funs.end())
		throw InvalidIndex(f.cid, g.cid);
#endif
	fit->second(f, g, S, T);
}
