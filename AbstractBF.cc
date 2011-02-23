#include "AbstractBF.hh"
#include "Dispatcher.hh"

double overlap(const AbstractBF& f, const AbstractBF& g)
{
	return Dispatcher::singleton().overlap(f, g);
}

double kineticEnergy(const AbstractBF& f, const AbstractBF& g)
{
	return Dispatcher::singleton().kineticEnergy(f, g);
}

void oneElectron(const AbstractBF& f, const AbstractBF& g, double *S, double *T)
{
	Dispatcher::singleton().oneElectron(f, g, S, T);
}
