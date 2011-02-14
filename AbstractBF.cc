#include "AbstractBF.hh"
#include "Dispatcher.hh"

double overlap(const AbstractBF& f, const AbstractBF& g)
{
	return Dispatcher::singleton().overlap(f, g);
}
