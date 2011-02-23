#ifndef ABSTRACTBFPAIR_HH
#define ABSTRACTBFPAIR_HH

#include "AbstractBF.hh"

class AbstractBFPair
{
	public:
		AbstractBFPair(const AbstractBF& f, const AbstractBF& g):
			_f(f), _g(g) {}

		const AbstractBF& f() const { return _f; }
		const AbstractBF& g() const { return _g; }

		virtual double overlap() const
		{
			double S, T;
			oneElectron(S, T);
			return S;
		}
		virtual double kineticEnergy() const
		{
			double S, T;
			oneElectron(S, T);
			return T;
		}
		virtual void oneElectron(double& S, double& T) const = 0;

	private:
		const AbstractBF& _f;
		const AbstractBF& _g;
};

#endif // ABSTRACTBFPAIR_HH
