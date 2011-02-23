#ifndef CGTOPAIR_HH
#define CGTOPAIR_HH

#include <Eigen/Dense>
#include "AbstractBFPair.hh"
#include "CGTO.hh"
#include "gaussint/gto_one_elec.hh"
#include "constants.hh"

class CGTOPair: public AbstractBFPair
{
	public:
		CGTOPair(const CGTO& f, const CGTO& g): AbstractBFPair(f, g) {}

		const CGTO& f() const
		{
			return static_cast< const CGTO& >(AbstractBFPair::f());
		}
		const CGTO& g() const
		{
			return static_cast< const CGTO& >(AbstractBFPair::g());
		}

		virtual double overlap() const;
		virtual double kineticEnergy() const;
		virtual void oneElectron(double& S, double& T) const;

		Eigen::ArrayXXd alpha() const
		{
			return f().widths().replicate(1, g().size());
		}
		Eigen::ArrayXXd beta() const
		{
			return g().widths().transpose().replicate(f().size(), 1);
		}
		Eigen::ArrayXXd asum() const
		{
			return alpha() + beta();
		}
		Eigen::ArrayXXd ared() const
		{
			return alpha()*beta() / asum();
		}
		Eigen::ArrayXXd exp_ared() const
		{
			return (-r().squaredNorm()*ared()).exp();
		}
		Eigen::Vector3d r() const
		{
			return g().center() - f().center();
		}

		static AbstractBFPair *create(const AbstractBF& f, const AbstractBF& g)
		{
			const CGTO& ff = dynamic_cast< const CGTO& >(f);
			const CGTO& gg = dynamic_cast< const CGTO& >(g);
			return new CGTOPair(ff, gg);
		}
};

#endif // CGTOPAIR_HH
