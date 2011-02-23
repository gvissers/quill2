#ifndef CGTOSPECPAIR_HH
#define CGTOSPECPAIR_HH

#include <Eigen/Dense>
#include "CGTOPair.hh"
#include "CGTOSpec.hh"
#include "gaussint/gto_one_elec.hh"
#include "constants.hh"

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
class CGTOSpecPair: public CGTOPair
{
	public:
		CGTOSpecPair(const CGTOSpec<lx1, ly1, lz1>& f,
			const CGTOSpec<lx2, ly2, lz2>& g):
			CGTOPair(f, g) {}

		const CGTOSpec<lx1, ly1, lz1>& f() const
		{
			return static_cast< const CGTOSpec<lx1, ly1, lz1>& >(AbstractBFPair::f());
		}
		const CGTOSpec<lx1, ly1, lz1>& g() const
		{
			return static_cast< const CGTOSpec<lx1, ly1, lz1>& >(AbstractBFPair::g());
		}

		double overlap() const;
		double kineticEnergy() const;
		void oneElectron(double& S, double& T) const;

		static AbstractBFPair *create(const AbstractBF& f, const AbstractBF& g)
		{
			const CGTOSpec<lx1, ly1, lz1>& ff = dynamic_cast< const CGTOSpec<lx1, ly1, lz1>& >(f);
			const CGTOSpec<lx2, ly2, lz2>& gg = dynamic_cast< const CGTOSpec<lx2, ly2, lz2>& >(g);
			return new CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>(ff, gg);
		}
};

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
double CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* generic_primitive_overlap(f().ls(), g().ls(), alpha(), beta(), asum(), exp_ared(), r()).matrix()
		* g().weights();
}

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
double CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* generic_primitive_kinetic(f().ls(), g().ls(), alpha(), beta(), asum(), ared(), exp_ared(), r()).matrix()
		* g().weights();
}

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
void CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::oneElectron(double &S, double &T) const
{
	Eigen::ArrayXXd Sp, Tp;
	generic_primitive_one_elec(f().ls(), g().ls(), alpha(), beta(), asum(), ared(), exp_ared(), r(), Sp, Tp);
	S = Constants::pi_sqrt_pi * f().weights().transpose() * Sp.matrix() * g().weights();
	T = Constants::pi_sqrt_pi * f().weights().transpose() * Tp.matrix() * g().weights();
}

#endif // CGTOSPECPAIR_HH
