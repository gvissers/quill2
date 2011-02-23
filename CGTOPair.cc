#include "CGTOPair.hh"
#include "gaussint/gto_one_elec.hh"

double CGTOPair::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* gto_overlap_primitive_generic(f().ls(), g().ls(),
			alpha(), beta(), asum(), exp_ared(), r()).matrix()
		* g().weights();
}

double CGTOPair::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* gto_kinetic_primitive_generic(f().ls(), g().ls(),
			alpha(), beta(), asum(), ared(), exp_ared(), r()).matrix()
		* g().weights();
}

void CGTOPair::oneElectron(double &S, double &T) const
{
	Eigen::ArrayXXd Sp, Tp;
	gto_one_elec_primitive_generic(f().ls(), g().ls(), alpha(), beta(),
		asum(), ared(), exp_ared(), r(), Sp, Tp);
	S = Constants::pi_sqrt_pi * f().weights().transpose()
		* Sp.matrix() * g().weights();
	T = Constants::pi_sqrt_pi * f().weights().transpose()
		* Tp.matrix() * g().weights();
}

