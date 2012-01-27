#ifndef GTO_ELEC_REP_HH
#define GTO_ELEC_REP_HH

#include <Eigen/Core>

std::vector<Eigen::ArrayXXd> gto_elec_rep_primitive_generic_1d(
	int lA, int lB, int lC, int lD,
	const Eigen::ArrayXXd& zeta, const Eigen::ArrayXXd& eta,
	double xA, double xB, double xC, double xD,
	const Eigen::ArrayXXd& P, const Eigen::ArrayXXd& Q,
	const Eigen::ArrayXXd& dPQ);

#endif // GTO_ELEC_REP_HH