#include <vector>
#include "gaussint/gto_nuc_attr.hh"
#include "constants.hh"
#include "limits.hh"
#include "gaussint/boys.hh"
#include <Exception.hh>

static Eigen::ArrayXXd gto_nuc_attr_primitive_generic_1d(int l1, int l2,
	double x1, double x2,
	const Eigen::ArrayXXd& asum,
	const Eigen::ArrayXXd& P, double xnuc,
	const Eigen::ArrayXXd& U)
{
#ifdef DEBUG
	if (l1 < l2)
		throw Li::Exception("l1 must not be smaller than l2");
#endif
	
	if (l1 == 0)
		return Fm(0, U);

	int lsum = l1+l2;
	std::vector< std::vector<Eigen::ArrayXXd> > As(lsum+1);

	Eigen::ArrayXXd dpa = P - x1;
	Eigen::ArrayXXd dpc = P - xnuc;
	// A_0,0
	for (int m = 0; m <= lsum; ++m)
		As[0].push_back(Fm(m, U));
	// A_1,0
	for (int m = 0; m <= lsum-1; ++m)
		As[1].push_back(dpa*As[0][m] - dpc*As[0][m+1]);
	// A_a,0
	for (int i1 = 1; i1 < lsum; ++i1)
	{
		for (int m = 0; m <= lsum-i1-1; ++m)
		{
			As[i1+1].push_back(dpa*As[i1][m] - dpc*As[i1][m+1]
				+ 0.5*i1*(As[i1-1][m]-As[i1-1][m+1])/asum);
		}
	}
	
	if (l2 == 0)
		return As[l1][0];

	double dab = x2 - x1;
	for (int i2 = 1; i2 <= l2; i2++)
	{
		for (int i1 = lsum; i1 >= l1+i2; i1--)
			As[i1][0] -= dab * As[i1-1][0];
	}
	
	return As[lsum][0];
}

Eigen::ArrayXXd gto_nuc_attr_primitive_generic(
	const Eigen::Vector3i& ls1, const Eigen::Vector3i& ls2,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::Vector3d& pos1, const Eigen::Vector3d& pos2,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge)
{
	Eigen::ArrayXXd P[3];
	Eigen::ArrayXXd U;
	Eigen::ArrayXXd A = Eigen::ArrayXXd::Zero(alpha.rows(), alpha.cols());

	for (int i = 0; i < 3; i++)
		P[i] = (alpha*pos1[i] + beta*pos2[i]) / asum;

	for (int inuc = 0; inuc < nuc_pos.cols(); inuc++)
	{
		Eigen::Vector3d pos = nuc_pos.col(inuc);
		U = asum * ((P[0] - pos[0]).square()
			+ (P[1] - pos[1]).square()
			+ (P[2] - pos[2]).square());

		Eigen::ArrayXXd Ap = Eigen::ArrayXXd::Ones(alpha.rows(), alpha.cols());
		for (int i = 0; i < 3; i++)
		{
			if (ls1[i] < ls2[i])
				Ap *= gto_nuc_attr_primitive_generic_1d(
					ls2[i], ls1[i], pos2[i], pos1[i], asum,
					P[i], pos[i], U);
			else
				Ap *= gto_nuc_attr_primitive_generic_1d(
					ls1[i], ls2[i], pos1[i], pos2[i], asum,
					P[i], pos[i], U);
		}
		
		A += -nuc_charge[inuc] * Ap;
	}
	
	return M_2_SQRTPI * exp_ared * A / asum;
}

double gto_nuc_attr_generic(const Eigen::Vector3i& ls1,
	const Eigen::VectorXd& weights1, const Eigen::VectorXd& widths1,
	const Eigen::Vector3d& pos1,
	const Eigen::Vector3i& ls2,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& pos2,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge)
{
	int n1 = widths1.size(), n2 = widths2.size();
	Eigen::ArrayXXd alpha = widths1.replicate(1, n2);
	Eigen::ArrayXXd beta = widths2.transpose().replicate(n1, 1);
	Eigen::ArrayXXd asum = alpha + beta;
	Eigen::ArrayXXd ared = (widths1 * widths2.transpose()).array() / asum;
	Eigen::Vector3d r = pos1 - pos2;
	Eigen::ArrayXXd exp_ared = (-r.squaredNorm() * ared).exp();

	return Constants::pi_sqrt_pi * weights1.transpose()
		* gto_nuc_attr_primitive_generic(ls1, ls2, alpha, beta, pos1, pos2,
			asum, exp_ared, nuc_pos, nuc_charge).matrix()
		* weights2;
}
