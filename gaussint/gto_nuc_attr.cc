#include <vector>
#include "gaussint/gto_nuc_attr.hh"
#include "constants.hh"
#include "limits.hh"
#include "gaussint/boys.hh"
#include <Exception.hh>

static std::vector<Eigen::ArrayXXd> gto_nuc_attr_primitive_generic_1d(
	int lA, int lB, double xA, double xB,
	const Eigen::ArrayXXd& asum,
	const Eigen::ArrayXXd& P, double xnuc)
{
#ifdef DEBUG
	if (lA < lB)
		throw Li::Exception("lA must not be smaller than lB");
#endif
	
	int rows = P.rows();
	int cols = P.cols();
	if (lA == 0)
		return std::vector<Eigen::ArrayXXd>(1, Eigen::ArrayXXd::Ones(rows, cols));

	int lsum = lA+lB;
	std::vector< std::vector<Eigen::ArrayXXd> > As(lsum+1);

	Eigen::ArrayXXd dPA = P - xA;
	Eigen::ArrayXXd dPC = P - xnuc;
	// A_0,0
	As[0].push_back(Eigen::ArrayXXd::Ones(rows, cols));
	// A_1,0
	As[1].push_back(dPA);
	As[1].push_back(-dPC);
	// A_a,0
	for (int iA = 1; iA < lsum; ++iA)
	{
		As[iA+1].push_back(dPA*As[iA][0] + 0.5*iA*As[iA-1][0]/asum);
		for (int m = 1; m < iA; ++m)
		{
			As[iA+1].push_back(dPA*As[iA][m] - dPC*As[iA][m-1]
				+ 0.5*iA*(As[iA-1][m]-As[iA-1][m-1])/asum);
		}
		As[iA+1].push_back(dPA*As[iA][iA] - dPC*As[iA][iA-1]
			- 0.5*iA*As[iA-1][iA-1]/asum);
		As[iA+1].push_back(-dPC*As[iA][iA]);
	}
	
	double dAB = xB - xA;
	for (int iB = 1; iB <= lB; iB++)
	{
		for (int iA = lsum; iA >= lA+iB; iA--)
		{
			for (int m = 0; m < iA; m++)
				As[iA][m] -= dAB * As[iA-1][m];
		}
	}
	
	return As.back();
}

Eigen::ArrayXXd gto_nuc_attr_primitive_generic(
	const Eigen::Vector3i& lsA, const Eigen::Vector3i& lsB,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::Vector3d& posA, const Eigen::Vector3d& posB,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge)
{
	Eigen::ArrayXXd P[3];
	Eigen::ArrayXXd U;
	Eigen::ArrayXXd A = Eigen::ArrayXXd::Zero(alpha.rows(), alpha.cols());
	Eigen::Vector3i lsAB = lsA + lsB;
	int lsum = lsAB.sum();

	for (int i = 0; i < 3; i++)
		P[i] = (alpha*posA[i] + beta*posB[i]) / asum;

	//for (int inuc = 0; inuc < nuc_pos.cols(); inuc++)
	int inuc = 0;
	{
		Eigen::Vector3d pos = nuc_pos.col(inuc);
		U = asum * ((P[0] - pos[0]).square()
			+ (P[1] - pos[1]).square()
			+ (P[2] - pos[2]).square());

		std::vector<Eigen::ArrayXXd> Ax =
			lsB.x() <= lsA.x()
			? gto_nuc_attr_primitive_generic_1d(lsA.x(), lsB.x(),
				posA.x(), posB.x(), asum, P[0], pos.x())
			: gto_nuc_attr_primitive_generic_1d(lsB.x(), lsA.x(),
				posB.x(), posA.x(), asum, P[0], pos.x());
		std::vector<Eigen::ArrayXXd> Ay =
			lsB.y() <= lsA.y()
			? gto_nuc_attr_primitive_generic_1d(lsA.y(), lsB.y(),
				posA.y(), posB.y(), asum, P[1], pos.y())
			: gto_nuc_attr_primitive_generic_1d(lsB.y(), lsA.y(),
				posB.y(), posA.y(), asum, P[1], pos.y());
		std::vector<Eigen::ArrayXXd> Az =
			lsB.z() <= lsA.z()
			? gto_nuc_attr_primitive_generic_1d(lsA.z(), lsB.z(),
				posA.z(), posB.z(), asum, P[2], pos.z())
			: gto_nuc_attr_primitive_generic_1d(lsB.z(), lsA.z(),
				posB.z(), posA.z(), asum, P[2], pos.z());

		for (int m = 0; m <= lsum; ++m)
		{
			Eigen::ArrayXXd Am = Eigen::ArrayXXd::Zero(alpha.rows(), alpha.cols());
			for (int ix = 0; ix <= std::min(m, lsAB.x()); ++ix)
			{
				for (int iy = 0; iy <= std::min(m-ix, lsAB.y()); ++iy)
				{
					int iz = m - ix - iy;
					if (iz <= lsAB.z())
						Am += Ax[ix] * Ay[iy] * Az[iz];
				}
			}
			A += -nuc_charge[inuc] * Am * Fm(m, U);
		}
	}
	
	return exp_ared * A / asum;
}

double gto_nuc_attr_generic(const Eigen::Vector3i& lsA,
	const Eigen::VectorXd& weights1, const Eigen::VectorXd& widths1,
	const Eigen::Vector3d& posA,
	const Eigen::Vector3i& lsB,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& posB,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge)
{
	int n1 = widths1.size(), n2 = widths2.size();
	Eigen::ArrayXXd alpha = widths1.replicate(1, n2);
	Eigen::ArrayXXd beta = widths2.transpose().replicate(n1, 1);
	Eigen::ArrayXXd asum = alpha + beta;
	Eigen::ArrayXXd ared = (widths1 * widths2.transpose()).array() / asum;
	Eigen::Vector3d r = posA - posB;
	Eigen::ArrayXXd exp_ared = (-r.squaredNorm() * ared).exp();

	return 2 * M_PI * weights1.transpose()
		* gto_nuc_attr_primitive_generic(lsA, lsB, alpha, beta, posA, posB,
			asum, exp_ared, nuc_pos, nuc_charge).matrix()
		* weights2;
}
