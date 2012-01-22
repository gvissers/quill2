#include <vector>
#include "gaussint/gto_nuc_attr.hh"
#include "constants.hh"
#include "limits.hh"
#include "gaussint/boys.hh"
#include <Exception.hh>

template<>
Eigen::ArrayXXd gto_nuc_attr_primitive_specialized<0, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::Vector3d& posA, const Eigen::Vector3d& posB,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge)
{
	int nr_nuc = nuc_pos.cols();
	int nr_rows = alpha.rows(), nr_cols = alpha.cols();
	Eigen::ArrayXXd U = Eigen::ArrayXXd::Zero(nr_rows, nr_cols*nr_nuc);
	Eigen::ArrayXXd dPC(nr_rows, nr_cols*nr_nuc);

	for (int i = 0; i < 3; i++)
	{
		Eigen::ArrayXXd P = (alpha*posA[i] + beta*posB[i]) / asum;
		for (int iC = 0; iC < nr_nuc; ++iC)
			dPC.block(0, iC*nr_cols, nr_rows, nr_cols) = P - nuc_pos(i, iC);
		U += dPC.square();
	}
	U *= asum.replicate(1, nr_nuc);

	Eigen::ArrayXXd Am = Fm(0, U);
	Eigen::ArrayXXd A = Eigen::ArrayXXd::Zero(nr_rows, nr_cols);
	for (int iC = 0; iC < nr_nuc; ++iC)
		A += -nuc_charge[iC] * Am.block(0, iC*nr_cols, nr_rows, nr_cols);
	
	return exp_ared * A / asum;
}

static std::vector<Eigen::ArrayXXd> gto_nuc_attr_primitive_generic_1d(
	int lA, int lB, double xA, double xB,
	const Eigen::ArrayXXd& theta,
	const Eigen::ArrayXXd& P, const Eigen::ArrayXXd& dPC)
{
	if (lA < lB)
	{
		std::swap(lA, lB);
		std::swap(xA, xB);
	}
	
	int rows = P.rows();
	int cols = dPC.cols();
	if (lA == 0)
		return std::vector<Eigen::ArrayXXd>(1, Eigen::ArrayXXd::Ones(rows, cols));

	int lsum = lA+lB;
	std::vector< std::vector<Eigen::ArrayXXd> > As(lsum+1);

	Eigen::ArrayXXd dPA = (P - xA).replicate(1, cols/P.cols());
		
	// A_0,0
	As[0].push_back(Eigen::ArrayXXd::Ones(rows, cols));
	// A_1,0
	As[1].push_back(dPA);
	As[1].push_back(-dPC);
	// A_a,0
	for (int iA = 1; iA < lsum; ++iA)
	{
		As[iA+1].push_back(dPA*As[iA][0] + 0.5*iA*As[iA-1][0]/theta);
		for (int m = 1; m < iA; ++m)
		{
			As[iA+1].push_back(dPA*As[iA][m] - dPC*As[iA][m-1]
				+ 0.5*iA*(As[iA-1][m]-As[iA-1][m-1])/theta);
		}
		As[iA+1].push_back(dPA*As[iA][iA] - dPC*As[iA][iA-1]
			- 0.5*iA*As[iA-1][iA-1]/theta);
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
	int nr_nuc = nuc_pos.cols();
	int nr_rows = alpha.rows(), nr_cols = alpha.cols();
	Eigen::ArrayXXd theta = asum.replicate(1, nr_nuc);
	Eigen::ArrayXXd U = Eigen::ArrayXXd::Zero(nr_rows, nr_cols*nr_nuc);
	Eigen::ArrayXXd dPC(nr_rows, nr_cols*nr_nuc);
	std::vector<Eigen::ArrayXXd> Axyz[3];

	for (int i = 0; i < 3; i++)
	{
		Eigen::ArrayXXd P = (alpha*posA[i] + beta*posB[i]) / asum;
		for (int iC = 0; iC < nr_nuc; ++iC)
			dPC.block(0, iC*nr_cols, nr_rows, nr_cols) = P - nuc_pos(i, iC);
		U += dPC.square();

		Axyz[i] = gto_nuc_attr_primitive_generic_1d(lsA[i], lsB[i],
			posA[i], posB[i], theta, P, dPC);
	}
	U *= theta;

	Eigen::Vector3i lsAB = lsA + lsB;
	int lsum = lsAB.sum();

	Eigen::ArrayXXd F = Fm(lsum, U);
	Eigen::ArrayXXd expmU = (-U).exp();
	Eigen::ArrayXXd Am = Axyz[0][lsAB.x()] * Axyz[1][lsAB.y()]
		* Axyz[2][lsAB.z()] * F;
	Eigen::ArrayXXd A = Eigen::ArrayXXd::Zero(nr_rows, nr_cols);
	for (int iC = 0; iC < nr_nuc; ++iC)
		A += -nuc_charge[iC] * Am.block(0, iC*nr_cols, nr_rows, nr_cols);

	for (int m = lsum-1; m >= 0; --m)
	{
		Am.setZero();
		for (int ix = 0; ix <= std::min(m, lsAB.x()); ++ix)
		{
			for (int iy = 0; iy <= std::min(m-ix, lsAB.y()); ++iy)
			{
				int iz = m - ix - iy;
				if (iz <= lsAB.z())
					Am += Axyz[0][ix] * Axyz[1][iy] * Axyz[2][iz];
			}
		}
		F = (expmU + 2*U*F) / (2*m+1);
		Am *= F;
		
		for (int iC = 0; iC < nr_nuc; ++iC)
			A += -nuc_charge[iC] * Am.block(0, iC*nr_cols, nr_rows, nr_cols);
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
