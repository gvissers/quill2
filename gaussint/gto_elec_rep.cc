#include <vector>
#include "gaussint/gto_elec_rep.hh"

std::vector<Eigen::ArrayXXd> gto_elec_rep_primitive_generic_1d(
	int lA, int lB, int lC, int lD,
	const Eigen::ArrayXXd& zeta, const Eigen::ArrayXXd& eta,
	double xA, double xB, double xC, double xD,
	const Eigen::ArrayXXd& P, const Eigen::ArrayXXd& Q,
	const Eigen::ArrayXXd& dPQ)
{
	if (lA < lB)
	{
		std::swap(lA, lB);
		std::swap(xA, xB);
	}
	if (lC < lD)
	{
		std::swap(lC, lD);
		std::swap(xC, xD);
	}
	
	int l1 = lA+lB, l2 = lC+lD, lsum = l1+l2;
	int nr_rows = P.rows(), nr_cols = P.cols();
	if (lsum == 0)
		return std::vector<Eigen::ArrayXXd>(1, Eigen::ArrayXXd::Ones(nr_rows, nr_cols));

	Eigen::ArrayXXd asum = zeta + eta;
	Eigen::ArrayXXd inv_eta = 0.5 * eta.inverse();
	Eigen::ArrayXXd inv_zeta = 0.5 * zeta.inverse();
	Eigen::ArrayXXd inv_sum = 0.5 * asum.inverse();
	Eigen::ArrayXXd inv_ez = inv_eta + inv_zeta;
	Eigen::ArrayXXd rho1 = inv_zeta * inv_zeta / inv_ez;
	Eigen::ArrayXXd rho2 = inv_eta * inv_eta / inv_ez;
	Eigen::ArrayXXd ared = (zeta*eta) / asum;
	double dAB = xB - xA;
	double dCD = xD - xC;
	Eigen::ArrayXXd dAP = P - xA;
	Eigen::ArrayXXd dCQ = Q - xC;
	Eigen::ArrayXXd dPW = eta * dPQ / asum;
	Eigen::ArrayXXd dQW = -zeta * dPQ / asum;
	
	std::vector< std::vector< std::vector<Eigen::ArrayXXd> > > coefs(
		l1+1, std::vector< std::vector<Eigen::ArrayXXd> >(l2+1));
	// C_0,0,0,0
	coefs[0][0].push_back(Eigen::ArrayXXd::Ones(nr_rows, nr_cols));
	if (l2 > 0)
	{
		// C_0,0,1,0
		coefs[0][1].push_back(dCQ);
		coefs[0][1].push_back(dQW);
		// C_0,0,c,0
		for (int iC = 1; iC < l2; ++iC)
		{
			coefs[0][iC+1].push_back(dCQ*coefs[0][iC][0]
				+ 0.5*iC*coefs[0][iC-1][0]/eta);
			for (int m = 1; m < iC; ++m)
			{
				coefs[0][iC+1].push_back(dCQ*coefs[0][iC][m]
					+ iC*inv_eta*coefs[0][iC-1][m]
					+ dQW*coefs[0][iC][m-1]
					- iC*rho2*coefs[0][iC-1][m-1]);
			}
			coefs[0][iC+1].push_back(dCQ*coefs[0][iC][iC]
				+ dQW*coefs[0][iC][iC-1]
				- iC*rho2*coefs[0][iC-1][iC-1]);
			coefs[0][iC+1].push_back(dQW*coefs[0][iC][iC]);
		}
	}
	
	if (l1 > 0)
	{
		// C_1,0,c,0
		for (int iC = 0; iC <= l2; ++iC)
		{
			coefs[1][iC].push_back(dAP*coefs[0][iC][0]);
			for (int m = 1; m <= iC; ++m)
			{
				coefs[1][iC].push_back(dAP*coefs[0][iC][m]
					+ dPW*coefs[0][iC][m-1]
					+ iC*inv_sum*coefs[0][iC-1][m-1]);
			}
			coefs[1][iC].push_back(dPW*coefs[0][iC][iC]);
		}
		
		// C_a,0,c,0
		for (int iA = 1; iA < l1; ++iA)
		{
			coefs[iA+1][0].push_back(dAP*coefs[iA][0][0]
				+ iA*inv_zeta*coefs[iA-1][0][0]);
			for (int m = 1; m < iA; ++m)
			{
				coefs[iA+1][0].push_back(dAP*coefs[iA][0][m]
					+ iA*inv_zeta*coefs[iA-1][0][m]
					+ dPW*coefs[iA][0][m-1]
					- iA*rho1*coefs[iA-1][0][m-1]);
			}
			coefs[iA+1][0].push_back(dAP*coefs[iA][0][iA]
				+ dPW*coefs[iA][0][iA-1]
				- iA*rho1*coefs[iA-1][0][iA-1]);
			coefs[iA+1][0].push_back(dPW*coefs[iA][0][iA]);
			
			for (int iC = 1; iC <= l2; ++iC)
			{
				coefs[iA+1][iC].push_back(dAP*coefs[iA][iC][0]
					+ iA*coefs[iA-1][iC][0]);
				for (int m = 1; m < iA+iC; ++m)
				{
					coefs[iA+1][iC].push_back(dAP*coefs[iA][iC][m]
						+ iA*inv_zeta*coefs[iA-1][iC][m]
						+ dPW*coefs[iA][iC][m-1]
						- iA*rho1*coefs[iA-1][iC][m-1]
						+ iC*inv_sum*coefs[iA][iC-1][m-1]);
				}
				coefs[iA+1][iC].push_back(dAP*coefs[iA][iC][iA+iC]
					+ dPW*coefs[iA][iC][iA+iC-1]
					- iA*rho1*coefs[iA-1][iC][iA+iC-1]
					+ iC*inv_sum*coefs[iA][iC-1][iA+iC-1]);
				coefs[iA+1][iC].push_back(dPW*coefs[iA][iC][iA+iC]);
			}
		}
	}

	// C_a,0,c,d
	for (int iA = lA; iA <= l1; ++iA)
	{
		for (int iD = 1; iD <= lD; ++iD)
		{
			for (int iC = l2; iC >= lC+iD; --iC)
			{
				for (int m = 0; m < iA+iC; ++m)
					coefs[iA][iC][m] += dCD*coefs[iA][iC-1][m];
			}
		}
	}
	
	// C_a,b,c,d
	for (int iB = 1; iB <= lB; ++iB)
	{
		for (int iA = l1; iA >= lA+iB; --iA)
		{
			for (int m = 0; m < iA+l2; ++m)
				coefs[iA][l2][m] += dAB*coefs[iA-1][l2][m];
		}
	}
	
	return coefs[l1][l2];
}
