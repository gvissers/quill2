#include "CGTOQuad.hh"
#include "boys.hh"

void CGTOQuad::elecRepPrim1d_psss(int i,
	const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
	Eigen::ArrayXXd& C0, Eigen::ArrayXXd& C1) const
{
#ifdef DEBUG
	if (lsum(i) != 1)
		throw Li::Exception("Total angular momentum should be 1");
#endif

	if (lAB(i) == 1)
	{
		double x = lA(i) == 1 ? centerA(i) : centerB(i);
		C0 = Pi - x;
		C1 = widthsCD() * (Qi - Pi) / widthsSum();
	}
	else
	{
		double x = lC(i) == 1 ? centerC(i) : centerD(i);
		C0 = Qi - x;
		C1 = -widthsAB() * (Qi - Pi) / widthsSum();
	}
}

void CGTOQuad::elecRepPrim1d_ppss(int i,
	const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
	Eigen::ArrayXXd& C0, Eigen::ArrayXXd& C1, Eigen::ArrayXXd& C2) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| ((lA(i) != 1 || lB(i) != 1) && (lC(i) != 1 || lD(i) != 1)))
		throw Li::Exception("angular momentum should be 1,1,0,0, or 0,0,1,1"); 
#endif
	Eigen::ArrayXXd inv_zeta = 0.5 * widthsAB().inverse();
	Eigen::ArrayXXd inv_eta = 0.5 * widthsCD().inverse();
	auto inv_ez = inv_eta + inv_zeta;

	if (lAB(i) > 0)
	{
		double xA = centerA(i), xB = centerB(i), dAB = xB - xA;
		Eigen::ArrayXXd dAP = Pi - xA;
		Eigen::ArrayXXd dPW = widthsCD() * (Qi - Pi) / widthsSum();
		auto rho1 = inv_zeta.square() / inv_ez;

		C0 = dAP.square() + inv_zeta - dAB*dAP;
		C1 = (2*dAP - dAB)*dPW - rho1;
		C2 = dPW.square();
	}
	else
	{
		double xC = centerC(i), xD = centerD(i), dCD = xD - xC;
		Eigen::ArrayXXd dCQ = Qi - xC;
		Eigen::ArrayXXd dQW = -widthsAB() * (Qi - Pi) / widthsSum();
		auto rho2 = inv_eta.square() / inv_ez;

		C0 = dCQ.square() + inv_eta - dCD*dCQ;
		C1 = (2*dCQ - dCD)*dQW - rho2;
		C2 = dQW.square();
	}
}

void CGTOQuad::elecRepPrim1d_psps(int i,
	const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
	Eigen::ArrayXXd& C0, Eigen::ArrayXXd& C1, Eigen::ArrayXXd& C2) const
{
	int lA = this->lA(i), lB = this->lB(i),
		lC = this->lC(i), lD = this->lD(i);
	int l1 = lA+lB, l2 = lC+lD;
#ifdef DEBUG
	if (l1 != 1 || l2 != 1)
		throw Li::Exception("Angular momentum should be 1 for both orbital pairs");
#endif

	double xA = lA == 1 ? centerA(i) : centerB(i);
	double xC = lC == 1 ? centerC(i) : centerD(i);

	Eigen::ArrayXXd asum = widthsSum();
	Eigen::ArrayXXd dAP = Pi - xA;
	Eigen::ArrayXXd dCQ = Qi - xC;
	Eigen::ArrayXXd dPQ = Qi - Pi;
	Eigen::ArrayXXd dPW = widthsCD() * dPQ / asum;
	Eigen::ArrayXXd dQW = -widthsAB() * dPQ / asum;

	C0 = dAP*dCQ;
	C1 = dAP*dQW + dPW*dCQ + 0.5*asum.inverse();
	C2 = dPW*dQW;
}

void CGTOQuad::elecRepPrim1d(int i,
	const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
	std::vector<Eigen::ArrayXXd>& Ci) const
{
	int lA = this->lA(i), lB = this->lB(i),
		lC = this->lC(i), lD = this->lD(i);
	int l1 = lA+lB, l2 = lC+lD, lsum = l1+l2;
	int nr_rows = p().size(), nr_cols = q().size();

	if (lsum == 0)
	{
		Ci.assign(1, Eigen::ArrayXXd::Ones(nr_rows, nr_cols));
		return;
	}
	else if (lsum == 1)
	{
		Ci.resize(2);
		elecRepPrim1d_psss(i, Pi, Qi, Ci[0], Ci[1]);
		return;
	}
	else if (lsum == 2)
	{
		if ((l1 == 2 && lA == 1) || (l2 == 2 && lC == 1))
		{
			Ci.resize(3);
			elecRepPrim1d_ppss(i, Pi, Qi, Ci[0], Ci[1], Ci[2]);
			return;
		}
		else if (l1 == 1 && l2 == 1)
		{
			Ci.resize(3);
			elecRepPrim1d_psps(i, Pi, Qi, Ci[0], Ci[1], Ci[2]);
			return;
		}
	}
	
	double xA = centerA(i), xB = centerB(i), xC = centerC(i), xD = centerD(i);
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

	Eigen::ArrayXXd asum = widthsSum();
	Eigen::ArrayXXd inv_zeta = 0.5 * widthsAB().inverse();
	Eigen::ArrayXXd inv_eta = 0.5 * widthsCD().inverse();
	Eigen::ArrayXXd inv_sum = 0.5 * asum.inverse();
	Eigen::ArrayXXd inv_ez = inv_eta + inv_zeta;
	Eigen::ArrayXXd rho1 = inv_zeta * inv_zeta / inv_ez;
	Eigen::ArrayXXd rho2 = inv_eta * inv_eta / inv_ez;
	double dAB = xB - xA;
	double dCD = xD - xC;
	Eigen::ArrayXXd dAP = Pi - xA;
	Eigen::ArrayXXd dCQ = Qi - xC;
	Eigen::ArrayXXd dPQ = Qi - Pi;
	Eigen::ArrayXXd dPW = widthsCD() * dPQ / asum;
	Eigen::ArrayXXd dQW = -widthsAB() * dPQ / asum;
	
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
				+ iC*inv_eta*coefs[0][iC-1][0]);
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
					+ iA*inv_zeta*coefs[iA-1][iC][0]);
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
					coefs[iA][iC][m] -= dCD*coefs[iA][iC-1][m];
			}
		}
	}
	
	// C_a,b,c,d
	for (int iB = 1; iB <= lB; ++iB)
	{
		for (int iA = l1; iA >= lA+iB; --iA)
		{
			for (int m = 0; m < iA+l2; ++m)
				coefs[iA][l2][m] -= dAB*coefs[iA-1][l2][m];
		}
	}
	
	Ci.swap(coefs[l1][l2]);
}

double CGTOQuad::electronRepulsion_ssss() const
{
	Eigen::ArrayXXd T = ((Q(0) - P(0)).square()
		+ (Q(1) - P(1)).square()
		+ (Q(2) - P(2)).square()) * widthsReduced();
	return weightsAB().transpose()
		* (KK() * Fm(0, T) / widthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion_psss(const Eigen::Vector3i& ls) const
{
	Eigen::ArrayXXd C0, C1;
	Eigen::ArrayXXd T = Eigen::ArrayXXd::Zero(p().size(), q().size());
	Eigen::ArrayXXd Pi, Qi;
	for (int i = 0; i < 3; i++)
	{
		Pi = P(i);
		Qi = Q(i);
		if (ls[i] == 1)
			elecRepPrim1d_psss(i, Pi, Qi, C0, C1);
		T += (Qi - Pi).square();
	}
	T *= widthsReduced();
	
	Eigen::ArrayXXd F = Fm(1, T);
	Eigen::ArrayXXd A = C1 * F + C0 * ((-T).exp() + 2*T*F);
	return weightsAB().transpose()
		* (KK() * A / widthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion() const
{
	const Eigen::Vector3i& lsA = p().f().ls();
	const Eigen::Vector3i& lsB = p().g().ls();
	const Eigen::Vector3i& lsC = q().f().ls();
	const Eigen::Vector3i& lsD = q().g().ls();
	Eigen::Vector3i ls = lsA + lsB + lsC + lsD;
	int lsum = ls.sum();

	if (lsum == 0)
		return electronRepulsion_ssss();
	else if (lsum == 1)
		return electronRepulsion_psss(ls);

	std::vector<Eigen::ArrayXXd> Axyz[3];
	Eigen::ArrayXXd T = Eigen::ArrayXXd::Zero(p().size(), q().size());
	Eigen::ArrayXXd Pi, Qi;
	for (int i = 0; i < 3; i++)
	{
		Pi = P(i);
		Qi = Q(i);
		elecRepPrim1d(i, Pi, Qi, Axyz[i]);
		T += (Qi - Pi).square();
	}
	T *= widthsReduced();

	Eigen::ArrayXXd F = Fm(lsum, T);
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd A = Axyz[0][ls.x()] * Axyz[1][ls.y()]
		* Axyz[2][ls.z()] * F;
	Eigen::ArrayXXd Am(p().size(), q().size());
	for (int m = lsum-1; m >= 0; --m)
	{
		Am.setZero();
		for (int ix = 0; ix <= std::min(m, ls.x()); ++ix)
		{
			for (int iy = 0; iy <= std::min(m-ix, ls.y()); ++iy)
			{
				int iz = m - ix - iy;
				if (iz <= ls.z())
					Am += Axyz[0][ix] * Axyz[1][iy] * Axyz[2][iz];
			}
		}
		F = (expmT + 2*T*F) / (2*m+1);
		A += Am * F;
	}

	return weightsAB().transpose()
		* (KK() * A / widthsSum().sqrt()).matrix()
		* weightsCD();
}