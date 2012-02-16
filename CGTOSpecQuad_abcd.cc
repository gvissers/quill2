#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 0>::electronRepulsion() const
{
	auto T = widthsReduced()
		* ((Q(0).replicate(p().size(), 1).colwise() - P(0)).square()
			+ (Q(1).replicate(p().size(), 1).colwise() - P(1)).square()
			+ (Q(2).replicate(p().size(), 1).colwise() - P(2)).square());
	return weightsAB().transpose() * (KKW() * T.boys(0)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 1, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced()
		* (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(1, expmT);
	Eigen::ArrayXXd Fc = expmT + 2*T*F;

	int i;
	for (i = 0; i < 3; i++)
	{
		if (lAB(i) == 1)
			break;
	}

	double x = lA(i) == 1 ? centerA(i) : centerB(i);
	F *= dPW(i);
	F += Fc.colwise() * dxP(i, x);

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced()
		* (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(1, expmT);
	Eigen::ArrayXXd Fc = expmT + 2*T*F;

	int i;
	for (i = 0; i < 3; i++)
	{
		if (lCD(i) == 1)
			break;
	}

	double x = lC(i) == 1 ? centerC(i) : centerD(i);
	F *= dQW(i);
	F += Fc.rowwise() * dxQ(i, x);

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 2, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced()
		* (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);
	Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);

	int i;
	for (i = 0; i < 3; ++i)
		if (lAB(i) > 0) break;

	double xA = centerA(i), xB = centerB(i);
	auto dPWi = dPW(i);
	if (lAB(i) == 1)
	{
		// (d_ij,s|s,s)
		int j = i+1;
		if (lAB(j) != 1) ++j;

		auto dAPi = dxP(i, lA(i) > 0 ? xA : xB);
		auto dAPj = dxP(j, lA(j) > 0 ? centerA(j) : centerB(j));
		auto dPWj = dPW(j);

		F *= dPWi * dPWj;
		F += Fc * (dPWi.colwise() * dAPj + dPWj.colwise() * dAPi);
		Fc = expmT + 2*T*Fc;
		F += Fc.colwise() * (dAPi * dAPj);
	}
	else if (lA(i) == 1)
	{
		// (p,p|s,s)
		double dAB = xB - xA;
		auto dAPi = dxP(i, xA);

		F *= dPWi.square();
		F += Fc * (dPWi.colwise() * (2*dAPi - dAB) - rho1());
		Fc = expmT + 2*T*Fc;
		F += Fc.colwise() * (dAPi.square() + hInvWidthsAB() - dAB*dAPi);
	}
	else
	{
		// (d_i^2,s|s,s)
		auto dAPi = dxP(i, lA(i) > 0 ? xA : xB);

		F *= dPWi.square();
		F += Fc * (dPWi.colwise() * (2*dAPi) - rho1());
		Fc = expmT + 2*T*Fc;
		F += Fc.colwise() * (dAPi.square() + hInvWidthsAB());
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 1, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced()
		* (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);
	Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);

	int i, j;
	for (i = 0; i < 3; ++i)
		if (lsum(i) > 0) break;

	if (lsum(i) == 2)
	{
		// (p_i,s|p_i,s)
		j = i;
	}
	else
	{
		// (p_i,s|p_j,s)
		j = i+1;
		if (lsum(j) != 1) ++j;
		if (lAB(i) != 1)
			std::swap(i, j);
	}
	
	double xA = lA(i) > 0 ? centerA(i) : centerB(i);
	double xC = lC(j) > 0 ? centerC(j) : centerD(j);
	auto dPWi = dPW(i);
	auto dQWj = dQW(j);
	auto dAPi = dxP(i, xA);
	auto dCQj = dxQ(j, xC);

	F *= dPWi * dQWj;
	if (i == j)
		F += Fc * (dQWj.colwise()*dAPi + dPWi.rowwise()*dCQj
			+ 0.5*invWidthsSum());
	else
		F += Fc * (dQWj.colwise()*dAPi + dPWi.rowwise()*dCQj);
	Fc = expmT + 2*T*Fc;
	F += (Fc.colwise() * dAPi).rowwise() * dCQj;

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 2>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced()
		* (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);
	Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);

	int i;
	for (i = 0; i < 3; ++i)
		if (lCD(i) > 0) break;

	double xC = centerC(i), xD = centerD(i);
	auto dQWi = dQW(i);
	if (lCD(i) == 1)
	{
		// (d_ij,s|s,s)
		int j = i+1;
		if (lCD(j) != 1) ++j;

		auto dCQi = dxQ(i, lC(i) > 0 ? xC : xD);
		auto dCQj = dxQ(j, lC(j) > 0 ? centerC(j) : centerD(j));
		auto dQWj = dQW(j);

		F *= dQWi * dQWj;
		F += Fc * (dQWi.rowwise() * dCQj + dQWj.rowwise() * dCQi);
		Fc = expmT + 2*T*Fc;
		F += Fc.rowwise() * (dCQi * dCQj);
	}
	else if (lC(i) == 1)
	{
		// (p,p|s,s)
		double dCD = xD - xC;
		auto dCQi = dxQ(i, xC);

		F *= dQWi.square();
		F += Fc * (dQWi.rowwise() * (2*dCQi - dCD) - rho2());
		Fc = expmT + 2*T*Fc;
		F += Fc.rowwise() * (dCQi.square() + hInvWidthsCD() - dCD*dCQi);
	}
	else
	{
		// (d_i^2,s|s,s)
		auto dCQi = dxQ(i, lC(i) > 0 ? xC : xD);

		F *= dQWi.square();
		F += Fc * (dQWi.rowwise() * (2*dCQi) - rho2());
		Fc = expmT + 2*T*Fc;
		F += Fc.rowwise() * (dCQi.square() + hInvWidthsCD());
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}
