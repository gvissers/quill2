#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 0>::electronRepulsion() const
{
	auto T = widthsReduced().rowwise() *
		((Q(0) - p().centerA(0)).square()
		+ (Q(1) - p().centerA(1)).square()
		+ (Q(2) - p().centerA(2)).square());
	return weightsAB().transpose() * (KKW() * T.boys(0)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 1, 0>::electronRepulsion() const
{
	int i;
	for (i = 0; i < 3; i++)
	{
		if (lAB(i) == 1)
			break;
	}

	auto T = widthsReduced().rowwise() *
		((Q(0) - p().centerA(0)).square()
		+ (Q(1) - p().centerA(1)).square()
		+ (Q(2) - p().centerA(2)).square());
	return weightsAB().transpose() * (KKW() * T.boys(1) * dPW(i)).matrix()
		* weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 1>::electronRepulsion() const
{
	//return CGTOQuad::electronRepulsion();
	Eigen::ArrayXXd T = widthsReduced().rowwise() *
		((Q(0) - p().centerA(0)).square()
		+ (Q(1) - p().centerA(1)).square()
		+ (Q(2) - p().centerA(2)).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(1, expmT);

	int i;
	for (i = 0; i < 3; i++)
	{
		if (lCD(i) == 1)
			break;
	}

	double x = lC(i) > 0 ? centerC(i) : centerD(i);
	F = F*dQW(i) + (expmT + 2*T*F).rowwise()*dxQ(i, x);
	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 2, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced().rowwise() *
		((Q(0) - p().centerA(0)).square()
		+ (Q(1) - p().centerA(1)).square()
		+ (Q(2) - p().centerA(2)).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);

	int i;
	for (i = 0; i < 3; ++i)
		if (lAB(i) > 0) break;

	auto dPWi = dPW(i);
	if (lAB(i) == 1)
	{
		// (s,s|d_ij,s)
		int j = i+1;
		if (lAB(j) != 1) ++j;

		F *= dPWi * dPW(j);
	}
	else
	{
		// (s,s|p,p) and (ss,d_i^2s)
		Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);

		F *= dPWi.square();
		F -= Fc * rho1();
		Fc = expmT + 2*T*Fc;
		F += Fc.colwise() * hInvWidthsAB();
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 1, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced().rowwise() *
		((Q(0) - p().centerA(0)).square()
		+ (Q(1) - p().centerA(1)).square()
		+ (Q(2) - p().centerA(2)).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);
	Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);

	int i;
	for (i = 0; i < 3; ++i)
		if (lsum(i) > 0) break;

	if (lsum(i) == 1)
	{
		// (p_i,s|p_j,s)
		int j = i+1;
		if (lsum(j) != 1) ++j;
		if (lCD(i) != 1)
			std::swap(i, j);
		double xC = lC(i) > 0 ? centerC(i) : centerD(i);
		auto dPWj = dPW(j);

		F *= dQW(i) * dPWj;
		F += (Fc * dPWj).rowwise() * dxQ(i, xC);
	}
	else
	{
		// (p_i,s|p_i,s)
		double xC = lC(i) > 0 ? centerC(i) : centerD(i);
		auto dPWi = dPW(i);

		F *= dQW(i) * dPWi;
		F += Fc * (dPWi.rowwise()*dxQ(i, xC) + 0.5*invWidthsSum());
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 2>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced().rowwise() *
		((Q(0) - p().centerA(0)).square()
		+ (Q(1) - p().centerA(1)).square()
		+ (Q(2) - p().centerA(2)).square());
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

