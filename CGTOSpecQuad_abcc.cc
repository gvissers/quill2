#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_ABCC, 0, 0>::electronRepulsion() const
{
	auto T = widthsReduced().colwise() *
		((P(0) - q().centerA(0)).square()
		+ (P(1) - q().centerA(1)).square()
		+ (P(2) - q().centerA(2)).square());
	return weightsAB().transpose() * (KKW() * T.boys(0)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_ABCC, 1, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced().colwise() *
		((P(0) - q().centerA(0)).square()
		+ (P(1) - q().centerA(1)).square()
		+ (P(2) - q().centerA(2)).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(1, expmT);

	for (int i = 0; i < 3; i++)
	{
		if (lAB(i) == 1)
		{
			double x = lA(i) == 1 ? centerA(i) : centerB(i);
			F = F*dPW(i) + (expmT + 2*T*F).colwise()*dxP(i, x);
			break;
		}
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_ABCC, 0, 1>::electronRepulsion() const
{
	int i;
	for (i = 0; i < 3; i++)
	{
		if (lCD(i) == 1)
			break;
	}

	auto T = widthsReduced().colwise() *
		((P(0) - q().centerA(0)).square()
		+ (P(1) - q().centerA(1)).square()
		+ (P(2) - q().centerA(2)).square());
	return weightsAB().transpose() * (KKW() * T.boys(1) * dQW(i)).matrix()
		* weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_ABCC, 2, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced().colwise() *
		((P(0) - q().centerA(0)).square()
		+ (P(1) - q().centerA(1)).square()
		+ (P(2) - q().centerA(2)).square());
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
double CGTOSpecQuad<CGTOQuad::POS_SYM_ABCC, 1, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced().colwise() *
		((P(0) - q().centerA(0)).square()
		+ (P(1) - q().centerA(1)).square()
		+ (P(2) - q().centerA(2)).square());
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
		if (lAB(i) != 1)
			std::swap(i, j);
		double xA = lA(i) > 0 ? centerA(i) : centerB(i);
		auto dQWj = dQW(j);
		
		F *= dPW(i) * dQWj;
		F += (Fc * dQWj).colwise() * dxP(i, xA);
	}
	else
	{
		// (p_i,s|p_i,s)
		double xA = lA(i) > 0 ? centerA(i) : centerB(i);
		auto dQWi = dQW(i);

		F *= dPW(i) * dQWi;
		F += Fc * (dQWi.colwise()*dxP(i, xA) + 0.5*invWidthsSum());
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_ABCC, 0, 2>::electronRepulsion() const
{
	Eigen::ArrayXXd T = widthsReduced().colwise() *
		((P(0) - q().centerA(0)).square()
		+ (P(1) - q().centerA(1)).square()
		+ (P(2) - q().centerA(2)).square());
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);

	int i;
	for (i = 0; i < 3; ++i)
		if (lCD(i) > 0) break;

	auto dQWi = dQW(i);
	if (lCD(i) == 1)
	{
		// (s,s|d_ij,s)
		int j = i+1;
		if (lCD(j) != 1) ++j;

		F *= dQWi * dQW(j);
	}
	else
	{
		// (s,s|p,p) and (ss,d_i^2s)
		Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);

		F *= dQWi.square();
		F -= Fc * rho2();
		Fc = expmT + 2*T*Fc;
		F += Fc.rowwise() * hInvWidthsCD();
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}