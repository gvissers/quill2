#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 0>::electronRepulsion() const
{
	return mulWeights(KKW() * Fm(0));
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 1, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd F1 = Fm(1);
	auto F0 = expmT() + 2*T()*F1;

	int i;
	for (i = 0; i < 3; i++)
	{
		if (lAB(i) == 1)
			break;
	}

	double x = lA(i) == 1 ? centerA(i) : centerB(i);
	F1 = F1*dPW(i) + F0.colwise() * dxP(i, x);

	return mulWeights(KKW() * F1);
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd F1 = Fm(1);
	auto F0 = expmT() + 2*T()*F1;

	int i;
	for (i = 0; i < 3; i++)
	{
		if (lCD(i) == 1)
			break;
	}

	double x = lC(i) == 1 ? centerC(i) : centerD(i);
	F1 = F1*dQW(i) + F0.rowwise() * dxQ(i, x);

	return mulWeights(KKW() * F1);
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 2, 0>::electronRepulsion() const
{
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
	Eigen::ArrayXXd F2 = Fm(2);
	Eigen::ArrayXXd F1 = (1.0/3) * (expmT + 2*T*F2);
	auto F0 = expmT + 2*T*F1;

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

		F2 *= dPWi * dPWj;
		F2 += F1 * (dPWi.colwise() * dAPj + dPWj.colwise() * dAPi);
		F2 += F0.colwise() * (dAPi * dAPj);
	}
	else if (lA(i) == 1)
	{
		// (p,p|s,s)
		double dAB = xB - xA;
		auto dAPi = dxP(i, xA);

		F2 *= dPWi.square();
		F2 += F1 * (dPWi.colwise() * (2*dAPi - dAB) - rho1());
		F2 += F0.colwise() * (dAPi.square() + hInvWidthsAB() - dAB*dAPi);
	}
	else
	{
		// (d_i^2,s|s,s)
		auto dAPi = dxP(i, lA(i) > 0 ? xA : xB);

		F2 *= dPWi.square();
		F2 += F1 * (dPWi.colwise() * (2*dAPi) - rho1());
		F2 += F0.colwise() * (dAPi.square() + hInvWidthsAB());
	}

	return mulWeights(KKW() * F2);
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 1, 1>::electronRepulsion() const
{
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
	Eigen::ArrayXXd F2 = Fm(2);
	Eigen::ArrayXXd F1 = (1.0/3) * (expmT + 2*T*F2);
	auto F0 = expmT + 2*T*F1;

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

	F2 *= dPWi * dQWj;
	if (i == j)
		F2 += F1 * (dQWj.colwise()*dAPi + dPWi.rowwise()*dCQj
			+ 0.5*invWidthsSum());
	else
		F2 += F1 * (dQWj.colwise()*dAPi + dPWi.rowwise()*dCQj);
	F2 += (F0.colwise() * dAPi).rowwise() * dCQj;

	return mulWeights(KKW() * F2);
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 2>::electronRepulsion() const
{
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
	Eigen::ArrayXXd F2 = Fm(2);
	Eigen::ArrayXXd F1 = (1.0/3) * (expmT + 2*T*F2);
	auto F0 = expmT + 2*T*F1;

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

		F2 *= dQWi * dQWj;
		F2 += F1 * (dQWi.rowwise() * dCQj + dQWj.rowwise() * dCQi);
		F2 += F0.rowwise() * (dCQi * dCQj);
	}
	else if (lC(i) == 1)
	{
		// (p,p|s,s)
		double dCD = xD - xC;
		auto dCQi = dxQ(i, xC);

		F2 *= dQWi.square();
		F2 += F1 * (dQWi.rowwise() * (2*dCQi - dCD) - rho2());
		F2 += F0.rowwise() * (dCQi.square() + hInvWidthsCD() - dCD*dCQi);
	}
	else
	{
		// (d_i^2,s|s,s)
		auto dCQi = dxQ(i, lC(i) > 0 ? xC : xD);

		F2 *= dQWi.square();
		F2 += F1 * (dQWi.rowwise() * (2*dCQi) - rho2());
		F2 += F0.rowwise() * (dCQi.square() + hInvWidthsCD());
	}

	return mulWeights(KKW() * F2);
}
