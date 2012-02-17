#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 0>::electronRepulsion() const
{
	return weightsAB().transpose() * (KKW() * Fm(0)).matrix() * weightsCD();
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

	return weightsAB().transpose() * (KKW() * Fm(1) * dPW(i)).matrix()
		* weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd F1 = Fm(1);
	auto F0 = expmT() + 2*T()*F1;

	int i;
	for (i = 0; i < 3; i++)
	{
		if (lCD(i) == 1)
			break;
	}

	double x = lC(i) > 0 ? centerC(i) : centerD(i);
	F1 = F1*dQW(i) + F0.rowwise()*dxQ(i, x);
	return weightsAB().transpose() * (KKW() * F1).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 2, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd F2 = Fm(2);

	int i;
	for (i = 0; i < 3; ++i)
		if (lAB(i) > 0) break;

	auto dPWi = dPW(i);
	if (lAB(i) == 1)
	{
		// (s,s|d_ij,s)
		int j = i+1;
		if (lAB(j) != 1) ++j;

		F2 *= dPWi * dPW(j);
	}
	else
	{
		// (s,s|p,p) and (ss,d_i^2s)
		const Eigen::ArrayXXd& T = this->T();
		const Eigen::ArrayXXd& expmT = this->expmT();
		auto F1 = (1.0/3) * (expmT + 2*T*F2);
		auto F0 = expmT + 2*T*F1;

		F2 = F2*dPWi.square() - F1 * rho1()
			+ F0.colwise() * hInvWidthsAB();
	}

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 1, 1>::electronRepulsion() const
{
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
	Eigen::ArrayXXd F2 = Fm(2);
	auto F1 = (1.0/3) * (expmT + 2*T*F2);

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

		F2 = F2 * dQW(i) * dPWj + (F1 * dPWj).rowwise() * dxQ(i, xC);
	}
	else
	{
		// (p_i,s|p_i,s)
		double xC = lC(i) > 0 ? centerC(i) : centerD(i);
		auto dPWi = dPW(i);

		F2 = F2 * dQW(i) * dPWi
			+ F1 * (dPWi.rowwise()*dxQ(i, xC) + 0.5*invWidthsSum());
	}

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 2>::electronRepulsion() const
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

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}

