#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 0, 0>::electronRepulsion() const
{
	return weightsAB().transpose() * (KKW() * Fm(0)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 1, 0>::electronRepulsion() const
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
	F1 = F1*dPW(i) + F0.colwise()*dxP(i, x);

	return weightsAB().transpose() * (KKW() * F1).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 0, 1>::electronRepulsion() const
{
	int i;
	for (i = 0; i < 3; i++)
	{
		if (lCD(i) == 1)
			break;
	}

	return weightsAB().transpose() * (KKW() * Fm(1) * dQW(i)).matrix()
		* weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 2, 0>::electronRepulsion() const
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

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 1, 1>::electronRepulsion() const
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
		if (lAB(i) != 1)
			std::swap(i, j);
		double xA = lA(i) > 0 ? centerA(i) : centerB(i);
		auto dQWj = dQW(j);
		
		F2 = F2 * dPW(i) * dQWj
			+ (F1 * dQWj).colwise() * dxP(i, xA);
	}
	else
	{
		// (p_i,s|p_i,s)
		double xA = lA(i) > 0 ? centerA(i) : centerB(i);
		auto dQWi = dQW(i);

		F2 = F2 * dPW(i) * dQWi
			+ F1 * (dQWi.colwise()*dxP(i, xA) + 0.5*invWidthsSum());
	}

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 0, 2>::electronRepulsion() const
{
	Eigen::ArrayXXd F2 = Fm(2);

	int i;
	for (i = 0; i < 3; ++i)
		if (lCD(i) > 0) break;

	auto dQWi = dQW(i);
	if (lCD(i) == 1)
	{
		// (s,s|d_ij,s)
		int j = i+1;
		if (lCD(j) != 1) ++j;

		F2 *= dQWi * dQW(j);
	}
	else
	{
		// (s,s|p,p) and (ss,d_i^2s)
		const Eigen::ArrayXXd& T = this->T();
		const Eigen::ArrayXXd& expmT = this->expmT();
		auto F1 = (1.0/3) * (expmT + 2*T*F2);
		auto F0 = expmT + 2*T*F1;

		F2 = F2 * dQWi.square() - F1 * rho2() + F0.rowwise() * hInvWidthsCD();
	}

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}