#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 0, 0>::electronRepulsion() const
{
	return weightsAB().transpose() * (KKW() * Fm(0)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 1, 0>::electronRepulsion() const
{
	int i;
	for (i = 0; i < 3; i++)
		if (lAB(i) == 1) break;
		
	auto C = (KKW() * invWidthsSum()).rowwise()
		* widthsCD() * (q().centerA(i)-p().centerA(i));
	return weightsAB().transpose() * (C * Fm(1)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 0, 1>::electronRepulsion() const
{
	int i;
	for (i = 0; i < 3; i++)
		if (lCD(i) == 1) break;

	auto C = (KKW() * invWidthsSum()).colwise()
		* widthsAB() * (p().centerA(i)-q().centerA(i));
	return weightsAB().transpose() * (C * Fm(1)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 2, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd F2 = Fm(2);

	int i;
	for (i = 0; i < 3; ++i)
		if (lAB(i) > 0) break;

	if (lAB(i) == 1)
	{
		int j = i+1;
		if (lAB(j) != 1) ++j;

		double dd = (q().centerA(i) - p().centerA(i))
			* (q().centerA(j) - p().centerA(j));
		F2 *= (invWidthsSum().rowwise() * widthsCD()).square() * dd;
	}
	else
	{
		const Eigen::ArrayXXd& T = this->T();
		const Eigen::ArrayXXd& expmT = this->expmT();
		auto F1 = (1.0/3) * (expmT + 2*T*F2);
		auto F0 = expmT + 2*T*F1;
		double d = q().centerA(i) - p().centerA(i);
		F2 = F2 * (invWidthsSum().rowwise() * widthsCD()).square() * d * d
			- F1 * rho1()
			+ F0.colwise() * hInvWidthsAB();
	}

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 1, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd F2 = Fm(2);

	int i;
	for (i = 0; i < 3; ++i)
		if (lsum(i) > 0) break;

	if (lsum(i) == 1)
	{
		int j = i+1;
		if (lsum(j) != 1) ++j;
		double dd = -(q().centerA(i) - p().centerA(i))
			* (q().centerA(j) - p().centerA(j));
		F2 *= (invWidthsSum().square().colwise() * widthsAB()).rowwise()
			* widthsCD() * dd;
	}
	else
	{
		auto F1 = (1.0/3) * (expmT() + 2*T()*F2);
		double d = q().centerA(i) - p().centerA(i);
		F2 = ((F2 * invWidthsSum().square()).colwise()
			* widthsAB()).rowwise() * widthsCD() * (-d * d)
			+ 0.5 * F1 * invWidthsSum();
	}

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 0, 2>::electronRepulsion() const
{
	Eigen::ArrayXXd F2 = Fm(2);

	int i;
	for (i = 0; i < 3; i++)
		if (lCD(i) > 0) break;

	if (lCD(i) == 1)
	{
		int j = i+1;
		if (lCD(j) != 1) ++j;

		double dd = (q().centerA(i) - p().centerA(i))
			* (q().centerA(j) - p().centerA(j));
		F2 *= (invWidthsSum().colwise() * widthsAB()).square() * dd;
	}
	else
	{
		const Eigen::ArrayXXd& T = this->T();
		const Eigen::ArrayXXd& expmT = this->expmT();
		auto F1 = (1.0/3) * (expmT + 2*T*F2);
		auto F0 = expmT + 2*T*F1;
		double d = q().centerA(i) - p().centerA(i);
		F2 = F2 * (invWidthsSum().colwise() * widthsAB()).square() * d * d
			- F1 * rho2()
			+ F0.rowwise() * hInvWidthsCD();
	}

	return weightsAB().transpose() * (KKW() * F2).matrix() * weightsCD();
}