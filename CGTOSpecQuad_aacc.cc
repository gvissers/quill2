#include "CGTOSpecQuad.hh"
#include "boys.hh"

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_AACC, 0, 0>::electronRepulsion() const
{
	auto T = widthsReduced()
		* (p().centerA() - q().centerA()).squaredNorm();
	return weightsAB().transpose() * (KKW() * T.boys(0)).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_AACC, 1, 0>::electronRepulsion() const
{
	int i;
	for (i = 0; i < 3; i++)
		if (lAB(i) == 1) break;
		
	auto T = widthsReduced()
		* (p().centerA() - q().centerA()).squaredNorm();
	auto C = (KKW() * T.boys(1) * invWidthsSum()).rowwise()
		* widthsCD() * (q().centerA(i)-p().centerA(i));
	return weightsAB().transpose() * C.matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_AACC, 0, 1>::electronRepulsion() const
{
	int i;
	for (i = 0; i < 3; i++)
		if (lCD(i) == 1) break;

	auto T = widthsReduced()
		* (p().centerA() - q().centerA()).squaredNorm();
	auto C = (KKW() * T.boys(1) * invWidthsSum()).colwise()
		* widthsAB() * (p().centerA(i)-q().centerA(i));
	return weightsAB().transpose() * C.matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_AACC, 2, 0>::electronRepulsion() const
{
	Eigen::ArrayXXd T = (p().centerA() - q().centerA()).squaredNorm()
		* widthsReduced();
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);

	int i;
	for (i = 0; i < 3; ++i)
		if (lAB(i) > 0) break;

	if (lAB(i) == 1)
	{
		int j = i+1;
		if (lAB(j) != 1) ++j;

		double dd = (q().centerA(i) - p().centerA(i))
			* (q().centerA(j) - p().centerA(j));
		F *= (invWidthsSum().rowwise() * widthsCD()).square() * dd;
	}
	else
	{
		double d = q().centerA(i) - p().centerA(i);
		Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);
		F *= (invWidthsSum().rowwise() * widthsCD()).square() * d * d;
		F -= Fc * rho1();
		Fc = expmT + 2*T*Fc;
		F += Fc.colwise() * hInvWidthsAB();
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_AACC, 1, 1>::electronRepulsion() const
{
	Eigen::ArrayXXd T = (p().centerA() - q().centerA()).squaredNorm()
		* widthsReduced();
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);

	int i;
	for (i = 0; i < 3; ++i)
		if (lsum(i) > 0) break;

	if (lsum(i) == 1)
	{
		int j = i+1;
		if (lsum(j) != 1) ++j;
		double dd = -(q().centerA(i) - p().centerA(i))
			* (q().centerA(j) - p().centerA(j));
		F *= (invWidthsSum().square().colwise() * widthsAB()).rowwise()
			* widthsCD() * dd;
	}
	else
	{
		double d = q().centerA(i) - p().centerA(i);
		Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);
		F *= (invWidthsSum().square().colwise() * widthsAB()).rowwise()
			* widthsCD() * (-d * d);
		F += 0.5 * Fc * invWidthsSum();
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOQuad::POS_SYM_AACC, 0, 2>::electronRepulsion() const
{
	Eigen::ArrayXXd T = (p().centerA() - q().centerA()).squaredNorm()
		* widthsReduced();
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd F = T.boys(2, expmT);

	int i;
	for (i = 0; i < 3; i++)
		if (lCD(i) > 0) break;

	if (lCD(i) == 1)
	{
		int j = i+1;
		if (lCD(j) != 1) ++j;

		double dd = (q().centerA(i) - p().centerA(i))
			* (q().centerA(j) - p().centerA(j));
		F *= (invWidthsSum().colwise() * widthsAB()).square() * dd;
	}
	else
	{
		double d = q().centerA(i) - p().centerA(i);
		Eigen::ArrayXXd Fc = (1.0/3) * (expmT + 2*T*F);
		F *= (invWidthsSum().colwise() * widthsAB()).square() * d * d;
		F -= Fc * rho2();
		Fc = expmT + 2*T*Fc;
		F += Fc.rowwise() * hInvWidthsCD();
	}

	return weightsAB().transpose() * (KKW() * F).matrix() * weightsCD();
}