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
