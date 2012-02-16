#include "CGTOSpecQuad.hh"

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 0, 0>::electronRepulsion() const
{
	return weightsAB().transpose() * KKW().matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 1, 0>::electronRepulsion() const
{
	return 0;
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 0, 1>::electronRepulsion() const
{
	return 0;
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 2, 0>::electronRepulsion() const
{
	for (int i = 0; i < 3; i++)
	{
		int lsum = p().lsum(i) + q().lsum(i);
		if (lsum == 0)
			continue;
		if (lsum == 1)
			return 0;
		break;
	}

	auto A = ((-1.0/3)*rho1()).colwise() + hInvWidthsAB();
	return weightsAB().transpose() * (KKW() * A).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 1, 1>::electronRepulsion() const
{
	for (int i = 0; i < 3; i++)
	{
		int lsum = p().lsum(i) + q().lsum(i);
		if (lsum == 0)
			continue;
		if (lsum == 1)
			return 0;
		break;
	}

	auto A = (0.5/3)*invWidthsSum();
	return weightsAB().transpose() * (KKW() * A).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 0, 2>::electronRepulsion() const
{
	for (int i = 0; i < 3; i++)
	{
		int lsum = p().lsum(i) + q().lsum(i);
		if (lsum == 0)
			continue;
		if (lsum == 1)
			return 0;
		break;
	}

	auto A = ((-1.0/3)*rho2()).rowwise() + hInvWidthsCD();
	return weightsAB().transpose() * (KKW() * A).matrix() * weightsCD();
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 2, 1>::electronRepulsion() const
{
	return 0;
}

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 1, 2>::electronRepulsion() const
{
	return 0;
}
