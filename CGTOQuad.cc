#include "CGTOQuad.hh"
#include "gaussint/gto_elec_rep.hh"
#include "gaussint/boys.hh"

double CGTOQuad::electronRepulsion() const
{
	const Eigen::Vector3i& lsA = p().f().ls();
	const Eigen::Vector3i& lsB = p().g().ls();
	const Eigen::Vector3i& lsC = q().f().ls();
	const Eigen::Vector3i& lsD = q().g().ls();

	Eigen::ArrayXXd T = Eigen::ArrayXXd::Zero(p().size(), q().size());
	std::vector<Eigen::ArrayXXd> Axyz[3];
	for (int i = 0; i < 3; i++)
	{
		Eigen::ArrayXXd Pi = P(i);
		Eigen::ArrayXXd Qi = Q(i);
		
		Axyz[i] = gto_elec_rep_primitive_generic_1d(
			lsA[i], lsB[i], lsC[i], lsD[i], widthsAB(), widthsCD(),
			p().f().center(i), p().g().center(i),
			q().f().center(i), q().g().center(i),
			Pi, Qi, Qi-Pi);
		T += (Qi-Pi).square();
	}
	T *= widthsReduced();

	Eigen::Vector3i ls = lsA + lsB + lsC + lsD;
	int lsum = ls.sum();

	Eigen::ArrayXXd F = Fm(lsum, T);
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd A = Axyz[0][ls.x()] * Axyz[1][ls.y()]
		* Axyz[2][ls.z()] * F;
	Eigen::ArrayXXd Am(p().size(), q().size());
	for (int m = lsum-1; m >= 0; --m)
	{
		Am.setZero();
		for (int ix = 0; ix <= std::min(m, ls.x()); ++ix)
		{
			for (int iy = 0; iy <= std::min(m-ix, ls.y()); ++iy)
			{
				int iz = m - ix - iy;
				if (iz <= ls.z())
					Am += Axyz[0][ix] * Axyz[1][iy] * Axyz[2][iz];
			}
		}
		F = (expmT + 2*T*F) / (2*m+1);
		A += Am * F;
	}

	return weightsAB().transpose() * (KK() * A / widthsSum().sqrt()).matrix() * weightsCD();
}