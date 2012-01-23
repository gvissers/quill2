#include <vector>
#include <Exception.hh>
#include "gaussint/gto_overlap.hh"

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d&)
{
	return exp_ared / (asum * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * beta * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.y() * beta * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.z() * beta * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.x() * alpha * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 1, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.y() * alpha * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.z() * alpha * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<2, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.x()*r.x()*beta.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 1, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.y() * beta.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.z() * beta.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5 - r.x()*r.x()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.x() * r.y() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.x() * r.z() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 2, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.y()*r.y()*beta.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.y() * r.z() * beta.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.y() * r.x() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5 - r.y()*r.y()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.y() * r.z() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 2, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.z()*r.z()*beta.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.z() * r.x() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.z() * r.y() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5 - r.z()*r.z()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 2, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.x()*r.x()*alpha.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 1, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.y() * alpha.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 0, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.z() * alpha.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 2, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.y()*r.y()*alpha.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 1, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.y() * r.z() * alpha.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 2>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.z()*r.z()*alpha.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

Eigen::ArrayXXd gto_overlap_primitive_1d(int lA, int lB, double x,
	const Eigen::ArrayXXd& beta, const Eigen::ArrayXXd& asum)
{
#ifdef DEBUG
	if (lA < lB)
		throw Li::Exception("lA must not be smaller than lB");
#endif

	int powsum = lA + lB;
	if (powsum == 0)
		return Eigen::ArrayXXd::Ones(beta.rows(), beta.cols());
	if (powsum == 1)
		return x * beta / asum;

	std::vector<Eigen::ArrayXXd> coefs;
	coefs.push_back(x*beta/asum);
	// Recursion in lA:
	// S_lA+1,0 = [lA/2 S_lA-1,0 + x beta S_lA,0] / (alpha + beta)
	coefs.push_back((0.5 + x*beta*coefs[0]) / asum);
	for (int i = 2; i < powsum; i++)
		coefs.push_back((0.5*i*coefs[i-2] + x*beta*coefs[i-1]) / asum);

	// recursion in l2:
	// S_lA,lB = S_lA+1,lB-1 - x S_lA,lB-1
	for (int j = 0; j < lB; j++)
	{
		for (int i = powsum-1; i >= lA+j; i--)
			coefs[i] -= x * coefs[i-1];
	}

	return coefs[powsum-1];
}

Eigen::ArrayXXd gto_overlap_primitive_generic(
	const Eigen::Vector3i& lsA, const Eigen::Vector3i& lsB,
        const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
        const Eigen::ArrayXXd& asum,
        const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	Eigen::ArrayXXd Sp = exp_ared / (asum * asum.sqrt());

	for (int i = 0; i < 3; i++)
	{
		if (lsA[i]+lsB[i] > 0)
		{
			if (lsA[i] < lsB[i])
				Sp *= gto_overlap_primitive_1d(lsA[i], lsB[i],
					-r[i], alpha, asum);
			else
				Sp *= gto_overlap_primitive_1d(lsA[i], lsB[i],
					r[i], beta, asum);
		}
	}

	return Sp;
}

double gto_overlap_generic(const Eigen::Vector3i& lsA,
	const Eigen::VectorXd& weightsA, const Eigen::VectorXd& widthsA,
	const Eigen::Vector3d& posA,
	const Eigen::Vector3i& lsB,
	const Eigen::VectorXd& weightsB, const Eigen::VectorXd& widthsB,
	const Eigen::Vector3d& posB)
{
	int nA = widthsA.size(), nB = widthsB.size();
	Eigen::ArrayXXd alpha = widthsA.replicate(1, nB);
	Eigen::ArrayXXd beta = widthsB.transpose().replicate(nA, 1);
	Eigen::ArrayXXd asum = alpha + beta;
	Eigen::ArrayXXd ared = (widthsA * widthsB.transpose()).array() / asum;
	Eigen::Vector3d r = posB - posA;
	Eigen::ArrayXXd exp_ared = (-r.squaredNorm() * ared).exp();

	return Constants::pi_sqrt_pi * weightsA.transpose()
		* gto_overlap_primitive_generic(lsA, lsB, alpha, beta, asum,
			exp_ared, r).matrix()
		* weightsB;
}
