#ifndef OVERLAP_HH
#define OVERLAP_HH

#include <Eigen/Dense>
#include <Exception.hh>
#include "constants.hh"

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
struct OverlapKernel
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (r.x() * beta * OverlapKernel<lx1-1, ly1, lz1, lx2, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * (lx1-1) * OverlapKernel<lx1-2, ly1, lz1, lx2, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * lx2 * OverlapKernel<lx1-1, ly1, lz1, lx2-1, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			) / asum;
	}
};
template <>
struct OverlapKernel<0, 0, 0, 0, 0, 0>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d&)
	{
        	return exp_ared / (asum*asum.sqrt());
	}
};
template <>
struct OverlapKernel<1, 0, 0, 0, 0, 0>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return r.x() * beta * exp_ared / (asum.square()*asum.sqrt());
	}
};
template <>
struct OverlapKernel<0, 1, 0, 0, 0, 0>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return r.y() * beta * exp_ared / (asum.square()*asum.sqrt());
	}
};
template <>
struct OverlapKernel<0, 0, 1, 0, 0, 0>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return r.z() * beta * exp_ared / (asum.square()*asum.sqrt());
	}
};
template <>
struct OverlapKernel<0, 0, 0, 1, 0, 0>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return -r.x() * alpha * exp_ared / (asum.square()*asum.sqrt());
	}
};
template <>
struct OverlapKernel<0, 0, 0, 0, 1, 0>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return -r.y() * alpha * exp_ared / (asum.square()*asum.sqrt());
	}
};
template <>
struct OverlapKernel<0, 0, 0, 0, 0, 1>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return -r.z() * alpha * exp_ared / (asum.square()*asum.sqrt());
	}
};
template <>
struct OverlapKernel<0, 0, 1, 0, 0, 1>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (0.5 - r.z()*r.z()*ared)
			* exp_ared / (asum.square()*asum.sqrt());
	}
};
template <int lz1>
struct OverlapKernel<0, 0, lz1, 0, 0, 0>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (r.z() * beta * OverlapKernel<0, 0, lz1-1, 0, 0, 0>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * (lz1-1) * OverlapKernel<0, 0, lz1-2, 0, 0, 0>::calc(alpha, beta, asum, ared, exp_ared, r)
			) / asum;
	}
};
template <int lz2>
struct OverlapKernel<0, 0, 0, 0, 0, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return OverlapKernel<0, 0, lz2, 0, 0, 0>::calc(beta, alpha, asum, ared, exp_ared, -r);
	}
};
template <int lz2>
struct OverlapKernel<0, 0, 1, 0, 0, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return OverlapKernel<0, 0, lz2, 0, 0, 0>::calc(beta, alpha, asum, ared, exp_ared, -r);
	}
};
template <int lz1, int lz2>
struct OverlapKernel<0, 0, lz1, 0, 0, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (r.z() * beta * OverlapKernel<0, 0, lz1-1, 0, 0, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * (lz1-1) * OverlapKernel<0, 0, lz1-2, 0, 0, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * lz2 * OverlapKernel<0, 0, lz1-1, 0, 0, lz2-1>::calc(alpha, beta, asum, ared, exp_ared, r)
			) / asum;
	}
};
template <int lz1, int lz2>
struct OverlapKernel<0, 1, lz1, 0, 0, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return r.y() * beta
			* OverlapKernel<0, 0, lz1, 0, 0, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			/ asum;
	}
};
template <int lz1, int lz2>
struct OverlapKernel<0, 0, lz1, 0, 1, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return -r.y() * alpha
			* OverlapKernel<0, 0, lz1, 0, 0, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			/ asum;
	}
};
template <int lz1, int lz2>
struct OverlapKernel<0, 1, lz1, 0, 1, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (0.5 - r.y()*r.y()*ared)
			* OverlapKernel<0, 0, lz1, 0, 0, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			/ asum;
	}
};
template <int ly1, int lz1, int lz2>
struct OverlapKernel<0, ly1, lz1, 0, 0, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (r.y() * beta * OverlapKernel<0, ly1-1, lz1, 0, 0, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * (ly1-1) * OverlapKernel<0, ly1-2, lz1, 0, 0, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			) / asum;
	}
};
template <int lz1, int ly2, int lz2>
struct OverlapKernel<0, 0, lz1, 0, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return OverlapKernel<0, ly2, lz2, 0, 0, lz1>::calc(beta, alpha, asum, ared, exp_ared, -r);
	}
};
template <int lz1, int ly2, int lz2>
struct OverlapKernel<0, 1, lz1, 0, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return OverlapKernel<0, ly2, lz2, 0, 1, lz1>::calc(beta, alpha, asum, ared, exp_ared, -r);
	}
};
template <int ly1, int lz1, int ly2, int lz2>
struct OverlapKernel<0, ly1, lz1, 0, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (r.y() * beta * OverlapKernel<0, ly1-1, lz1, 0, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * (ly1-1) * OverlapKernel<0, ly1-2, lz1, 0, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * ly2 * OverlapKernel<0, ly1-1, lz1, 0, ly2-1, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			) / asum;
	}
};
template <int ly1, int lz1, int ly2, int lz2>
struct OverlapKernel<1, ly1, lz1, 0, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return r.x() * beta
			* OverlapKernel<0, ly1, lz1, 0, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			/ asum;
	}
};
template <int ly1, int lz1, int ly2, int lz2>
struct OverlapKernel<0, ly1, lz1, 1, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return -r.x() * alpha
			* OverlapKernel<0, ly1, lz1, 0, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			/ asum;
	}
};
template <int ly1, int lz1, int ly2, int lz2>
struct OverlapKernel<1, ly1, lz1, 1, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (0.5 - r.x()*r.x()*ared)
			* OverlapKernel<0, ly1, lz1, 0, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			/ asum;
	}
};
template <int lx1, int ly1, int lz1, int ly2, int lz2>
struct OverlapKernel<lx1, ly1, lz1, 0, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return (r.x() * beta * OverlapKernel<lx1-1, ly1, lz1, 0, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			+ 0.5 * (lx1-1) * OverlapKernel<lx1-2, ly1, lz1, 0, ly2, lz2>::calc(alpha, beta, asum, ared, exp_ared, r)
			) / asum;
	}
};
template <int ly1, int lz1, int lx2, int ly2, int lz2>
struct OverlapKernel<0, ly1, lz1, lx2, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return OverlapKernel<lx2, ly2, lz2, 0, ly1, lz1>::calc(beta, alpha, asum, ared, exp_ared, -r);
	}
};
template <int ly1, int lz1, int lx2, int ly2, int lz2>
struct OverlapKernel<1, ly1, lz1, lx2, ly2, lz2>
{
	static Eigen::ArrayXXd calc(
		const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
		const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
		const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
	{
		return OverlapKernel<lx2, ly2, lz2, 1, ly1, lz1>::calc(beta, alpha, asum, ared, exp_ared, -r);
	}
};

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
double gto_overlap(const Eigen::VectorXd& weights1,
	const Eigen::VectorXd& widths1, const Eigen::Vector3d& pos1,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& pos2)
{
	int n1 = widths1.size(), n2 = widths2.size();
	Eigen::ArrayXXd alpha = widths1.replicate(1, n2);
	Eigen::ArrayXXd beta = widths2.transpose().replicate(n1, 1);
	Eigen::ArrayXXd asum = alpha + beta;
        Eigen::ArrayXXd ared = (widths1 * widths2.transpose()).array() / asum;
	Eigen::Vector3d r = pos2 - pos1;
	Eigen::ArrayXXd exp_ared = (-r.squaredNorm()*ared).exp();

	return Constants::pi_sqrt_pi * weights1.transpose()
		* OverlapKernel<lx1, ly1, lz1, lx2, ly2, lz2>::calc(
			alpha, beta, asum, ared, exp_ared, r).matrix()
		* weights2;
}

#endif // OVERLAP_HH
