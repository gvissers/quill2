#include <vector>
#include <Exception.hh>
#include "gaussint/gto_one_elec.hh"
#include "constants.hh"

// N|l> = l|l-1>
// d/dx N|l> = l d/dx |l-1>
// N d/dx |l> = N (l|l-1> - 2 alpha|l+1>)
// = l(l-1)|l-2> - (l+1) 2 alpha |l>
// = l [(l-1)|l-2> - 2 alpha|l>] - 2 alpha|l>
// = l d/dx |l-1> - 2 alpha|l>
// = d/dx N|l> - 2 alpha|l>

// <la+1|lb> = (x_p-x_a)<la|lb> + 1/2asum <N la|lb> + 1/2asum <la|N lb>
// <la+1|d^2/dx^2|lb>
// = (x_p-x_a)<la|d^2/dx^2|lb> + 1/2asum <N la|d^2/dx^2|lb>
//   + 1/2asum<la|N d^2/dx^2|lb>
// = (x_p-x_a)<la|d^2/dx^2|lb> + 1/2asum <N la|d^2/dx^2|lb>
//   + 1/2asum <la|d/dx N d/dx|lb> - beta/asum<la|d/dx|lb>
// = (x_p-x_a)<la|d^2/dx^2|lb> + 1/2asum <N la|d^2/dx^2|lb>
//   + 1/2asum <la|d^2/dx^2 N|lb> - 2beta/asum<la|d/dx|lb>
// = (x_p-x_a)<la|d^2/dx^2|lb> + 1/2asum <N la|d^2/dx^2|lb>
//   + 1/2asum <la|d^2/dx^2 N|lb> + 2beta/asum<d/dx la|lb>
// = (x_p-x_a)<la|d^2/dx^2|lb> + 1/2asum <N la|d^2/dx^2|lb>
//   + 1/2asum <la|d^2/dx^2 N|lb> + 2la beta/asum <la-1|lb>
//   - 4 xi <la+1|lb>
// = (x_p-x_a)/asum<la|d^2/dx^2|lb> + la/2asum <la-1|d^2/dx^2|lb>
//   + lb/2asum <la|d^2/dx^2|lb-1>
//   - 4 xi{<la+1|lb> - la/2alpha <la-1|lb>}
// = {beta(x_b - x_a) <la|d^2/dx^2|lb> + la/2 <la-1|d^2/dx^2|lb>
//   + lb/2 <la|d^2/dx^2|lb-1>} / asum
//   - 4 xi<la+1|lb> + 2la beta/asum <la-1|lb>

// T_la+1,lb = [x beta T_la,lb + la/2 T_la-1,lb + lb/2 T_la,lb-1] / asum
//   + 2 ared S_la+1,lb - la beta/asum S_la-1,lb

void gto_one_elec_primitive_generic_1d(int l1, int l2,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	double x, Eigen::ArrayXXd& Sp, Eigen::ArrayXXd& Tp)
{
#ifdef DEBUG
	if (l1 < l2)
		throw Li::Exception("l1 must not be smaller than l2");
#endif

	int powsum = l1 + l2;
	if (powsum == 0)
		return;
	if (powsum == 1)
	{
		Sp *= x*beta/asum;
		Tp = x*beta*Tp/asum + 2*ared*Sp;
		return;
	}

	std::vector<Eigen::ArrayXXd> cS, cT;
	// Recursion in l1:
	// S_l+1,0 = [x beta S_l,0 + l/2 S_l-1,0] / asum
	// T_l+1,0 = [x beta T_l,0 + l/2 T_l-1,0 - l beta S_l-1,0] / asum + 2 ared Sl+1,0

	cS.push_back(Sp);
	cS.push_back(x*beta*Sp/asum);
	cT.push_back(Tp);
	cT.push_back(x*beta*Tp/asum + 2*ared*cS.back());
	for (int i1 = 1; i1 < l1; i1++)
	{
		cS.push_back((x*beta*cS[i1] + 0.5*i1*cS[i1-1]) / asum);
		cT.push_back((x*beta*cT[i1] + 0.5*i1*cT[i1-1] - i1*beta*cS[i1-1]) / asum + 2*ared*cS.back());
	}

	if (l2 > 0)
	{
		// S_l1,l2+1 = [-x alpha S_l1,l2 + l2/2 S_l1,l2-1 + l1/2 S_l1-1,l2] / asum
		// T_l1,l2+1 = [-x alpha T_l1,l2 + l2/2 T_l1,l2-1 + l1/2 T_l1-1,l2] / asum
		//           + 2 ared S_l1,l2+1 - l2 alpha S_l1,l2-1 / asum
		std::vector<Eigen::ArrayXXd> cS_prev, cT_prev;

		cS_prev.swap(cS);
		cT_prev.swap(cT);

		cS.push_back(-x*alpha*cS_prev[0]/asum);
		cT.push_back(-x*alpha*cT_prev[0]/asum + 2*ared*cS[0]);
		for (int i1 = 1; i1 <= l1; i1++)
		{
			cS.push_back((-x*alpha*cS_prev[i1] + 0.5*i1*cS_prev[i1-1]) / asum);
			cT.push_back((-x*alpha*cT_prev[i1] + 0.5*i1*cT_prev[i1-1]) / asum + 2*ared*cS[i1]);
		}

		for (int i2 = 1; i2 < l2; i2++)
		{
			cS_prev.swap(cS);
			cT_prev.swap(cT);

			cT[0] = (-x*alpha*cT_prev[0] + 0.5*i2*cT[0] - i2*alpha*cS[0]) / asum;
			cS[0] = (-x*alpha*cS_prev[0] + 0.5*i2*cS[0]) / asum;
			cT[0] += 2*ared*cS[0];
			for (int i1 = 1; i1 <= l1; i1++)
			{
				cT[i1] = (-x*alpha*cT_prev[i1] + 0.5*i2*cT[i1] + 0.5*i1*cT_prev[i1-1] - i2*alpha*cS[i1]) / asum;
				cS[i1] = (-x*alpha*cS_prev[i1] + 0.5*i2*cS[i1] + 0.5*i1*cS_prev[i1-1]) / asum;
				cT[i1] += 2*ared*cS[i1];
			}
		}
	}

	Sp = cS.back();
	Tp = cT.back();
}

void gto_one_elec_primitive_generic(
	const Eigen::Vector3i& ls1, const Eigen::Vector3i& ls2,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r,
	Eigen::ArrayXXd& Sp, Eigen::ArrayXXd& Tp)
{
	Sp = exp_ared / (asum * asum.sqrt());
	Tp = ared * (3 - 2*r.squaredNorm()*ared) * Sp;

	for (int i = 0; i < 3; i++)
	{
		if (ls1[i] + ls2[i] == 0)
			continue;

		if (ls1[i] < ls2[i])
		{
			gto_one_elec_primitive_generic_1d(ls2[i], ls1[i],
				beta, alpha, asum, ared, -r[i], Sp, Tp);
		}
		else
		{
			gto_one_elec_primitive_generic_1d(ls1[i], ls2[i],
				alpha, beta, asum, ared, r[i], Sp, Tp);
		}
	}
}

void gto_one_elec_generic(const Eigen::Vector3i& lsA,
	const Eigen::VectorXd& weightsA, const Eigen::VectorXd& widthsA,
	const Eigen::Vector3d& posA,
	const Eigen::Vector3i& lsB,
	const Eigen::VectorXd& weightsB, const Eigen::VectorXd& widthsB,
	const Eigen::Vector3d& posB,
	double *S, double *T)
{
	int nA = widthsA.size(), nB = widthsB.size();
	Eigen::ArrayXXd alpha = widthsA.replicate(1, nB);
	Eigen::ArrayXXd beta = widthsB.transpose().replicate(nA, 1);
	Eigen::ArrayXXd asum = alpha + beta;
	Eigen::ArrayXXd ared = (widthsA * widthsB.transpose()).array() / asum;
	Eigen::Vector3d r = posB - posA;
	Eigen::ArrayXXd exp_ared = (-r.squaredNorm() * ared).exp();
	Eigen::ArrayXXd Sp, Tp;

	gto_one_elec_primitive_generic(lsA, lsB, alpha, beta, asum, ared,
		exp_ared, r, Sp, Tp);

	*S = Constants::pi_sqrt_pi
		* weightsA.transpose() * Sp.matrix() * weightsB;
	*T = Constants::pi_sqrt_pi
		* weightsA.transpose() * Tp.matrix() * weightsB;
}

