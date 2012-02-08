#include "CGTOPair.hh"
#include "Dispatcher.hh"
#include "boys.hh"

const size_t CGTOPair::cid = Dispatcher::singleton().classID<CGTOPair>();

void CGTOPair::overlapPrim1D(int i, Eigen::ArrayXXd& Sp) const
{
	int lAB = lsum(i);
	if (lAB == 0)
		return;

	int lA = this->lA(i), lB = this->lB(i);
	double x = r(i);
	bool swapAB = lA < lB;

	if (swapAB)
	{
		std::swap(lA, lB);
		x = -x;
	}
	const Eigen::ArrayXXd& beta = swapAB ? widthsA() : widthsB();
	const Eigen::ArrayXXd& asum = widthsSum();

	if (lAB == 1)
	{
		Sp *= x*beta/asum;
		return;
	}

	std::vector<Eigen::ArrayXXd> coefs;
	coefs.push_back(x*beta/asum);
	// Recursion in lA:
	// S_lA+1,0 = [lA/2 S_lA-1,0 + x beta S_lA,0] / (alpha + beta)
	coefs.push_back((0.5 + x*beta*coefs[0]) / asum);
	for (int iA = 2; iA < lAB; i++)
	{
		coefs.push_back((0.5*iA*coefs[iA-2] + x*beta*coefs[iA-1])
			/ asum);
	}

	// recursion in lB:
	// S_lA,lB = S_lA+1,lB-1 - x S_lA,lB-1
	for (int iB = 0; iB < lB; ++iB)
	{
		for (int iA = lAB-1; iA >= lA+iB; --iA)
			coefs[iA] -= x*coefs[iA-1];
	}

	Sp *= coefs.back();
}

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

void CGTOPair::oneElecPrim1D(int i, Eigen::ArrayXXd& Sp,
	Eigen::ArrayXXd& Tp) const
{
	int lAB = lsum(i);
	if (lAB == 0)
		return;
	
	int lA = this->lA(i), lB = this->lB(i);
	double x = r(i);
	bool swapAB = lA < lB;
	
	if (swapAB)
	{
		std::swap(lA, lB);
		x = -x;
	}
	const Eigen::ArrayXXd& alpha = swapAB ? widthsB() : widthsA();
	const Eigen::ArrayXXd& beta = swapAB ? widthsA() : widthsB();
	const Eigen::ArrayXXd& asum = widthsSum();
	const Eigen::ArrayXXd& ared = widthsReduced();
	
	if (lAB == 1)
	{
		Sp *= x*beta/asum;
		Tp = x*beta*Tp/asum + 2*ared*Sp;
		return;
	}

	std::vector<Eigen::ArrayXXd> cS, cT;
	// Recursion in lA:
	// S_l+1,0 = [x beta S_l,0 + l/2 S_l-1,0] / asum
	// T_l+1,0 = [x beta T_l,0 + l/2 T_l-1,0 - l beta S_l-1,0] / asum + 2 ared Sl+1,0
	cS.push_back(Sp);
	cS.push_back(x*beta*Sp/asum);
	cT.push_back(Tp);
	cT.push_back(x*beta*Tp/asum + 2*ared*cS.back());

	for (int iA = 1; iA < lA; ++iA)
	{
		cS.push_back((x*beta*cS[iA] + 0.5*iA*cS[iA-1]) / asum);
		cT.push_back((x*beta*cT[iA] + 0.5*iA*cT[iA-1] - iA*beta*cS[iA-1]) / asum
			+ 2*ared*cS.back());
	}

	if (lB > 0)
	{
		// S_l1,l2+1 = [-x alpha S_l1,l2 + l2/2 S_l1,l2-1 + l1/2 S_l1-1,l2] / asum
		// T_l1,l2+1 = [-x alpha T_l1,l2 + l2/2 T_l1,l2-1 + l1/2 T_l1-1,l2] / asum
		//           + 2 ared S_l1,l2+1 - l2 alpha S_l1,l2-1 / asum
		std::vector<Eigen::ArrayXXd> cS_prev, cT_prev;

		cS_prev.swap(cS);
		cT_prev.swap(cT);

		cS.push_back(-x*alpha*cS_prev[0]/asum);
		cT.push_back(-x*alpha*cT_prev[0]/asum + 2*ared*cS[0]);
		for (int iA = 1; iA <= lA; ++iA)
		{
			cS.push_back((-x*alpha*cS_prev[iA] + 0.5*iA*cS_prev[iA-1]) / asum);
			cT.push_back((-x*alpha*cT_prev[iA] + 0.5*iA*cT_prev[iA-1]) / asum
				+ 2*ared*cS[iA]);
		}

		for (int iB = 1; iB < lB; ++iB)
		{
			cS_prev.swap(cS);
			cT_prev.swap(cT);

			cT[0] = (-x*alpha*cT_prev[0] + 0.5*iB*cT[0] - iB*alpha*cS[0]) / asum;
			cS[0] = (-x*alpha*cS_prev[0] + 0.5*iB*cS[0]) / asum;
			cT[0] += 2*ared*cS[0];
			for (int iA = 1; iA <= lA; ++iA)
			{
				cT[iA] = (-x*alpha*cT_prev[iA] + 0.5*iB*cT[iA]
					+ 0.5*iA*cT_prev[iA-1] - iB*alpha*cS[iA]) / asum;
				cS[iA] = (-x*alpha*cS_prev[iA] + 0.5*iB*cS[iA]
					+ 0.5*iA*cS_prev[iA-1]) / asum;
				cT[iA] += 2*ared*cS[iA];
			}
		}
	}

	Sp = cS.back();
	Tp = cT.back();
}

void CGTOPair::nucAttrPrim1D(int i, const Eigen::ArrayXXd& theta,
	const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& dPC,
	std::vector<Eigen::ArrayXXd>& res) const
{
	int rows = theta.rows();
	int cols = theta.cols();
	int lAB = lsum(i);
	if (lAB == 0)
	{
		res.assign(1, Eigen::ArrayXXd::Ones(rows, cols));
		return;
	}
	
	int lA = this->lA(i), lB = this->lB(i);
	double xA = centerA(i), xB = centerB(i);
	if (lA < lB)
	{
		std::swap(lA, lB);
		std::swap(xA, xB);
	}

	Eigen::ArrayXXd dPA = (Pi - xA).replicate(1, cols/Pi.cols());
	if (lAB == 1)
	{
		res.resize(2);
		res[0] = dPA;
		res[1] = -dPC;
		return;
	}

	std::vector< std::vector<Eigen::ArrayXXd> > As(lAB+1);
	// A_0,0
	As[0].push_back(Eigen::ArrayXXd::Ones(rows, cols));
	// A_1,0
	As[1].push_back(dPA);
	As[1].push_back(-dPC);
	// A_a,0
	for (int iA = 1; iA < lAB; ++iA)
	{
		As[iA+1].push_back(dPA*As[iA][0] + 0.5*iA*As[iA-1][0]/theta);
		for (int m = 1; m < iA; ++m)
		{
			As[iA+1].push_back(dPA*As[iA][m] - dPC*As[iA][m-1]
				+ 0.5*iA*(As[iA-1][m]-As[iA-1][m-1])/theta);
		}
		As[iA+1].push_back(dPA*As[iA][iA] - dPC*As[iA][iA-1]
			- 0.5*iA*As[iA-1][iA-1]/theta);
		As[iA+1].push_back(-dPC*As[iA][iA]);
	}

	double dAB = xB - xA;
	for (int iB = 1; iB <= lB; iB++)
	{
		for (int iA = lAB; iA >= lA+iB; iA--)
		{
			for (int m = 0; m < iA; m++)
				As[iA][m] -= dAB * As[iA-1][m];
		}
	}

	res.swap(As.back());
}

double CGTOPair::overlap() const
{
	Eigen::ArrayXXd Sp = gaussReduced() / (widthsSum() * widthsSum().sqrt());

	for (int i = 0; i < 3; i++)
		overlapPrim1D(i, Sp);

	return Constants::pi_sqrt_pi
		* f().weights().transpose() * Sp.matrix() * g().weights();
}

double CGTOPair::kineticEnergy() const
{
	Eigen::ArrayXXd Sp = Eigen::ArrayXXd::Zero(f().size(), g().size());
	Eigen::ArrayXXd Tp = widthsReduced()
		* (3 - 2*r().squaredNorm()*widthsReduced()) * Sp;

	for (int i = 0; i < 3; i++)
		oneElecPrim1D(i, Sp, Tp);

	return Constants::pi_sqrt_pi
		* f().weights().transpose() * Tp.matrix() * g().weights();
}

void CGTOPair::oneElectron(double &S, double &T) const
{
	Eigen::ArrayXXd Sp = gaussReduced() / (widthsSum() * widthsSum().sqrt());
	Eigen::ArrayXXd Tp = widthsReduced()
		* (3 - 2*r().squaredNorm()*widthsReduced()) * Sp;

	for (int i = 0; i < 3; i++)
		oneElecPrim1D(i, Sp, Tp);

	S = Constants::pi_sqrt_pi
		* f().weights().transpose() * Sp.matrix() * g().weights();
	T = Constants::pi_sqrt_pi
		* f().weights().transpose() * Tp.matrix() * g().weights();
}

double CGTOPair::nuclearAttraction(const Eigen::MatrixXd& nuc_pos,
	const Eigen::VectorXd& nuc_charge) const
{
	int nr_nuc = nuc_pos.cols();
	int nr_rows = f().size(), nr_cols = g().size();
	Eigen::ArrayXXd theta = widthsSum().replicate(1, nr_nuc);
	Eigen::ArrayXXd Pi;
	Eigen::ArrayXXd U = Eigen::ArrayXXd::Zero(nr_rows, nr_cols*nr_nuc);
	Eigen::ArrayXXd dPC(nr_rows, nr_cols*nr_nuc);
	std::vector<Eigen::ArrayXXd> Axyz[3];

	for (int i = 0; i < 3; i++)
	{
		Pi = P(i);
		for (int iC = 0; iC < nr_nuc; ++iC)
			dPC.block(0, iC*nr_cols, nr_rows, nr_cols) = Pi - nuc_pos(i, iC);
		U += dPC.square();

		nucAttrPrim1D(i, theta, Pi, dPC, Axyz[i]);
	}
	U *= theta;

	Eigen::Vector3i lsAB = f().ls() + g().ls();
	int lsABsum = lsAB.sum();

	Eigen::ArrayXXd F = Fm(lsABsum, U);
	Eigen::ArrayXXd expmU = (-U).exp();
	Eigen::ArrayXXd Am = Axyz[0][lsAB.x()] * Axyz[1][lsAB.y()]
		* Axyz[2][lsAB.z()] * F;
	Eigen::ArrayXXd A = Eigen::ArrayXXd::Zero(nr_rows, nr_cols);
	for (int iC = 0; iC < nr_nuc; ++iC)
		A += -nuc_charge[iC] * Am.block(0, iC*nr_cols, nr_rows, nr_cols);

	for (int m = lsABsum-1; m >= 0; --m)
	{
		Am.setZero();
		for (int ix = 0; ix <= std::min(m, lsAB.x()); ++ix)
		{
			for (int iy = 0; iy <= std::min(m-ix, lsAB.y()); ++iy)
			{
				int iz = m - ix - iy;
				if (iz <= lsAB.z())
				{
					Am += Axyz[0][ix] * Axyz[1][iy]
						* Axyz[2][iz];
				}
			}
		}
		F = (expmU + 2*U*F) / (2*m+1);
		Am *= F;

		for (int iC = 0; iC < nr_nuc; ++iC)
		{
			A += -nuc_charge[iC]
				* Am.block(0, iC*nr_cols, nr_rows, nr_cols);
		}
	}

	return 2 * M_PI * f().weights().transpose()
		* (gaussReduced() * A / widthsSum()).matrix()
		* g().weights();
}
