#include <Exception.hh>
#include "CGTOShellQuad.hh"
#include "boys.hh"

class Fms
{
	static const int pkt_size = Eigen::internal::packet_traits<double>::size;

public:
	typedef Eigen::ArrayXXd::AlignedMapType Block;
	typedef Eigen::ArrayXXd::ConstAlignedMapType ConstBlock;

	Fms(int maxm, const Eigen::ArrayXXd& T, const Eigen::ArrayXXd& expmT,
		const Eigen::ArrayXXd& KKW):
		_n1(T.rows()), _n2(T.cols()),
		_nn(pkt_size * ((_n1*_n2+pkt_size-1) / pkt_size)),
		_F(_nn, maxm+1)
	{
		(*this)(maxm) = T.boys(maxm, expmT);
		for (int m = maxm-1; m >= 0; --m)
			(*this)(m) = (expmT + 2*T*(*this)(m+1)) / (2*m+1);
		for (int m = maxm; m >= 0; --m)
			(*this)(m) *= KKW;
	}

	Block operator()(int m)
	{
		return Eigen::ArrayXXd::MapAligned(_F.data()+m*_nn, _n1, _n2);
	}
	ConstBlock operator()(int m) const
	{
		return Eigen::ArrayXXd::MapAligned(_F.data()+m*_nn, _n1, _n2);
	}

private:
	int _n1;
	int _n2;
	int _nn;
	Eigen::ArrayXXd _F;
};

CGTOShellQuad::CGTOShellQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD):
	_pAB(pAB), _pCD(pCD),
	_inv_widths_sum((widthsAB().replicate(1, pCD.size()).rowwise()
		+ widthsCD()).inverse()),
	_dPQ(_pAB.size(), 3*_pCD.size()),
	_lAB(pAB.lsum()), _lCD(pCD.lsum()), _lsum(_lAB+_lCD),
	_ints((_lAB+1)*(_lAB+1)*(_lAB+1), (_lCD+1)*(_lCD+1)*(_lCD+1)),
	_have_eri(false)
{
	if (_pAB.samePositionId())
	{
		if (_pCD.samePositionId())
		{
			_pos_sym = _pAB.positionIdA() == _pCD.positionIdA()
				? POS_SYM_AAAA : POS_SYM_AACC;
		}
		else
		{
			_pos_sym = POS_SYM_AACD;
		}
	}
	else if (_pCD.samePositionId())
	{
		_pos_sym = POS_SYM_ABCC;
	}
	else
	{
		_pos_sym = POS_SYM_ABCD;
	}

	_dPQ.block(0, 0, _pAB.size(), _pCD.size())
		= Q(0).replicate(_pAB.size(), 1).colwise() - P(0);
	_dPQ.block(0, _pCD.size(), _pAB.size(), _pCD.size())
		= Q(1).replicate(_pAB.size(), 1).colwise() - P(1);
	_dPQ.block(0, 2*_pCD.size(), _pAB.size(), _pCD.size())
		= Q(2).replicate(_pAB.size(), 1).colwise() - P(2);
	
	_T = widthsReduced() * (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	_expmT = (-_T).exp();
}

void CGTOShellQuad::elecRepPrim1d_abcd(int i, EriCoefs& coefs) const
{
	int l1 = _pAB.lsum(), l2 = _pCD.lsum();

	int idx = 0;
	// C_0,0,0,0
	coefs(idx++).setOnes();
	if (_lsum == 1)
	{
		if (l2 == 1)
		{
			// C_0,0,1,0
			coefs(idx++).rowwise() = dxQ(i, _pCD.centerA(i));
			coefs(idx) = dQW(i);
		}
		else
		{
			// C_1,0,0,0
			coefs(idx++).colwise() = dxP(i, _pAB.centerA(i));
			coefs(idx) = dPW(i);
		}
		return;
	}
	if (_lsum == 2)
	{
		if (l2 == 2)
		{
			auto dCQi = dxQ(i, _pCD.centerA(i));
			auto dQWi = dQW(i);

			// C_0,0,1,0
			coefs(idx++).rowwise() = dCQi;
			coefs(idx++) = dQWi;
			// C_0,0,2,0
			coefs(idx++).rowwise() = dCQi.square() + hInvWidthsCD();
			coefs(idx++) = dQWi.rowwise()*dCQi*2 - rho2();
			coefs(idx) = dQWi.square();
		}
		else if (l2 == 1)
		{
			auto dCQi = dxQ(i, _pCD.centerA(i));
			auto dQWi = dQW(i);
			auto dAPi = dxP(i, _pAB.centerA(i));
			auto dPWi = dPW(i);

			// C_0,0,1,0
			coefs(idx++).rowwise() = dCQi;
			coefs(idx++) = dQWi;
			// C_1,0,0,0
			coefs(idx++).colwise() = dAPi;
			coefs(idx++) = dPWi;
			// C_1,0,1,0
			coefs(idx++) = dCQi.replicate(_pAB.size(), 1).colwise()*dAPi;
			coefs(idx++) = dQWi.colwise()*dAPi + dPWi.rowwise()*dCQi
				+ 0.5*invWidthsSum();
			coefs(idx) = dPWi*dQWi;
		}
		else
		{
			auto dAPi = dxP(i, _pAB.centerA(i));
			auto dPWi = dPW(i);

			// C_1,0,0,0
			coefs(idx++).colwise() = dAPi;
			coefs(idx++) = dPWi;
			// C_2,0,0,0
			coefs(idx++).colwise() = dAPi.square() + hInvWidthsAB();
			coefs(idx++) = dPWi.colwise()*dAPi*2 - rho1();
			coefs(idx) = dPWi.square();
		}
		return;
	}

	if (l2 > 0)
	{
		auto dCQi = dxQ(i, _pCD.centerA(i));
		auto inv_eta = hInvWidthsCD();
		auto dQWi = dQW(i);
		auto rho2 = this->rho2();

		// C_0,0,1,0
		coefs(idx++).rowwise() = dCQi;
		coefs(idx++) = dQWi;

		// C_0,0,c,0
		for (int iC = 1; iC < l2; ++iC)
		{
			coefs(idx).rowwise() = coefs(idx-iC-1).row(0)*dCQi
				+ coefs(idx-2*iC-1).row(0)*(iC*inv_eta);
			++idx;
			for (int m = 1; m < iC; ++m, ++idx)
			{
				coefs(idx) = coefs(idx-iC-1).rowwise()*dCQi
					+ coefs(idx-2*iC-1).rowwise()*(iC*inv_eta)
					+ dQWi*coefs(idx-iC-2)
					- iC*rho2*coefs(idx-2*iC-2);
			}
			coefs(idx) = coefs(idx-iC-1).rowwise()*dCQi
				+ dQWi*coefs(idx-iC-2)
				- iC*rho2*coefs(idx-2*iC-2);
			++idx;
			coefs(idx) = dQWi*coefs(idx-iC-2);
			++idx;
		}
	}

	if (l1 > 0)
	{
		auto dAPi = dxP(i, _pAB.centerA(i));
		auto inv_zeta = hInvWidthsAB();
		Eigen::ArrayXXd dPWi = dPW(i);
		auto rho1 = this->rho1();
		auto inv_sum = 0.5 * invWidthsSum();

		int n2 = (l2+1)*(l2+2)/2;
		// C_1,0,c,0
		for (int iC = 0; iC <= l2; ++iC, ++n2)
		{
			coefs(idx) = coefs(idx-n2).colwise()*dAPi;
			++idx;
			for (int m = 1; m <= iC; ++m, ++idx)
			{
				coefs(idx) = coefs(idx-n2).colwise()*dAPi
					+ dPWi*coefs(idx-n2-1)
					+ iC*inv_sum*coefs(idx-n2-iC-1);
			}
			coefs(idx) = dPWi*coefs(idx-n2-1);
			++idx;
		}

		// C_a,0,c,0
		for (int iA = 1; iA < l1; ++iA)
		{
			coefs(idx) = coefs(idx-n2).colwise()*dAPi
				+ coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta);
			++idx;
			for (int m = 1; m < iA; ++m, ++idx)
			{
				coefs(idx) = coefs(idx-n2).colwise()*dAPi
					+ coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta)
					+ dPWi*coefs(idx-n2-1)
					- iA*rho1*coefs(idx-2*n2+l2);
			}
			coefs(idx) = coefs(idx-n2).colwise()*dAPi
				+ dPWi*coefs(idx-n2-1)
				- iA*rho1*coefs(idx-2*n2+l2);
			++idx;
			coefs(idx) = dPWi*coefs(idx-n2-1);
			++idx;
			++n2;

			for (int iC = 1; iC <= l2; ++iC, ++n2)
			{
				coefs(idx) = coefs(idx-n2).colwise()*dAPi
					+ coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta);
				++idx;
				for (int m = 1; m < iA+iC; ++m, ++idx)
				{
					coefs(idx) = coefs(idx-n2).colwise()*dAPi
						+ coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta)
						+ dPWi*coefs(idx-n2-1)
						- iA*rho1*coefs(idx-2*n2+l2)
						+ iC*inv_sum*coefs(idx-n2-iC-iA-1);
				}
				coefs(idx) = coefs(idx-n2).colwise()*dAPi
					+ dPWi*coefs(idx-n2-1)
					- iA*rho1*coefs(idx-2*n2+l2)
					+ iC*inv_sum*coefs(idx-n2-iC-iA-1);
				++idx;
				coefs(idx) = dPWi*coefs(idx-n2-1);
				++idx;
			}
		}
	}
}

double CGTOShellQuad::eri_xx(int lx1, int lx2, const EriCoefs& Cx, const Fms& fms,
	Eigen::ArrayXXd& Ctot) const
{
	int m = lx1 + lx2;
#ifdef DEBUG
	if (m == 0)
		throw Li::Exception("Total angular momentum must be greater than 0");
#endif

	EriCoefs::AllMBlock Cxl = Cx.allM(lx1, lx2);
	Ctot = Cxl(m) * fms(m);
	for (--m; m >= 0; --m)
		Ctot += Cxl(m) * fms(m);

	return mulWeights(Ctot);
}

double CGTOShellQuad::eri_xx(int lx1, int ly1, int lx2, int ly2,
	const EriCoefs& Cx, const EriCoefs& Cy,
	const Fms& fms, Eigen::ArrayXXd& Cm, Eigen::ArrayXXd& Ctot) const
{
	int lx = lx1 + lx2, ly = ly1 + ly2;

	if (lx == 0)
		return eri_xx(ly1, ly2, Cy, fms, Ctot);
	if (ly == 0)
		return eri_xx(lx1, lx2, Cx, fms, Ctot);

	int m = lx + ly;
	EriCoefs::AllMBlock Cxl = Cx.allM(lx1, lx2);
	EriCoefs::AllMBlock Cyl = Cy.allM(ly1, ly2);
	Ctot = Cxl(lx) * Cyl(ly) * fms(m);
	for (--m; m > 1; --m)
	{
		Cm.setZero();
		int mxmax = std::min(m, lx);
		int mxmin = std::max(m-ly, 0);
		for (int mx = mxmax; mx >= mxmin; --mx)
		{
			int my = m-mx;
			Cm += Cxl(mx) * Cyl(my);
		}
		Ctot += Cm * fms(m);
	}
	Ctot += (Cxl(1) * Cyl(0) + Cxl(0) * Cyl(1)) * fms(1);
	Ctot += Cxl(0) * Cyl(0) * fms(0);

	return mulWeights(Ctot);
}

double CGTOShellQuad::eri_xx(int lx1, int ly1, int lz1, int lx2, int ly2, int lz2,
	const EriCoefs& Cx, const EriCoefs& Cy, const EriCoefs& Cz,
	const Fms& fms, Eigen::ArrayXXd& Cm, Eigen::ArrayXXd& Ctot) const
{
	int lx = lx1 + lx2, ly = ly1 + ly2, lz = lz1 + lz2;
	if (lx == 0)
		return eri_xx(ly1, lz1, ly2, lz2, Cy, Cz, fms, Cm, Ctot);
	if (ly == 0)
		return eri_xx(lx1, lz1, lx2, lz2, Cx, Cz, fms, Cm, Ctot);
	if (lz == 0)
		return eri_xx(lx1, ly1, lx2, ly2, Cx, Cy, fms, Cm, Ctot);
	
	int m = lx + ly + lz;
	EriCoefs::AllMBlock Cxl = Cx.allM(lx1, lx2);
	EriCoefs::AllMBlock Cyl = Cy.allM(ly1, ly2);
	EriCoefs::AllMBlock Czl = Cz.allM(lz1, lz2);
	Ctot = Cxl(lx) * Cyl(ly) * Czl(lz) * fms(m);
	for (--m; m > 1; --m)
	{
		Cm.setZero();
		int mxmax = std::min(m, lx);
		for (int mx = mxmax; mx >= 0; --mx)
		{
			int mymax = std::min(m-mx, ly);
			int mymin = std::max(m-mx-lz, 0);
			for (int my = mymax; my >= mymin; --my)
			{
				int mz = m-mx-my;
				Cm += Cxl(mx) * Cyl(my) * Czl(mz);
			}
		}
		Ctot += Cm * fms(m);
	}
	Ctot += (Cxl(1) * Cyl(0) * Czl(0)
		+ Cxl(0) * Cyl(1) * Czl(0)
		+ Cxl(0) * Cyl(0) * Czl(1)) * fms(1);
	Ctot += Cxl(0) * Cyl(0) * Czl(0) * fms(0);

	return mulWeights(Ctot);
}

double CGTOShellQuad::eri_10(const EriCoefs::AllMBlock& Cxl,
	const Fms& fms) const
{
	return mulWeights(Cxl(1)*fms(1) + Cxl(0)*fms(0));
}

double CGTOShellQuad::eri_20(const EriCoefs::AllMBlock& Cxl,
	const Fms& fms) const
{
	return mulWeights(Cxl(2)*fms(2) + Cxl(1)*fms(1) + Cxl(0)*fms(0));
}

double CGTOShellQuad::eri_11(const EriCoefs::AllMBlock& Cxl,
	const EriCoefs::AllMBlock& Cyl,
	const Fms& fms, Eigen::ArrayXXd& Ctot) const
{
	Ctot = Cxl(1)*Cyl(1)*fms(2)
		+ (Cxl(1)*Cyl(0)+Cxl(0)*Cyl(1))*fms(1)
		+ Cxl(0)*Cyl(0)*fms(0);
	return mulWeights(Ctot);
}

void CGTOShellQuad::setEri() const
{
	Fms fms(_lsum, _T, _expmT, KKW());
	
	eri(0, 0, 0, 0, 0, 0) = mulWeights(fms(0));
	if (_lsum == 0)
	{
		_have_eri = true;
		return;
	}

	int n1 = _pAB.size(), n2 = _pCD.size();
	int lsum1 = _pAB.lsum(), lsum2 = _pCD.lsum();
	EriCoefs Cx(lsum1, lsum2, n1, n2);
	EriCoefs Cy(lsum1, lsum2, n1, n2);
	EriCoefs Cz(lsum1, lsum2, n1, n2);

	elecRepPrim1d_abcd(0, Cx);
	elecRepPrim1d_abcd(1, Cy);
	elecRepPrim1d_abcd(2, Cz);

	if (lsum1 >= 1)
	{
		eri(1, 0, 0, 0, 0, 0) = eri_10(Cx.allM(1, 0), fms);
		eri(0, 1, 0, 0, 0, 0) = eri_10(Cy.allM(1, 0), fms);
		eri(0, 0, 1, 0, 0, 0) = eri_10(Cz.allM(1, 0), fms);
		if (_lsum == 1)
		{
			_have_eri = true;
			return;
		}
	}
	if (lsum2 >= 1)
	{
		eri(0, 0, 0, 1, 0, 0) = eri_10(Cx.allM(0, 1), fms);
		eri(0, 0, 0, 0, 1, 0) = eri_10(Cy.allM(0, 1), fms);
		eri(0, 0, 0, 0, 0, 1) = eri_10(Cz.allM(0, 1), fms);
		if (_lsum == 1)
		{
			_have_eri = true;
			return;
		}
	}

	Eigen::ArrayXXd Ctot(n1, n2);
	if (lsum1 >= 2)
	{
		eri(2, 0, 0, 0, 0, 0) = eri_20(Cx.allM(2, 0), fms);
		eri(1, 1, 0, 0, 0, 0) = eri_11(Cx.allM(1, 0), Cy.allM(1, 0), fms, Ctot);
		eri(1, 0, 1, 0, 0, 0) = eri_11(Cx.allM(1, 0), Cz.allM(1, 0), fms, Ctot);
		eri(0, 2, 0, 0, 0, 0) = eri_20(Cy.allM(2, 0), fms);
		eri(0, 1, 1, 0, 0, 0) = eri_11(Cy.allM(1, 0), Cz.allM(1, 0), fms, Ctot);
		eri(0, 0, 2, 0, 0, 0) = eri_20(Cz.allM(2, 0), fms);
		if (_lsum == 2)
		{
			_have_eri = true;
			return;
		}
	}
	if (lsum1 >= 1 && lsum2 >= 1)
	{
		eri(1, 0, 0, 1, 0, 0) = eri_20(Cx.allM(1, 1), fms);
		eri(1, 0, 0, 0, 1, 0) = eri_11(Cx.allM(1, 0), Cy.allM(0, 1), fms, Ctot);
		eri(1, 0, 0, 0, 0, 1) = eri_11(Cx.allM(1, 0), Cz.allM(0, 1), fms, Ctot);
		eri(0, 1, 0, 1, 0, 0) = eri_11(Cy.allM(1, 0), Cx.allM(0, 1), fms, Ctot);
		eri(0, 1, 0, 0, 1, 0) = eri_20(Cy.allM(1, 1), fms);
		eri(0, 1, 0, 0, 0, 1) = eri_11(Cy.allM(1, 0), Cz.allM(0, 1), fms, Ctot);
		eri(0, 0, 1, 1, 0, 0) = eri_11(Cz.allM(1, 0), Cx.allM(0, 1), fms, Ctot);
		eri(0, 0, 1, 0, 1, 0) = eri_11(Cz.allM(1, 0), Cy.allM(0, 1), fms, Ctot);
		eri(0, 0, 1, 0, 0, 1) = eri_20(Cz.allM(1, 1), fms);
		if (_lsum == 2)
		{
			_have_eri = true;
			return;
		}
	}
	if (lsum2 >= 2)
	{
		eri(0, 0, 0, 2, 0, 0) = eri_20(Cx.allM(0, 2), fms);
		eri(0, 0, 0, 1, 1, 0) = eri_11(Cx.allM(0, 1), Cy.allM(0, 1), fms, Ctot);
		eri(0, 0, 0, 1, 0, 1) = eri_11(Cx.allM(0, 1), Cz.allM(0, 1), fms, Ctot);
		eri(0, 0, 0, 0, 2, 0) = eri_20(Cy.allM(0, 2), fms);
		eri(0, 0, 0, 0, 1, 1) = eri_11(Cy.allM(0, 1), Cz.allM(0, 1), fms, Ctot);
		eri(0, 0, 0, 0, 0, 2) = eri_20(Cz.allM(0, 2), fms);
		if (_lsum == 2)
		{
			_have_eri = true;
			return;
		}
	}

	Eigen::ArrayXXd Cm(n1, n2);
	for (int l = 3; l <= _lsum; ++l)
	{
		for (int l1 = lsum1; l1 >= std::max(l-lsum2, 0); --l1)
		{
			int l2 = l - l1;
			for (int lx1 = l1; lx1 >= 0; --lx1)
			{
				for (int lx2 = l2; lx2 >= 0; --lx2)
				{
					for (int ly1 = l1-lx1; ly1 >= 0; --ly1)
					{
						int lz1 = l1-lx1-ly1;
						for (int ly2 = l2-lx2; ly2 >= 0; --ly2)
						{
							int lz2 = l2-lx2-ly2;
							eri(lx1, ly1, lz1, lx2, ly2, lz2) = eri_xx(lx1, ly1, lz1, lx2, ly2, lz2, Cx, Cy, Cz, fms, Cm, Ctot);
						}
					}
				}
			}
		}
	}

	_have_eri = true;
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1,
	int lx2, int ly2, int lzC, int lzD) const
{
	if (lzD == 0)
		return eri(lx1, ly1, lz1, lx2, ly2, lzC);
	return eri(lx1, ly1, lz1, lx2, ly2, lzC+1, lzD-1)
		- _pCD.dAB(2) * eri(lx1, ly1, lz1, lx2, ly2, lzC, lzD-1);
	
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1,
	int lx2, int lyC, int lzC, int lyD, int lzD) const
{
	if (lyD == 0)
		return eri(lx1, ly1, lz1, lx2, lyC, lzC, lzD);
	return eri(lx1, ly1, lz1, lx2, lyC+1, lzC, lyD-1, lzD)
		- _pCD.dAB(1) * eri(lx1, ly1, lz1, lx2, lyC, lzC, lyD-1, lzD);
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1,
	int lxC, int lyC, int lzC, int lxD, int lyD, int lzD) const
{
	if (_pos_sym == POS_SYM_AACD || _pos_sym == POS_SYM_ABCD)
	{
		if (lxD == 0)
			return eri(lx1, ly1, lz1, lxC, lyC, lzC, lyD, lzD);
		return eri(lx1, ly1, lz1, lxC+1, lyC, lzC, lxD-1, lyD, lzD)
			- _pCD.dAB(0) * eri(lx1, ly1, lz1, lxC, lyC, lzC, lxD-1, lyD, lzD);
	}
	return eri(lx1, ly1, lz1, lxC+lxD, lyC+lyD, lzC+lzD);
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lzA, int lzB,
	int lxC, int lyC, int lzC, int lxD, int lyD, int lzD) const
{
	if (lzB == 0)
		return eri(lx1, ly1, lzA, lxC, lyC, lzC, lxD, lyD, lzD);
	return eri(lx1, ly1, lzA+1, lzB-1, lxC, lyC, lzC, lxD, lyD, lzD)
		- _pAB.dAB(2) * eri(lx1, ly1, lzA, lzB-1, lxC, lyC, lzC, lxD, lyD, lzD);
}

inline double CGTOShellQuad::eri(int lx1, int lyA, int lzA, int lyB, int lzB,
	int lxC, int lyC, int lzC, int lxD, int lyD, int lzD) const
{
	if (lyB == 0)
		return eri(lx1, lyA, lzA, lzB, lxC, lyC, lzC, lxD, lyD, lzD);
	return eri(lx1, lyA+1, lzA, lyB-1, lzB, lxC, lyC, lzC, lxD, lyD, lzD)
		- _pAB.dAB(1) * eri(lx1, lyA, lzA, lyB-1, lzB, lxC, lyC, lzC, lxD, lyD, lzD);
}

double CGTOShellQuad::eri(int lxA, int lyA, int lzA, int lxB, int lyB, int lzB,
	int lxC, int lyC, int lzC, int lxD, int lyD, int lzD) const
{
	if (_pos_sym == POS_SYM_ABCC || _pos_sym == POS_SYM_ABCD)
	{
		if (lxB == 0)
			return eri(lxA, lyA, lzA, lyB, lzB, lxC, lyC, lzC, lxD, lyD, lzD);
		if (lxB == 1)
			return eri(lxA+1, lyA, lzA, lyB, lzB, lxC, lyC, lzC, lxD, lyD, lzD)
				- _pAB.dAB(0) * eri(lxA, lyA, lzA, lyB, lzB, lxC, lyC, lzC, lxD, lyD, lzD);
		return eri(lxA+1, lyA, lzA, lxB-1, lyB, lzB, lxC, lyC, lzC, lxD, lyD, lzD)
			- _pAB.dAB(0) * eri(lxA, lyA, lzA, lxB-1, lyB, lzB, lxC, lyC, lzC, lxD, lyD, lzD);
	}
	return eri(lxA+lxB, lyA+lyB, lzA+lzB, lxC, lyC, lzC, lxD, lyD, lzD);
}
