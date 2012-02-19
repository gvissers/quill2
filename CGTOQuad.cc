#include "CGTOQuad.hh"
#include "BFQuadPool.hh"
#include "boys.hh"

class EriCoefs
{
public:
	EriCoefs(int l1, int l2, int p_size, int q_size):
		_imax(l1+1), _jmax(l2+1), _p_size(p_size), _q_size(q_size),
		_C(p_size, (_imax*_jmax*(_imax+_jmax)*q_size)/2) {}

	int index(int i, int j, int m) const
	{
		return (i*_jmax*(i+_jmax) + j*(2*i+j+1)) / 2 + m;
	}
	
	Eigen::Block<Eigen::ArrayXXd> operator()(int i, int j, int m)
	{
		return (*this)(index(i,j,m));
	}
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int i, int j,
		int m) const
	{
		return (*this)(index(i,j,m));
	}
	Eigen::Block<Eigen::ArrayXXd> operator()(int idx)
	{
		return _C.block(0, idx*_q_size, _p_size, _q_size);
	}
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int idx) const
	{
		return _C.block(0, idx*_q_size, _p_size, _q_size);
	}

	Eigen::Block<Eigen::ArrayXXd> operator()(int i, int j, int m, int count)
	{
		return _C.block(0, (index(i,j,m))*_q_size, _p_size,
			count*_q_size);
	}
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int i, int j,
		int m, int count) const
	{
		return _C.block(0, (index(i,j,m))*_q_size, _p_size,
			count*_q_size);
	}

private:
	int _imax;
	int _jmax;
	int _p_size;
	int _q_size;
	Eigen::ArrayXXd _C;
};

class FmCoefs
{
public:
	FmCoefs(int mmax, int p_size, int q_size):
		_p_size(p_size), _q_size(q_size), _C(p_size, (mmax+1)*q_size),
		_mmin(0), _mmax(0)
	{
		block(0).setOnes();
	}

	int maxM() const { return _mmax; }
	int minM() const { return _mmin; }

	Eigen::Block<Eigen::ArrayXXd> operator[](int m)
	{
		return block(m);
	}
	const Eigen::Block<const Eigen::ArrayXXd> operator[](int m) const
	{
		return block(m);
	}

	Eigen::Block<Eigen::ArrayXXd> block(int m, int count=1)
	{
		return _C.block(0, (m-_mmin)*_q_size, _p_size, count*_q_size);
	}
	const Eigen::Block<const Eigen::ArrayXXd> block(int m, int count=1) const
	{
		return _C.block(0, (m-_mmin)*_q_size, _p_size, count*_q_size);
	}

	template <typename Derived0, typename Derived1>
	void multiplyCol(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1)
	{
		block(_mmax+1) = block(_mmax) * C1;
		for (int m = _mmax; m > _mmin; --m)
			block(m) = block(m).colwise()*C0 + block(m-1)*C1;
		block(_mmin).colwise() *= C0;
		++_mmax;
	}
	template <typename Derived0, typename Derived1>
	void multiplyRow(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1)
	{
		block(_mmax+1) = block(_mmax) * C1;
		for (int m = _mmax; m > _mmin; --m)
			block(m) = block(m).rowwise()*C0 + block(m-1)*C1;
		block(_mmin).rowwise() *= C0;
		++_mmax;
	}
	template <typename Derived0, typename Derived1>
	void multiply(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1)
	{
		block(_mmax+1) = block(_mmax) * C1;
		for (int m = _mmax; m > _mmin; --m)
			block(m) = block(m)*C0 + block(m-1)*C1;
		block(_mmin) *= C0;
		++_mmax;
	}
	template <typename Derived1>
	void multiply_noC0(const Eigen::ArrayBase<Derived1>& C1)
	{
		block(_mmin, _mmax-_mmin+1) *= C1.replicate(1, _mmax-_mmin+1);
		++_mmin;
		++_mmax;
	}
	template <typename Derived0, typename Derived1, typename Derived2>
	void multiplyCol(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2)
	{
		if (_mmax == _mmin)
		{
			block(_mmin+2) = block(_mmin)*C2;
			block(_mmin+1) = block(_mmin)*C1;
			block(_mmin).colwise() *= C0;
		}
		else
		{
			block(_mmax+2) = block(_mmax)*C2;
			block(_mmax+1) = block(_mmax)*C1 + block(_mmax-1)*C2;
			for (int m = _mmax; m > _mmin+1; --m)
			{
				block(m) = block(m).colwise()*C0 + block(m-1)*C1
					+ block(m-2)*C2;
			}
			block(_mmin+1) = block(_mmin+1).colwise()*C0 + block(_mmin)*C1;
			block(_mmin).colwise() *= C0;
		}
		_mmax += 2;
	}
	template <typename Derived0, typename Derived1, typename Derived2>
	void multiplyRow(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2)
	{
		if (_mmax == _mmin)
		{
			block(_mmin+2) = block(_mmin)*C2;
			block(_mmin+1) = block(_mmin)*C1;
			block(_mmin).rowwise() *= C0;
		}
		else
		{
			block(_mmax+2) = block(_mmax)*C2;
			block(_mmax+1) = block(_mmax)*C1 + block(_mmax-1)*C2;
			for (int m = _mmax; m > _mmin+1; --m)
			{
				block(m) = block(m).rowwise()*C0 + block(m-1)*C1
					+ block(m-2)*C2;
			}
			block(_mmin+1) = block(_mmin+1).rowwise()*C0 + block(_mmin)*C1;
			block(_mmin).rowwise() *= C0;
		}
		_mmax += 2;
	}
	template <typename Derived0, typename Derived1, typename Derived2>
	void multiply(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2)
	{
		if (_mmax == _mmin)
		{
			block(_mmin+2) = block(_mmin)*C2;
			block(_mmin+1) = block(_mmin)*C1;
			block(_mmin) *= C0;
		}
		else
		{
			block(_mmax+2) = block(_mmax)*C2;
			block(_mmax+1) = block(_mmax)*C1 + block(_mmax-1)*C2;
			for (int m = _mmax; m > _mmin+1; --m)
			{
				block(m) = block(m)*C0 + block(m-1)*C1
					+ block(m-2)*C2;
			}
			block(_mmin+1) = block(_mmin+1)*C0 + block(_mmin)*C1;
			block(_mmin) *= C0;
		}
		_mmax += 2;
	}
	template <typename Derived1, typename Derived2>
	void multiply_noC0(const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2)
	{
		if (_mmax == _mmin)
		{
			block(_mmin+1) = block(_mmin)*C2;
			block(_mmin) *= C1;
		}
		else
		{
			block(_mmax+1) = block(_mmax)*C2;
			for (int m = _mmax; m > _mmin; --m)
				block(m) = block(m)*C1 + block(m-1)*C2;
			block(_mmin) *= C1;
		}
		++_mmin;
		_mmax += 2;
	}
	void multiply(const EriCoefs& coefs, int l1, int l2)
	{
		int lsum = l1+l2;
		for (int m = lsum+_mmax; m >= _mmin; --m)
		{
			// i >= 0
			// i <= lsum
			// m-i >= _mmin => i <= m-_mmin
			// m-i <= _mmax => i >= m-_mmax
			int i = std::max(0, m-_mmax);
			block(m) = coefs(l1, l2, i) * block(m-i);
			for (++i; i <= std::min(lsum, m-_mmin); ++i)
				block(m) += coefs(l1, l2, i) * block(m-i);
		}
		_mmax += lsum;
	}

private:
	int _p_size;
	int _q_size;
	Eigen::ArrayXXd _C;
	int _mmin;
	int _mmax;
};

static const CGTOPair& firstPair(const CGTOPair& pp, const CGTOPair& qq)
{
	return pp.ishellPair() < qq.ishellPair() ? qq : pp;
}

static const CGTOPair& secondPair(const CGTOPair& pp, const CGTOPair& qq)
{
	return pp.ishellPair() < qq.ishellPair() ? pp : qq;
}

CGTOQuad::CGTOQuad(const CGTOPair& pp, const CGTOPair& qq):
	AbstractBFQuad(firstPair(pp, qq), secondPair(pp, qq)),
	_ishell_quad(CGTOShellList::pairIndex(p().ishellPair(), q().ishellPair())),
	_shell_quad(CGTOShellList::singleton().quad(_ishell_quad)) {}

void CGTOQuad::elecRepPrim1d_aacc_psss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 1)
		throw Li::Exception("Not a (p,s,s,s) quad");
#endif

	double dPQi = q().centerA(i) - p().centerA(i);
	if (lAB(i) == 1)
		Cm.multiply_noC0(invWidthsSum().rowwise() * (dPQi * widthsCD()));
	else
		Cm.multiply_noC0(invWidthsSum().colwise() * (-dPQi * widthsAB()));
}

void CGTOQuad::elecRepPrim1d_abcc_psss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 1)
		throw Li::Exception("Not a (p,s,s,s) quad");
#endif

	if (lAB(i) == 1)
	{
		double x = lA(i) == 1 ? centerA(i) : centerB(i);
		Cm.multiplyCol(dxP(i, x), dPW(i));
	}
	else
	{
		Cm.multiply_noC0(dQW(i));
	}
}

void CGTOQuad::elecRepPrim1d_aacd_psss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 1)
		throw Li::Exception("Not a (p,s,s,s) quad");
#endif

	if (lAB(i) == 1)
	{
		Cm.multiply_noC0(dPW(i));
	}
	else
	{
		double x = lC(i) == 1 ? centerC(i) : centerD(i);
		Cm.multiplyRow(dxQ(i, x), dQW(i));
	}
}

void CGTOQuad::elecRepPrim1d_abcd_psss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 1)
		throw Li::Exception("Not a (p,s,s,s) quad");
#endif

	if (lAB(i) == 1)
	{
		double x = lA(i) == 1 ? centerA(i) : centerB(i);
		Cm.multiplyCol(dxP(i, x), dPW(i));
	}
	else
	{
		double x = lC(i) == 1 ? centerC(i) : centerD(i);
		Cm.multiplyRow(dxQ(i, x), dQW(i));
	}
}

void CGTOQuad::elecRepPrim1d_abcc_ppss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| ((lA(i) != 1 || lB(i) != 1) && (lC(i) != 1 || lD(i) != 1)))
		throw Li::Exception("Not a (p,p,s,s) quad");
#endif

	if (lAB(i) > 0)
	{
		double xA = centerA(i), xB = centerB(i), dAB = xB - xA;
		auto dAPi = dxP(i, xA);
		auto dPWi = dPW(i);

		Cm.multiplyCol(dAPi.square() + hInvWidthsAB() - dAB*dAPi,
			dPWi.colwise()*(2*dAPi - dAB) - rho1(),
			dPWi.square());
	}
	else
	{
		Cm.multiplyRow(hInvWidthsCD(), -rho2(), dQW(i).square());
	}
}

void CGTOQuad::elecRepPrim1d_aacd_ppss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| ((lA(i) != 1 || lB(i) != 1) && (lC(i) != 1 || lD(i) != 1)))
		throw Li::Exception("Not a (p,p,s,s) quad");
#endif

	if (lAB(i) > 0)
	{
		Cm.multiplyCol(hInvWidthsAB(), -rho1(), dPW(i).square());
	}
	else
	{
		double xC = centerC(i), xD = centerD(i), dCD = xD - xC;
		auto dCQi = dxQ(i, xC);
		auto dQWi = dQW(i);

		Cm.multiplyRow(dCQi.square() + hInvWidthsCD() - dCD*dCQi,
			dQWi.rowwise()*(2*dCQi - dCD) - rho2(),
			dQWi.square());
	}
}

void CGTOQuad::elecRepPrim1d_abcd_ppss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| ((lA(i) != 1 || lB(i) != 1) && (lC(i) != 1 || lD(i) != 1)))
		throw Li::Exception("Not a (p,p,s,s) quad");
#endif

	if (lAB(i) > 0)
	{
		double xA = centerA(i), xB = centerB(i), dAB = xB - xA;
		auto dAPi = dxP(i, xA);
		auto dPWi = dPW(i);

		Cm.multiplyCol(dAPi.square() + hInvWidthsAB() - dAB*dAPi,
			dPWi.colwise()*(2*dAPi - dAB) - rho1(),
			dPWi.square());
	}
	else
	{
		double xC = centerC(i), xD = centerD(i), dCD = xD - xC;
		auto dCQi = dxQ(i, xC);
		auto dQWi = dQW(i);

		Cm.multiplyRow(dCQi.square() + hInvWidthsCD() - dCD*dCQi,
			dQWi.rowwise()*(2*dCQi - dCD) - rho2(),
			dQWi.square());
	}
}

void CGTOQuad::elecRepPrim1d_abcc_dsss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| (lA(i) != 2 && lB(i) != 2 && lC(i) != 2 && lD(i) != 2))
		throw Li::Exception("Not a (d,s,s,s) quad");
#endif

	if (lAB(i) > 0)
	{
		double x = lA(i) > 0 ? centerA(i) : centerB(i);
		auto dAPi = dxP(i, x);
		auto dPWi = dPW(i);

		Cm.multiplyCol(dAPi.square() + hInvWidthsAB(),
			dPWi.colwise()*(2*dAPi) - rho1(),
			dPWi.square());
	}
	else
	{
		Cm.multiplyRow(hInvWidthsCD(), -rho2(), dQW(i).square());
	}
}

void CGTOQuad::elecRepPrim1d_aacd_dsss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| (lA(i) != 2 && lB(i) != 2 && lC(i) != 2 && lD(i) != 2))
		throw Li::Exception("Not a (d,s,s,s) quad");
#endif

	if (lAB(i) > 0)
	{
		Cm.multiplyCol(hInvWidthsAB(), -rho1(), dPW(i).square());
	}
	else
	{
		double x = lC(i) > 0 ? centerC(i) : centerD(i);
		auto dCQi = dxQ(i, x);
		auto dQWi = dQW(i);

		Cm.multiplyRow(dCQi.square() + hInvWidthsCD(),
			dQWi.rowwise()*(2*dCQi) - rho2(),
			dQWi.square());
	}
}

void CGTOQuad::elecRepPrim1d_abcd_dsss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| (lA(i) != 2 && lB(i) != 2 && lC(i) != 2 && lD(i) != 2))
		throw Li::Exception("Not a (d,s,s,s) quad");
#endif

	if (lAB(i) > 0)
	{
		double x = lA(i) > 0 ? centerA(i) : centerB(i);
		auto dAPi = dxP(i, x);
		auto dPWi = dPW(i);

		Cm.multiplyCol(dAPi.square() + hInvWidthsAB(),
			dPWi.colwise()*(2*dAPi) - rho1(),
			dPWi.square());
	}
	else
	{
		double x = lC(i) > 0 ? centerC(i) : centerD(i);
		auto dCQi = dxQ(i, x);
		auto dQWi = dQW(i);

		Cm.multiplyRow(dCQi.square() + hInvWidthsCD(),
			dQWi.rowwise()*(2*dCQi) - rho2(),
			dQWi.square());
	}
}

void CGTOQuad::elecRepPrim1d_abcc_psps(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	int l1 = lA(i)+lB(i), l2 = lC(i)+lD(i);
	if (l1 != 1 || l2 != 1)
		throw Li::Exception("Not a (p,s,p,s) quad");
#endif

	double xA = lA(i) == 1 ? centerA(i) : centerB(i);
	auto dQWi = dQW(i);

	Cm.multiply_noC0(dQWi.colwise()*dxP(i, xA) + 0.5*invWidthsSum(),
		dPW(i)*dQWi);
}

void CGTOQuad::elecRepPrim1d_aacd_psps(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	int l1 = lA(i)+lB(i), l2 = lC(i)+lD(i);
	if (l1 != 1 || l2 != 1)
		throw Li::Exception("Not a (p,s,p,s) quad");
#endif

	double xC = lC(i) == 1 ? centerC(i) : centerD(i);
	auto dPWi = dPW(i);

	Cm.multiply_noC0(dPWi.rowwise()*dxQ(i, xC) + 0.5*invWidthsSum(),
		dPWi*dQW(i));
}

void CGTOQuad::elecRepPrim1d_abcd_psps(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	int l1 = lA(i)+lB(i), l2 = lC(i)+lD(i);
	if (l1 != 1 || l2 != 1)
		throw Li::Exception("Not a (p,s,p,s) quad");
#endif

	double xA = lA(i) == 1 ? centerA(i) : centerB(i);
	double xC = lC(i) == 1 ? centerC(i) : centerD(i);
	auto dAPi = dxP(i, xA);
	auto dCQi = dxQ(i, xC);
	auto dPWi = dPW(i);
	auto dQWi = dQW(i);

	Cm.multiply(dAPi.replicate(1, q().size()).rowwise()*dCQi,
		dQWi.colwise()*dAPi + dPWi.rowwise()*dCQi + 0.5*invWidthsSum(),
		dPWi*dQWi);
}

void CGTOQuad::elecRepPrim1d_aaaa(int i, FmCoefs& Cm) const
{
	int lA = this->lA(i), lB = this->lB(i),
		lC = this->lC(i), lD = this->lD(i);
	int l1 = lA+lB, l2 = lC+lD, lsum = l1+l2;

#ifdef DEBUG
	if (lsum % 2 == 1)
		throw Li::Exception("Odd total angular momentum not allowed");
#endif

	if (lsum == 0)
	{
		return;
	}
	if (lsum == 2)
	{
		if (l1 == 0)
			Cm.multiplyRow(hInvWidthsCD(), -rho2());
		else if (l1 == 1)
			Cm.multiply_noC0(0.5*invWidthsSum());
		else
 			Cm.multiplyCol(hInvWidthsAB(), -rho1());
		return;
	}
	if (lsum == 4)
	{
		if (l1 == 0)
		{
			auto inv_eta = hInvWidthsCD();
			auto rho2 = this->rho2();
			Cm.multiplyRow(3*inv_eta.square(),
				rho2.rowwise()*(-6*inv_eta), 3*rho2.square());
		}
		else if (l1 == 1)
		{
			Cm.multiply_noC0(invWidthsSum().rowwise()*(1.5*hInvWidthsCD()),
				-1.5*invWidthsSum()*rho2());
		}		
		else if (l1 == 2)
		{
			auto inv_zeta = hInvWidthsAB();
			auto inv_eta = hInvWidthsCD();
			auto rho1 = this->rho1();
			auto rho2 = this->rho2();
			Cm.multiply(inv_zeta.replicate(1, q().size()).rowwise()*inv_eta,
				rho2.colwise()*(-inv_zeta) - rho1.rowwise()*inv_eta,
				rho1*rho2 + 0.5*invWidthsSum().square());
		}
		else if (l1 == 3)
		{
			Cm.multiply_noC0(invWidthsSum().colwise()*(1.5*hInvWidthsAB()),
				-1.5*invWidthsSum()*rho1());
		}
		else
		{
			auto inv_zeta = hInvWidthsAB();
			auto rho1 = this->rho1();
			Cm.multiplyCol(3*inv_zeta.square(),
				rho1.colwise()*(-6*inv_zeta), 3*rho1.square());
		}
		return;
	}

	EriCoefs coefs(l1, l2, p().size(), q().size());
	// C_0,0,0,0
	coefs(0).setOnes();

	// C_0,0,c,0
	if (l2 > 1)
	{
		auto inv_eta = hInvWidthsCD();
		auto rho2 = this->rho2();
		
		coefs(0, 2, 0).rowwise() = hInvWidthsCD();
		coefs(0, 2, 1) = -rho2;
		coefs(0, 2, 2).setZero();

		for (int iC = 3; iC < l2; iC += 2)
		{
			coefs(0, iC+1, 0) = coefs(0, iC-1, 0).rowwise()*(iC*inv_eta);
			for (int m = 1; m < iC-1; ++m)
			{
				coefs(0, iC+1, m) = coefs(0, iC-1, m).rowwise()*(iC*inv_eta)
					- iC*rho2*coefs(0, iC-1, m-1);
			}
			coefs(0, iC+1, iC-1) = -iC*rho2*coefs(0, iC-1, iC-2);
			coefs(0, iC+1, iC).setZero();
			coefs(0, iC+1, iC+1).setZero();
		}
	}

	if (l1 > 0)
	{
		auto inv_sum = 0.5 * invWidthsSum();

		// C_1,0,c,0
		if (l2 > 0)
		{
			coefs(1, 1, 0).setZero();
			coefs(1, 1, 1) = inv_sum;
			coefs(1, 1, 2).setZero();
			for (int iC = 3; iC <= l2; iC += 2)
			{
				coefs(1, iC, 0).setZero();
				for (int m = 1; m < iC; ++m)
					coefs(1, iC, m) = iC*inv_sum*coefs(0, iC-1, m-1);
				coefs(1, iC, iC).setZero();
				coefs(1, iC, iC+1).setZero();
			}
		}

		if (l1 > 1)
		{
			auto inv_zeta = hInvWidthsAB();
			auto rho1 = this->rho1();
			
			coefs(2, 0, 0).colwise() = inv_zeta;
			coefs(2, 0, 1) = -rho1;
			coefs(2, 0, 2).setZero();

			for (int iC = 2; iC <= l2; iC += 2)
			{
				coefs(2, iC, 0) = coefs(0, iC, 0).colwise()*inv_zeta;
				coefs(2, iC, 1) = coefs(0, iC, 1).colwise()*inv_zeta
					- rho1*coefs(0, iC, 0);
				for (int m = 2; m < iC; ++m)
				{
					coefs(2, iC, m) = coefs(0, iC, m).colwise()*inv_zeta
						- rho1*coefs(0, iC, m-1)
						+ iC*inv_sum*coefs(1, iC-1, m-1);
				}
				coefs(2, iC, iC) = -rho1*coefs(0, iC, iC-1)
					+ iC*inv_sum*coefs(1, iC-1, iC-1);
				coefs(2, iC, iC+1).setZero();
				coefs(2, iC, iC+2).setZero();
			}

			// C_a,0,c,0
			for (int iA = 2; iA < l1; ++iA)
			{
				if (iA%2 == 1)
				{
					coefs(iA+1, 0, 0) = coefs(iA-1, 0, 0).colwise()*(iA*inv_zeta);
					for (int m = 1; m < iA-1; ++m)
					{
						coefs(iA+1, 0, m) = coefs(iA-1, 0, m).colwise()*(iA*inv_zeta)
							- iA*rho1*coefs(iA-1, 0, m-1);
					}
					coefs(iA+1, 0, iA-1) = -iA*rho1*coefs(iA-1, 0, iA-2);
					coefs(iA+1, 0, iA).setZero();
					coefs(iA+1, 0, iA+1).setZero();
				}

				for (int iC = 1+iA%2; iC <= l2; iC += 2)
				{
					coefs(iA+1, iC, 0) = coefs(iA-1, iC, 0).colwise()*(iA*inv_zeta);
					for (int m = 1; m < iA+iC-1; ++m)
					{
						coefs(iA+1, iC, m) = coefs(iA-1, iC, m).colwise()*(iA*inv_zeta)
							- iA*rho1*coefs(iA-1, iC, m-1)
							+ iC*inv_sum*coefs(iA, iC-1, m-1);
					}
					coefs(iA+1, iC, iA+iC-1) = -iA*rho1*coefs(iA-1, iC, iA+iC-2)
						+ iC*inv_sum*coefs(iA, iC-1, iA+iC-2);
					coefs(iA+1, iC, iA+iC).setZero();
					coefs(iA+1, iC, iA+iC+1).setZero();
				}
			}
		}
	}

	Cm.multiply(coefs, l1, l2);
}

void CGTOQuad::elecRepPrim1d_aacc(int i, FmCoefs& Cm) const
{
	int lA = this->lA(i), lB = this->lB(i),
		lC = this->lC(i), lD = this->lD(i);
	int l1 = lA+lB, l2 = lC+lD, lsum = l1+l2;

	if (lsum == 0)
	{
		return;
	}
	if (lsum == 1)
	{
		elecRepPrim1d_aacc_psss(i, Cm);
		return;
	}

	double dPQi = q().centerA(i) - p().centerA(i);
	if (lsum == 2)
	{
		if (l1 == 0)
		{
			Cm.multiplyRow(hInvWidthsCD(), -rho2(),
				((-dPQi * invWidthsSum()).colwise() * widthsAB()).square());
		}
		else if (l1 == 1)
		{
			Cm.multiply_noC0(0.5*invWidthsSum(),
				(invWidthsSum().square().colwise()
					* widthsAB()).rowwise() * widthsCD()
					* (-dPQi*dPQi));
		}
		else
		{
			Cm.multiplyCol(hInvWidthsAB(), -rho1(),
				((dPQi * invWidthsSum()).rowwise() * widthsCD()).square());
		}
		return;
	}

	EriCoefs coefs(l1, l2, p().size(), q().size());
	// C_0,0,0,0
	coefs(0, 0, 0).setOnes();
	if (l2 > 0)
	{
		auto inv_eta = hInvWidthsCD();
		auto dQWi = -dPQi * (invWidthsSum().colwise() * widthsAB());
		auto rho2 = this->rho2();

		// C_0,0,1,0
		coefs(0, 1, 0).setZero();
		coefs(0, 1, 1) = dQWi;
		// C_0,0,c,0
		for (int iC = 1; iC < l2; ++iC)
		{
			coefs(0, iC+1, 0) = coefs(0, iC-1, 0).rowwise()*(iC*inv_eta);
			for (int m = 1; m < iC; ++m)
			{
				coefs(0, iC+1, m) = coefs(0, iC-1, m).rowwise()*(iC*inv_eta)
					+ dQWi*coefs(0, iC, m-1)
					- iC*rho2*coefs(0, iC-1, m-1);
			}
			coefs(0, iC+1, iC) = dQWi*coefs(0, iC, iC-1)
				- iC*rho2*coefs(0, iC-1, iC-1);
			coefs(0, iC+1, iC+1) = dQWi*coefs(0, iC, iC);
		}
	}

	if (l1 > 0)
	{
		auto inv_zeta = hInvWidthsAB();
		auto dPWi = dPQi * (invWidthsSum().rowwise() * widthsCD());
		auto rho1 = this->rho1();
		auto inv_sum = 0.5*invWidthsSum();
		
		// C_1,0,c,0
		for (int iC = 0; iC <= l2; ++iC)
		{
			coefs(1, iC, 0).setZero();
			for (int m = 1; m <= iC; ++m)
			{
				coefs(1, iC, m) = dPWi*coefs(0, iC, m-1)
					+ iC*inv_sum*coefs(0, iC-1, m-1);
			}
			coefs(1, iC, iC+1) = dPWi*coefs(0, iC, iC);
		}

		// C_a,0,c,0
		for (int iA = 1; iA < l1; ++iA)
		{
			coefs(iA+1, 0, 0) = coefs(iA-1, 0, 0).colwise()*(iA*inv_zeta);
			for (int m = 1; m < iA; ++m)
			{
				coefs(iA+1, 0, m) = coefs(iA-1, 0, m).colwise()*(iA*inv_zeta)
					+ dPWi*coefs(iA, 0, m-1)
					- iA*rho1*coefs(iA-1, 0, m-1);
			}
			coefs(iA+1, 0, iA) = dPWi*coefs(iA, 0, iA-1)
				- iA*rho1*coefs(iA-1, 0, iA-1);
			coefs(iA+1, 0, iA+1) = dPWi*coefs(iA, 0, iA);

			for (int iC = 1; iC <= l2; ++iC)
			{
				coefs(iA+1, iC, 0) = coefs(iA-1, iC, 0).colwise()*(iA*inv_zeta);
				for (int m = 1; m < iA+iC; ++m)
				{
					coefs(iA+1, iC, m) = coefs(iA-1, iC, m).colwise()*(iA*inv_zeta)
						+ dPWi*coefs(iA, iC, m-1)
						- iA*rho1*coefs(iA-1, iC, m-1)
						+ iC*inv_sum*coefs(iA, iC-1, m-1);
				}
				coefs(iA+1, iC, iA+iC) = dPWi*coefs(iA, iC, iA+iC-1)
					- iA*rho1*coefs(iA-1, iC, iA+iC-1)
					+ iC*inv_sum*coefs(iA, iC-1, iA+iC-1);
				coefs(iA+1, iC, iA+iC+1) = dPWi*coefs(iA, iC, iA+iC);
			}
		}
	}

	Cm.multiply(coefs, l1, l2);
}

void CGTOQuad::elecRepPrim1d_abcc(int i, FmCoefs& Cm) const
{
	int lA = this->lA(i), lB = this->lB(i),
		lC = this->lC(i), lD = this->lD(i);
	int l1 = lA+lB, l2 = lC+lD, lsum = l1+l2;

	if (lsum == 0)
	{
		return;
	}
	else if (lsum == 1)
	{
		elecRepPrim1d_abcc_psss(i, Cm);
		return;
	}
	else if (lsum == 2)
	{
		if (l1 == 1)
			elecRepPrim1d_abcc_psps(i, Cm);
		else if (lA == 2 || lB == 2 || lC == 2 || lD == 2)
			elecRepPrim1d_abcc_dsss(i, Cm);
		else
			elecRepPrim1d_abcc_ppss(i, Cm);
		return;
	}

	double xA = centerA(i), xB = centerB(i);
	if (lA < lB)
	{
		std::swap(lA, lB);
		std::swap(xA, xB);
	}

	double dAB = xB - xA;

	EriCoefs coefs(l1, l2, p().size(), q().size());
	int idx = 0;
	// C_0,0,0,0
	coefs(idx++).setOnes();
	if (l2 > 0)
	{
		auto inv_eta = hInvWidthsCD();
		auto dQWi = dQW(i);
		auto rho2 = this->rho2();

		// C_0,0,1,0
		coefs(idx++).setZero();
		coefs(idx++) = dQWi;

		// C_0,0,c,0
		for (int iC = 1; iC < l2; ++iC)
		{
			coefs(idx) = coefs(idx-2*iC-1).rowwise()*(iC*inv_eta);
			++idx;
			for (int m = 1; m < iC; ++m, ++idx)
			{
				coefs(idx) = coefs(idx-2*iC-1).rowwise()*(iC*inv_eta)
					+ dQWi*coefs(idx-iC-2)
					- iC*rho2*coefs(idx-2*iC-2);
			}
			coefs(idx) = dQWi*coefs(idx-iC-2)
				- iC*rho2*coefs(idx-2*iC-2);
			++idx;
			coefs(idx) = dQWi*coefs(idx-iC-2);
			++idx;
		}
	}

	if (l1 > 0)
	{
		ColArray dAPi = dxP(i, xA);
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

	// C_a,b,c,d
	for (int iB = 1; iB <= lB; ++iB)
	{
		for (int iA = l1; iA >= lA+iB; --iA)
		{
			for (int m = 0; m < iA+l2; ++m)
				coefs(iA, l2, m) -= dAB*coefs(iA-1, l2, m);
		}
	}

	Cm.multiply(coefs, l1, l2);
}

void CGTOQuad::elecRepPrim1d_aacd(int i, FmCoefs& Cm) const
{
	int lA = this->lA(i), lB = this->lB(i),
		lC = this->lC(i), lD = this->lD(i);
	int l1 = lA+lB, l2 = lC+lD, lsum = l1+l2;

	if (lsum == 0)
	{
		return;
	}
	if (lsum == 1)
	{
		elecRepPrim1d_aacd_psss(i, Cm);
		return;
	}
	if (lsum == 2)
	{
		if (l1 == 1)
			elecRepPrim1d_aacd_psps(i, Cm);
		else if (lA == 2 || lB == 2 || lC == 2 || lD == 2)
			elecRepPrim1d_aacd_dsss(i, Cm);
		else
			elecRepPrim1d_aacd_ppss(i, Cm);
		return;
	}

	double xC = centerC(i), xD = centerD(i);
	if (lC < lD)
	{
		std::swap(lC, lD);
		std::swap(xC, xD);
	}

	double dCD = xD - xC;

	EriCoefs coefs(l1, l2, p().size(), q().size());
	int idx = 0;
	// C_0,0,0,0
	coefs(idx++).setOnes();
	if (l2 > 0)
	{
		auto dCQi = dxQ(i, xC);
		auto inv_eta = hInvWidthsCD();
		auto dQWi = dQW(i);
		auto rho2 = this->rho2();

		// C_0,0,1,0
		coefs(idx++).rowwise() = dCQi;
		coefs(idx++) = dQWi;
		
		// C_0,0,c,0
		for (int iC = 1; iC < l2; ++iC)
		{
			coefs(idx) = coefs(idx-iC-1).rowwise()*dCQi
				+ coefs(idx-2*iC-1).rowwise()*(iC*inv_eta);
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
		auto inv_zeta = hInvWidthsAB();
		Eigen::ArrayXXd dPWi = dPW(i);
		auto rho1 = this->rho1();
		auto inv_sum = 0.5 * invWidthsSum();

		int n2 = (l2+1)*(l2+2)/2;
		// C_1,0,c,0
		for (int iC = 0; iC <= l2; ++iC, ++n2)
		{
			coefs(idx++).setZero();
			for (int m = 1; m <= iC; ++m, ++idx)
			{
				coefs(idx) = dPWi*coefs(idx-n2-1)
					+ iC*inv_sum*coefs(idx-n2-iC-1);
			}
			coefs(idx) = dPWi*coefs(idx-n2-1);
			++idx;
		}

		// C_a,0,c,0
		for (int iA = 1; iA < l1; ++iA)
		{
			coefs(idx) = coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta);
			++idx;
			for (int m = 1; m < iA; ++m, ++idx)
			{
				coefs(idx) = coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta)
					+ dPWi*coefs(idx-n2-1)
					- iA*rho1*coefs(idx-2*n2+l2);
			}
			coefs(idx) = dPWi*coefs(idx-n2-1)
				- iA*rho1*coefs(idx-2*n2+l2);
			++idx;
			coefs(idx) = dPWi*coefs(idx-n2-1);
			++idx;
			++n2;

			for (int iC = 1; iC <= l2; ++iC, ++n2)
			{
				coefs(idx) = coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta);
				++idx;
				for (int m = 1; m < iA+iC; ++m, ++idx)
				{
					coefs(idx) = coefs(idx-2*n2+l2+1).colwise()*(iA*inv_zeta)
						+ dPWi*coefs(idx-n2-1)
						- iA*rho1*coefs(idx-2*n2+l2)
						+ iC*inv_sum*coefs(idx-n2-iC-iA-1);
				}
				coefs(idx) = dPWi*coefs(idx-n2-1)
					- iA*rho1*coefs(idx-2*n2+l2)
					+ iC*inv_sum*coefs(idx-n2-iC-iA-1);
				++idx;
				coefs(idx) = dPWi*coefs(idx-n2-1);
				++idx;
			}
		}
	}

	// C_a,0,c,d
	for (int iD = 1; iD <= lD; ++iD)
	{
		for (int iC = l2; iC >= lC+iD; --iC)
		{
			for (int m = 0; m < l1+iC; ++m)
				coefs(l1, iC, m) -= dCD*coefs(l1, iC-1, m);
		}
	}

	Cm.multiply(coefs, l1, l2);
}

void CGTOQuad::elecRepPrim1d_abcd(int i, FmCoefs& Cm) const
{
	int lA = this->lA(i), lB = this->lB(i),
		lC = this->lC(i), lD = this->lD(i);
	int l1 = lA+lB, l2 = lC+lD, lsum = l1+l2;

	if (lsum == 0)
	{
		return;
	}
	else if (lsum == 1)
	{
		elecRepPrim1d_abcd_psss(i, Cm);
		return;
	}
	else if (lsum == 2)
	{
		if (l1 == 1)
			elecRepPrim1d_abcd_psps(i, Cm);
		else if (lA == 2 || lB == 2 || lC == 2 || lD == 2)
			elecRepPrim1d_abcd_dsss(i, Cm);
		else
			elecRepPrim1d_abcd_ppss(i, Cm);
		return;
	}

	double xA = centerA(i), xB = centerB(i), xC = centerC(i), xD = centerD(i);
	if (lA < lB)
	{
		std::swap(lA, lB);
		std::swap(xA, xB);
	}
	if (lC < lD)
	{
		std::swap(lC, lD);
		std::swap(xC, xD);
	}

	double dAB = xB - xA;
	double dCD = xD - xC;
	
	EriCoefs coefs(l1, l2, p().size(), q().size());
	int idx = 0;
	// C_0,0,0,0
	coefs(idx++).setOnes();
	if (l2 > 0)
	{
		auto dCQi = dxQ(i, xC);
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
		auto dAPi = dxP(i, xA);
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

	// C_a,0,c,d
	for (int iA = lA; iA <= l1; ++iA)
	{
		for (int iD = 1; iD <= lD; ++iD)
		{
			for (int iC = l2; iC >= lC+iD; --iC)
			{
				for (int m = 0; m < iA+iC; ++m)
					coefs(iA, iC, m) -= dCD*coefs(iA, iC-1, m);
			}
		}
	}

	// C_a,b,c,d
	for (int iB = 1; iB <= lB; ++iB)
	{
		for (int iA = l1; iA >= lA+iB; --iA)
		{
			for (int m = 0; m < iA+l2; ++m)
				coefs(iA, l2, m) -= dAB*coefs(iA-1, l2, m);
		}
	}

	Cm.multiply(coefs, l1, l2);
}

double CGTOQuad::electronRepulsion_aaaa() const
{
	for (int i = 0; i < 3; i++)
	{
		if (lsum(i) % 2 == 1)
			return 0;
	}

	int ltot = lsum(0) + lsum(1) + lsum(2);
	FmCoefs Cm(ltot, p().size(), q().size());
	for (int i = 0; i < 3; i++)
		elecRepPrim1d_aaaa(i, Cm);

	int m = Cm.maxM();
	Eigen::ArrayXXd A = Cm[m] / (2*m+1);
	for (--m; m >= Cm.minM(); --m)
		A += Cm[m] / (2*m+1);

	return mulWeights(KKW() * A);
}

double CGTOQuad::electronRepulsion_aacc() const
{
	int ltot = lsum(0) + lsum(1) + lsum(2);
	FmCoefs Cm(ltot, p().size(), q().size());
	for (int i = 0; i < 3; i++)
		elecRepPrim1d_aacc(i, Cm);

	int m = Cm.maxM();
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
 	Eigen::ArrayXXd F = Fm(m);
	Eigen::ArrayXXd A = Cm[m] * F;
	for (--m; m >= Cm.minM(); --m)
	{
		F = (expmT + 2*T*F) / (2*m+1);
		A += Cm[m] * F;
	}

	return mulWeights(KKW() * A);
}

double CGTOQuad::electronRepulsion_abcc() const
{
	int ltot = lsum(0) + lsum(1) + lsum(2);
	FmCoefs Cm(ltot, p().size(), q().size());
	for (int i = 0; i < 3; i++)
		elecRepPrim1d_abcc(i, Cm);

	int m = Cm.maxM();
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
 	Eigen::ArrayXXd F = Fm(m);
	Cm[m] *= F;
	for (--m; m >= Cm.minM(); --m)
	{
		F = (expmT + 2*T*F) / (2*m+1);
		Cm[m] = Cm[m] * F + Cm[m+1];
	}

	return mulWeights(KKW() * Cm[Cm.minM()]);
}

double CGTOQuad::electronRepulsion_aacd() const
{
	int ltot = lsum(0) + lsum(1) + lsum(2);
	FmCoefs Cm(ltot, p().size(), q().size());
	for (int i = 0; i < 3; i++)
		elecRepPrim1d_aacd(i, Cm);

	int m = Cm.maxM();
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
 	Eigen::ArrayXXd F = Fm(m);
	Cm[m] *= F;
	for (--m; m >= Cm.minM(); --m)
	{
		F = (expmT + 2*T*F) / (2*m+1);
		Cm[m] = Cm[m] * F + Cm[m+1];
	}

	return mulWeights(KKW() * Cm[Cm.minM()]);
}

double CGTOQuad::electronRepulsion_abcd() const
{
	int ltot = lsum(0) + lsum(1) + lsum(2);
	FmCoefs Cm(ltot, p().size(), q().size());
	for (int i = 0; i < 3; i++)
		elecRepPrim1d_abcd(i, Cm);

	int m = Cm.maxM();
	const Eigen::ArrayXXd& T = this->T();
	const Eigen::ArrayXXd& expmT = this->expmT();
 	Eigen::ArrayXXd F = Fm(m);
	Cm[m] *= F;
	for (--m; m >= Cm.minM(); --m)
	{
		F = (expmT + 2*T*F) / (2*m+1);
		Cm[m] = Cm[m] * F + Cm[m+1];
	}

	return mulWeights(KKW() * Cm[Cm.minM()]);
}

double CGTOQuad::electronRepulsion() const
{
	switch (_shell_quad.positionSymmetry())
	{
		case CGTOShellQuad::POS_SYM_AAAA:
			return electronRepulsion_aaaa();
		case CGTOShellQuad::POS_SYM_AACC:
			return electronRepulsion_aacc();
		case CGTOShellQuad::POS_SYM_ABCC:
			return electronRepulsion_abcc();
		case CGTOShellQuad::POS_SYM_AACD:
			return electronRepulsion_aacd();
		default:
			return electronRepulsion_abcd();
	}
}

AbstractBFQuad* CGTOQuad::create(const AbstractBFPair& p,
	const AbstractBFPair& q, BFQuadPool& pool)
{
	try
	{
#ifdef DEBUG
		const CGTOPair& pp = dynamic_cast< const CGTOPair& >(p);
		const CGTOPair& qq = dynamic_cast< const CGTOPair& >(q);
#else
		const CGTOPair& pp = static_cast< const CGTOPair& >(p);
		const CGTOPair& qq = static_cast< const CGTOPair& >(q);
#endif
		return new(pool) CGTOQuad(pp, qq);
	}
	catch (const std::bad_cast&)
	{
		throw Li::Exception("Invalid basis function pair type");
	}
}
