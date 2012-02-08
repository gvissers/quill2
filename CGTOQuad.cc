#include "CGTOQuad.hh"
#include "boys.hh"

class EriCoefs
{
public:
	EriCoefs(int l1, int l2, int p_size, int q_size):
		_imax(l1+1), _jmax(l2+1), _p_size(p_size), _q_size(q_size),
		_C(p_size, (_imax*_jmax*(_imax+_jmax)*q_size)/2) {}

	int llIndex(int i, int j) const
	{
		return (i*_jmax*(i+_jmax) + j*(2*i+j+1)) / 2;
	}
	
	Eigen::Block<Eigen::ArrayXXd> operator()(int i, int j, int m)
	{
		return (*this)(llIndex(i,j) + m);
	}
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int i, int j,
		int m) const
	{
		return (*this)(llIndex(i,j) + m);
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
		return _C.block(0, (llIndex(i,j) + m)*_q_size, _p_size,
			count*_q_size);
	}
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int i, int j,
		int m, int count) const
	{
		return _C.block(0, (llIndex(i,j) + m)*_q_size, _p_size,
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
		_m(0)
	{
		block(0).setOnes();
	}

	int maxM() const { return _m; }

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
		return _C.block(0, m*_q_size, _p_size, count*_q_size);
	}
	const Eigen::Block<const Eigen::ArrayXXd> block(int m, int count=1) const
	{
		return _C.block(0, m*_q_size, _p_size, count*_q_size);
	}

	void multiplyCol(const CGTOQuad::ColArray& C0, const Eigen::ArrayXXd& C1)
	{
		block(_m+1) = block(_m) * C1;
		for (int m = _m; m > 0; --m)
			block(m) = block(m).colwise()*C0 + block(m-1)*C1;
		block(0).colwise() *= C0;
		++_m;
	}
	void multiplyRow(const CGTOQuad::RowArray& C0, const Eigen::ArrayXXd& C1)
	{
		block(_m+1) = block(_m) * C1;
		for (int m = _m; m > 0; --m)
			block(m) = block(m).rowwise()*C0 + block(m-1)*C1;
		block(0).rowwise() *= C0;
		++_m;
	}
	void multiply(const Eigen::ArrayXXd& C0, const Eigen::ArrayXXd& C1)
	{
		block(_m+1) = block(_m) * C1;
		for (int m = _m; m > 0; --m)
			block(m) = block(m)*C0 + block(m-1)*C1;
		block(0) *= C0;
		++_m;
	}
	void multiply_noC0(const Eigen::ArrayXXd& C1)
	{
		for (int m = _m; m >= 0; --m)
			block(m+1) = block(m)*C1;
		block(0).setZero();
		++_m;
	}
	void multiplyCol(const CGTOQuad::ColArray& C0, const Eigen::ArrayXXd& C1,
		const Eigen::ArrayXXd& C2)
	{
		if (_m == 0)
		{
			block(2) = block(0)*C2;
			block(1) = block(0)*C1;
			block(0).colwise() *= C0;
		}
		else
		{
			block(_m+2) = block(_m)*C2;
			block(_m+1) = block(_m)*C1 + block(_m-1)*C2;
			for (int m = _m; m > 1; --m)
			{
				block(m) = block(m).colwise()*C0 + block(m-1)*C1
					+ block(m-2)*C2;
			}
			block(1) = block(1).colwise()*C0 + block(0)*C1;
			block(0).colwise() *= C0;
		}
		_m += 2;
	}
	void multiplyRow(const CGTOQuad::RowArray& C0, const Eigen::ArrayXXd& C1,
		const Eigen::ArrayXXd& C2)
	{
		if (_m == 0)
		{
			block(2) = block(0)*C2;
			block(1) = block(0)*C1;
			block(0).rowwise() *= C0;
		}
		else
		{
			block(_m+2) = block(_m)*C2;
			block(_m+1) = block(_m)*C1 + block(_m-1)*C2;
			for (int m = _m; m > 1; --m)
			{
				block(m) = block(m).rowwise()*C0 + block(m-1)*C1
					+ block(m-2)*C2;
			}
			block(1) = block(1).rowwise()*C0 + block(0)*C1;
			block(0).rowwise() *= C0;
		}
		_m += 2;
	}
	void multiply(const Eigen::ArrayXXd& C0, const Eigen::ArrayXXd& C1,
		const Eigen::ArrayXXd& C2)
	{
		if (_m == 0)
		{
			block(2) = block(0)*C2;
			block(1) = block(0)*C1;
			block(0) *= C0;
		}
		else
		{
			block(_m+2) = block(_m)*C2;
			block(_m+1) = block(_m)*C1 + block(_m-1)*C2;
			for (int m = _m; m > 1; --m)
			{
				block(m) = block(m)*C0 + block(m-1)*C1
					+ block(m-2)*C2;
			}
			block(1) = block(1)*C0 + block(0)*C1;
			block(0) *= C0;
		}
		_m += 2;
	}
	void multiply_noC0(const Eigen::ArrayXXd& C1, const Eigen::ArrayXXd& C2)
	{
		if (_m == 0)
		{
			block(2) = block(0)*C2;
			block(1) = block(0)*C1;
			block(0).setZero();
		}
		else
		{
			block(_m+2) = block(_m)*C2;
			for (int m = _m; m > 0; --m)
				block(m+1) = block(m)*C1 + block(m-1)*C2;
			block(1) = block(0)*C1;
			block(0).setZero();
		}
		_m += 2;
	}
	void multiply(const EriCoefs& coefs, int l1, int l2)
	{
		int lsum = l1+l2;
		Eigen::ArrayXXd C(_p_size, _q_size);
		for (int m = lsum+_m; m >= 0; --m)
		{
			C.setZero();
			// i >= 0
			// i <= lsum
			// m-i >= 0 => i <= m
			// m-i <= _m => i >= m-_m
			for (int i = std::max(0, m-_m); i <= std::min(lsum, m); ++i)
				C += coefs(l1, l2, i) * block(m-i);
			block(m) = C;
		}
		_m += lsum;
	}

private:
	int _p_size;
	int _q_size;
	Eigen::ArrayXXd _C;
	int _m;
};

int n_abcd=0, n_abcc=0, n_aacd=0, n_aacc=0, n_aaaa=0;

CGTOQuad::CGTOQuad(const CGTOPair& p, const CGTOPair& q):
	AbstractBFQuad(p, q),
	_inv_widths_sum((widthsAB().replicate(1, q.size()).rowwise()
		+ widthsCD()).inverse())
{
	if (p.samePositionId())
	{
		if (q.samePositionId())
			_pos_sym = p.positionIdA() == q.positionIdA()
				? POS_SYM_AAAA : POS_SYM_AACC;
		else
			_pos_sym = POS_SYM_AACD;
	}
	else if (q.samePositionId())
	{
		_pos_sym = POS_SYM_ABCC;
	}
	else
	{
		_pos_sym = POS_SYM_ABCD;
	}

switch(_pos_sym)
{
	case POS_SYM_AAAA: n_aaaa++; break;
	case POS_SYM_AACC: n_aacc++; break;
	case POS_SYM_AACD: n_aacd++; break;
	case POS_SYM_ABCC: n_abcc++; break;
	case POS_SYM_ABCD: n_abcd++; break;
}
}

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

void CGTOQuad::elecRepPrim1d_abcd_psss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 1)
		throw Li::Exception("Not a (p,s,s,s) quad");
#endif

	if (lAB(i) == 1)
	{
		double x = lA(i) == 1 ? centerA(i) : centerB(i);
		Cm.multiplyCol(dxP(i, x),
			(dPQ(i) * invWidthsSum()).rowwise() * widthsCD());
	}
	else
	{
		double x = lC(i) == 1 ? centerC(i) : centerD(i);
		Cm.multiplyRow(dxQ(i, x),
			(dPQ(i) * invWidthsSum()).colwise() * (-widthsAB()));
	}
}

void CGTOQuad::elecRepPrim1d_abcd_ppss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| ((lA(i) != 1 || lB(i) != 1) && (lC(i) != 1 || lD(i) != 1)))
		throw Li::Exception("Not a (p,p,s,s) quad");
#endif

	ColArray inv_zeta = 0.5 * widthsAB().inverse();
	RowArray inv_eta = 0.5 * widthsCD().inverse();
	auto inv_ez = inv_zeta.replicate(1, q().size()).rowwise() + inv_eta;

	if (lAB(i) > 0)
	{
		double xA = centerA(i), xB = centerB(i), dAB = xB - xA;
		ColArray dAPi = dxP(i, xA);
		Eigen::ArrayXXd dPW = (dPQ(i) * invWidthsSum()).rowwise() * widthsCD();
		auto rho1 = inv_ez.inverse().colwise() * inv_zeta.square();

		Cm.multiplyCol(dAPi.square() + inv_zeta - dAB*dAPi,
			dPW.colwise()*(2*dAPi - dAB) - rho1,
			dPW.square());
	}
	else
	{
		double xC = centerC(i), xD = centerD(i), dCD = xD - xC;
		RowArray dCQi = dxQ(i, xC);
		Eigen::ArrayXXd dQW = (dPQ(i) * invWidthsSum()).colwise() * (-widthsAB());
		auto rho2 = inv_ez.inverse().rowwise() * inv_eta.square();

		Cm.multiplyRow(dCQi.square() + inv_eta - dCD*dCQi,
			dQW.rowwise()*(2*dCQi - dCD) - rho2,
			dQW.square());
	}
}

void CGTOQuad::elecRepPrim1d_abcd_dsss(int i, FmCoefs& Cm) const
{
#ifdef DEBUG
	if (lsum(i) != 2
		|| (lA(i) != 2 && lB(i) != 2 && lC(i) != 2 && lD(i) != 2))
		throw Li::Exception("Not a (d,s,s,s) quad");
#endif

	ColArray inv_zeta = 0.5 * widthsAB().inverse();
	RowArray inv_eta = 0.5 * widthsCD().inverse();
	auto inv_ez = inv_zeta.replicate(1, q().size()).rowwise() + inv_eta;

	if (lAB(i) > 0)
	{
		double x = lA(i) > 0 ? centerA(i) : centerB(i);
		ColArray dAPi = dxP(i, x);
		Eigen::ArrayXXd dPW = (dPQ(i) * invWidthsSum()).rowwise() * widthsCD();
		auto rho1 = inv_zeta.square().replicate(1, q().size()) / inv_ez;

		Cm.multiplyCol(dAPi.square() + inv_zeta,
			dPW.colwise()*(2*dAPi) - rho1,
			dPW.square());
	}
	else
	{
		double x = lC(i) > 0 ? centerC(i) : centerD(i);
		RowArray dCQi = dxQ(i, x);
		Eigen::ArrayXXd dQW = (dPQ(i) * invWidthsSum()).colwise() * (-widthsAB());
		auto rho2 = inv_ez.inverse().rowwise() * inv_eta.square();

		Cm.multiplyRow(dCQi.square() + inv_eta,
			dQW.rowwise()*(2*dCQi) - rho2,
			dQW.square());
	}
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

	ColArray dAPi = dxP(i, xA);
	RowArray dCQi = dxQ(i, xC);
	Eigen::ArrayXXd dPQi = dPQ(i);
	Eigen::ArrayXXd dPW = (dPQi * invWidthsSum()).rowwise() * widthsCD();
	Eigen::ArrayXXd dQW = (dPQi * invWidthsSum()).colwise() * (-widthsAB());

	Cm.multiply(dAPi.matrix()*dCQi.matrix(),
		dQW.colwise()*dAPi + dPW.rowwise()*dCQi + 0.5*invWidthsSum(),
		dPW*dQW);
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
		return;
	
	Eigen::ArrayXXd inv_zeta = 0.5 * widthsAB().inverse().replicate(1, q().size());
	Eigen::ArrayXXd inv_eta = 0.5 * widthsCD().inverse().replicate(p().size(), 1);
	Eigen::ArrayXXd inv_ez = inv_eta + inv_zeta;
	if (lsum == 2)
	{
		if (l1 == 0)
			Cm.multiply(inv_eta, -inv_eta.square() / inv_ez);
		else if (l1 == 1)
			Cm.multiply_noC0(0.5*invWidthsSum());
		else
 			Cm.multiply(inv_zeta, -inv_zeta.square() / inv_ez);
		return;
	}
	if (lsum == 4)
	{
		if (l1 == 0)
		{
			const Eigen::ArrayXXd rho2 = -inv_eta.square() / inv_ez;
			Cm.multiply(3*inv_eta.square(), 6*inv_eta*rho2,
				 3*rho2.square());
		}
		else if (l1 == 1)
		{
			const Eigen::ArrayXXd rho2 = -inv_eta.square() / inv_ez;
			Cm.multiply_noC0(1.5*invWidthsSum()*inv_eta,
				1.5*invWidthsSum()*rho2);
		}		
		else if (l1 == 2)
		{
			const Eigen::ArrayXXd rho1 = -inv_zeta.square() / inv_ez;
			const Eigen::ArrayXXd rho2 = -inv_eta.square() / inv_ez;
			Cm.multiply(inv_zeta*inv_eta,
				inv_zeta*rho2 + inv_eta*rho1,
				rho1*rho2 + 0.5*invWidthsSum().square());
		}
		else if (l1 == 3)
		{
			const Eigen::ArrayXXd rho1 = -inv_zeta.square() / inv_ez;
			Cm.multiply_noC0(1.5*invWidthsSum()*inv_zeta,
				1.5*invWidthsSum()*rho1);
		}
		else
		{
			const Eigen::ArrayXXd rho1 = -inv_zeta.square() / inv_ez;
			Cm.multiply(3*inv_zeta.square(), 6*inv_zeta*rho1,
				3*rho1.square());
		}
		return;
	}

	Eigen::ArrayXXd zero = Eigen::ArrayXXd::Zero(p().size(), q().size());

	EriCoefs coefs(l1, l2, p().size(), q().size());
	// C_0,0,0,0
	coefs(0, 0, 0).setOnes();

	// C_0,0,c,0
	if (l2 > 1)
	{
		coefs(0, 2, 0) = inv_eta;
		coefs(0, 2, 1) = -inv_eta.square() / inv_ez;
		coefs(0, 2, 2).setZero();

		auto rho2 = coefs(0, 2, 1);
		for (int iC = 3; iC < l2; iC += 2)
		{
			coefs(0, iC+1, 0) = iC*inv_eta*coefs(0, iC-1, 0);
			for (int m = 1; m < iC-1; ++m)
			{
				coefs(0, iC+1, m) = iC*inv_eta*coefs(0, iC-1, m)
					+ iC*rho2*coefs(0, iC-1, m-1);
			}
			coefs(0, iC+1, iC-1) = iC*rho2*coefs(0, iC-1, iC-2);
			coefs(0, iC+1, iC).setZero();
			coefs(0, iC+1, iC+1).setZero();
		}
	}

	if (l1 > 0)
	{
		Eigen::ArrayXXd inv_sum = 0.5 * invWidthsSum();

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
			coefs(2, 0, 0) = inv_zeta;
			coefs(2, 0, 1) = -inv_zeta.square() / inv_ez;
			coefs(2, 0, 2).setZero();

			auto rho1 = coefs(2, 0, 1);
			for (int iC = 2; iC <= l2; iC += 2)
			{
				coefs(2, iC, 0) = inv_zeta*coefs(0, iC, 0);
				coefs(2, iC, 1) = inv_zeta*coefs(0, iC, 1)
					+ rho1*coefs(0, iC, 0);
				for (int m = 2; m < iC; ++m)
				{
					coefs(2, iC, m) = inv_zeta*coefs(0, iC, m)
						+ rho1*coefs(0, iC, m-1)
						+ iC*inv_sum*coefs(1, iC-1, m-1);
				}
				coefs(2, iC, iC) = rho1*coefs(0, iC, iC-1)
					+ iC*inv_sum*coefs(1, iC-1, iC-1);
				coefs(2, iC, iC+1).setZero();
				coefs(2, iC, iC+2).setZero();
			}

			// C_a,0,c,0
			for (int iA = 2; iA < l1; ++iA)
			{
				if (iA%2 == 1)
				{
					coefs(iA+1, 0, 0) = iA*inv_zeta*coefs(iA-1, 0, 0);
					for (int m = 1; m < iA-1; ++m)
					{
						coefs(iA+1, 0, m) = iA*inv_zeta*coefs(iA-1, 0, m)
							+ iA*rho1*coefs(iA-1, 0, m-1);
					}
					coefs(iA+1, 0, iA-1) = iA*rho1*coefs(iA-1, 0, iA-2);
					coefs(iA+1, 0, iA).setZero();
					coefs(iA+1, 0, iA+1).setZero();
				}

				for (int iC = 1+iA%2; iC <= l2; iC += 2)
				{
					coefs(iA+1, iC, 0) = iA*inv_zeta*coefs(iA-1, iC, 0);
					for (int m = 1; m < iA+iC-1; ++m)
					{
						coefs(iA+1, iC, m) = iA*inv_zeta*coefs(iA-1, iC, m)
							+ iA*rho1*coefs(iA-1, iC, m-1)
							+ iC*inv_sum*coefs(iA, iC-1, m-1);
					}
					coefs(iA+1, iC, iA+iC-1) = iA*rho1*coefs(iA-1, iC, iA+iC-2)
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
	else if (lsum == 1)
	{
		elecRepPrim1d_aacc_psss(i, Cm);
		return;
	}

	ColArray inv_zeta = 0.5 * widthsAB().inverse();
	RowArray inv_eta = 0.5 * widthsCD().inverse();
	double dPQi = q().centerA(i) - p().centerA(i);
	if (lsum == 2)
	{
		if (l1 == 0)
		{
			Cm.multiplyRow(inv_eta,
				-((invWidthsSum().colwise() * widthsAB()).rowwise() * inv_eta),
				((-dPQi * invWidthsSum()).colwise() * widthsAB()).square());
		}
		else if (l1 == 1)
		{
			Cm.multiply_noC0(0.5*invWidthsSum(),
				-dPQi*dPQi*invWidthsSum().square()*widthsProduct());
		}
		else if (l1 == 2)
		{
			Cm.multiplyCol(inv_zeta,
				-((invWidthsSum().colwise() * inv_zeta).rowwise() * widthsCD()),
				((dPQi * invWidthsSum()).rowwise() * widthsCD()).square());
		}
		return;
	}

	EriCoefs coefs(l1, l2, p().size(), q().size());
	// C_0,0,0,0
	coefs(0, 0, 0).setOnes();
	if (l2 > 0)
	{
		Eigen::ArrayXXd dQW = -dPQi * (invWidthsSum().colwise() * widthsAB());
		Eigen::ArrayXXd rho2 = (invWidthsSum().colwise() * widthsAB())
			.rowwise() * inv_eta;

		// C_0,0,1,0
		coefs(0, 1, 0).setZero();
		coefs(0, 1, 1) = dQW;
		// C_0,0,c,0
		for (int iC = 1; iC < l2; ++iC)
		{
			coefs(0, iC+1, 0) = coefs(0, iC-1, 0).rowwise()*(iC*inv_eta);
			for (int m = 1; m < iC; ++m)
			{
				coefs(0, iC+1, m) = coefs(0, iC-1, m).rowwise()*(iC*inv_eta)
					+ dQW*coefs(0, iC, m-1)
					- iC*rho2*coefs(0, iC-1, m-1);
			}
			coefs(0, iC+1, iC) = dQW*coefs(0, iC, iC-1)
				- iC*rho2*coefs(0, iC-1, iC-1);
			coefs(0, iC+1, iC+1) = dQW*coefs(0, iC, iC);
		}
	}

	if (l1 > 0)
	{
		Eigen::ArrayXXd dPW = dPQi * (invWidthsSum().rowwise() * widthsCD());
		Eigen::ArrayXXd rho1 = (invWidthsSum().colwise() * inv_zeta)
			.rowwise() * widthsCD();
		Eigen::ArrayXXd inv_sum = 0.5*invWidthsSum();
		
		// C_1,0,c,0
		for (int iC = 0; iC <= l2; ++iC)
		{
			coefs(1, iC, 0).setZero();
			for (int m = 1; m <= iC; ++m)
			{
				coefs(1, iC, m) = dPW*coefs(0, iC, m-1)
					+ iC*inv_sum*coefs(0, iC-1, m-1);
			}
			coefs(1, iC, iC+1) = dPW*coefs(0, iC, iC);
		}

		// C_a,0,c,0
		for (int iA = 1; iA < l1; ++iA)
		{
			coefs(iA+1, 0, 0) = coefs(iA-1, 0, 0).colwise()*(iA*inv_zeta);
			for (int m = 1; m < iA; ++m)
			{
				coefs(iA+1, 0, m) = coefs(iA-1, 0, m).colwise()*(iA*inv_zeta)
					+ dPW*coefs(iA, 0, m-1)
					- iA*rho1*coefs(iA-1, 0, m-1);
			}
			coefs(iA+1, 0, iA) = dPW*coefs(iA, 0, iA-1)
				- iA*rho1*coefs(iA-1, 0, iA-1);
			coefs(iA+1, 0, iA+1) = dPW*coefs(iA, 0, iA);

			for (int iC = 1; iC <= l2; ++iC)
			{
				coefs(iA+1, iC, 0) = coefs(iA-1, iC, 0).colwise()*(iA*inv_zeta);
				for (int m = 1; m < iA+iC; ++m)
				{
					coefs(iA+1, iC, m) = coefs(iA-1, iC, m).colwise()*(iA*inv_zeta)
						+ dPW*coefs(iA, iC, m-1)
						- iA*rho1*coefs(iA-1, iC, m-1)
						+ iC*inv_sum*coefs(iA, iC-1, m-1);
				}
				coefs(iA+1, iC, iA+iC) = dPW*coefs(iA, iC, iA+iC-1)
					- iA*rho1*coefs(iA-1, iC, iA+iC-1)
					+ iC*inv_sum*coefs(iA, iC-1, iA+iC-1);
				coefs(iA+1, iC, iA+iC+1) = dPW*coefs(iA, iC, iA+iC);
			}
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
	Eigen::ArrayXXd dPQi = dPQ(i);

	EriCoefs coefs(l1, l2, p().size(), q().size());
	// C_0,0,0,0
	coefs(0, 0, 0).setOnes();
	if (l2 > 0)
	{
		RowArray dCQi = dxQ(i, xC);
		RowArray inv_eta = 0.5 * widthsCD().inverse();
		Eigen::ArrayXXd dQW = (invWidthsSum() * dPQi).colwise() * (-widthsAB());
		Eigen::ArrayXXd rho2 = (invWidthsSum().colwise() * widthsAB())
			.rowwise() * inv_eta;

		// C_0,0,1,0
		coefs(0, 1, 0).rowwise() = dCQi;
		coefs(0, 1, 1) = dQW;
		// C_0,0,c,0
		for (int iC = 1; iC < l2; ++iC)
		{
			coefs(0, iC+1, 0) = coefs(0, iC, 0).rowwise()*dCQi
				+ coefs(0, iC-1, 0).rowwise()*(iC*inv_eta);
			for (int m = 1; m < iC; ++m)
			{
				coefs(0, iC+1, m) = coefs(0, iC, m).rowwise()*dCQi
					+ coefs(0, iC-1, m).rowwise()*(iC*inv_eta)
					+ dQW*coefs(0, iC, m-1)
					- iC*rho2*coefs(0, iC-1, m-1);
			}
			coefs(0, iC+1, iC) = coefs(0, iC, iC).rowwise()*dCQi
				+ dQW*coefs(0, iC, iC-1)
				- iC*rho2*coefs(0, iC-1, iC-1);
			coefs(0, iC+1, iC+1) = dQW*coefs(0, iC, iC);
		}
	}

	if (l1 > 0)
	{
		ColArray dAPi = dxP(i, xA);
		ColArray inv_zeta = 0.5 * widthsAB().inverse();
		Eigen::ArrayXXd dPW = (invWidthsSum() * dPQi).rowwise() * widthsCD();
		Eigen::ArrayXXd rho1 = (invWidthsSum().colwise() * inv_zeta)
			.rowwise() * widthsCD();
		Eigen::ArrayXXd inv_sum = 0.5 * invWidthsSum();
		
		// C_1,0,c,0
		for (int iC = 0; iC <= l2; ++iC)
		{
			coefs(1, iC, 0) = coefs(0, iC, 0).colwise()*dAPi;
			for (int m = 1; m <= iC; ++m)
			{
				coefs(1, iC, m) = coefs(0, iC, m).colwise()*dAPi
					+ dPW*coefs(0, iC, m-1)
					+ iC*inv_sum*coefs(0, iC-1, m-1);
			}
			coefs(1, iC, iC+1) = dPW*coefs(0, iC, iC);
		}

		// C_a,0,c,0
		for (int iA = 1; iA < l1; ++iA)
		{
			coefs(iA+1, 0, 0) = coefs(iA, 0, 0).colwise()*dAPi
				+ coefs(iA-1, 0, 0).colwise()*(iA*inv_zeta);
			for (int m = 1; m < iA; ++m)
			{
				coefs(iA+1, 0, m) = coefs(iA, 0, m).colwise()*dAPi
					+ coefs(iA-1, 0, m).colwise()*(iA*inv_zeta)
					+ dPW*coefs(iA, 0, m-1)
					- iA*rho1*coefs(iA-1, 0, m-1);
			}
			coefs(iA+1, 0, iA) = coefs(iA, 0, iA).colwise()*dAPi
				+ dPW*coefs(iA, 0, iA-1)
				- iA*rho1*coefs(iA-1, 0, iA-1);
			coefs(iA+1, 0, iA+1) = dPW*coefs(iA, 0, iA);

			for (int iC = 1; iC <= l2; ++iC)
			{
				coefs(iA+1, iC, 0) = coefs(iA, iC, 0).colwise()*dAPi
					+ coefs(iA-1, iC, 0).colwise()*(iA*inv_zeta);
				for (int m = 1; m < iA+iC; ++m)
				{
					coefs(iA+1, iC, m) = coefs(iA, iC, m).colwise()*dAPi
						+ coefs(iA-1, iC, m).colwise()*(iA*inv_zeta)
						+ dPW*coefs(iA, iC, m-1)
						- iA*rho1*coefs(iA-1, iC, m-1)
						+ iC*inv_sum*coefs(iA, iC-1, m-1);
				}
				coefs(iA+1, iC, iA+iC) = coefs(iA, iC, iA+iC).colwise()*dAPi
					+ dPW*coefs(iA, iC, iA+iC-1)
					- iA*rho1*coefs(iA-1, iC, iA+iC-1)
					+ iC*inv_sum*coefs(iA, iC-1, iA+iC-1);
				coefs(iA+1, iC, iA+iC+1) = dPW*coefs(iA, iC, iA+iC);
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

double CGTOQuad::electronRepulsion_aaaa_ssss() const
{
	return weightsAB().transpose()
		* (KK() * invWidthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion_aacc_ssss() const
{
	Eigen::ArrayXXd T = (p().centerA() - q().centerA()).squaredNorm()
		* widthsReduced();
	return weightsAB().transpose()
		* (KK() * Fm(0, T) * invWidthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion_abcd_ssss() const
{
	Eigen::ArrayXXd T = (dPQ(0).square() + dPQ(1).square() + dPQ(2).square())
		* widthsReduced();
	return weightsAB().transpose()
		* (KK() * Fm(0, T) * invWidthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion_abcd_psss() const
{
	Eigen::ArrayXXd T = (dPQ(0).square() + dPQ(1).square() + dPQ(2).square())
		* widthsReduced();
	Eigen::ArrayXXd F1 = Fm(1, T);
	Eigen::ArrayXXd F0 = (-T).exp() + 2*T*F1;
	for (int i = 0; i < 3; i++)
	{
		if (lAB(i) == 1)
		{
			double x = lA(i) == 1 ? centerA(i) : centerB(i);
			F0.colwise() *= dxP(i, x);
			F1 *= dPQ(i).rowwise() * widthsCD() * invWidthsSum();
			break;
		}
		if (lCD(i) == 1)
		{
			double x = lC(i) == 1 ? centerC(i) : centerD(i);
			F0.rowwise() *= dxQ(i, x);
			F1 *= dPQ(i).colwise() * (-widthsAB()) * invWidthsSum();
			break;
		}
	}

	return weightsAB().transpose()
		* (KK() * (F1+F0) * invWidthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion_aacc_psss() const
{
	Eigen::ArrayXXd T = (p().centerA() - q().centerA()).squaredNorm()
		* widthsReduced();
	Eigen::ArrayXXd C = KK() * Fm(1, T) * invWidthsSum() * invWidthsSum().sqrt();
	for (int i = 0; i < 3; i++)
	{
		if (lAB(i) == 1)
		{
			C.rowwise() *= widthsCD() * (q().centerA(i)-p().centerA(i));
			break;
		}
		if (lCD(i) == 1)
		{
			C.colwise() *= widthsAB() * (p().centerA(i)-q().centerA(i));
			break;
		}
	}

	return weightsAB().transpose() * C.matrix() * weightsCD();
}

double CGTOQuad::electronRepulsion_aaaa() const
{
	const Eigen::Vector3i& lsA = p().f().ls();
	const Eigen::Vector3i& lsB = p().g().ls();
	const Eigen::Vector3i& lsC = q().f().ls();
	const Eigen::Vector3i& lsD = q().g().ls();
	Eigen::Vector3i ls = lsA + lsB + lsC + lsD;
	int lsum = ls.sum();

	if (lsum == 0)
		return electronRepulsion_aaaa_ssss();

	for (int i = 0; i < 3; i++)
	{
		if (ls[i] % 2 == 1)
			return 0;
	}

	FmCoefs Cm(lsum, p().size(), q().size());
	for (int i = 0; i < 3; i++)
		elecRepPrim1d_aaaa(i, Cm);

	int m = Cm.maxM();
	Eigen::ArrayXXd A = Cm[m] / (2*m+1);
	for (--m; m >= 0; --m)
		A += Cm[m] / (2*m+1);

	return weightsAB().transpose()
		* (KK() * A * invWidthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion_aacc() const
{
	const Eigen::Vector3i& lsA = p().f().ls();
	const Eigen::Vector3i& lsB = p().g().ls();
	const Eigen::Vector3i& lsC = q().f().ls();
	const Eigen::Vector3i& lsD = q().g().ls();
	Eigen::Vector3i ls = lsA + lsB + lsC + lsD;

	int lsum = ls.sum();

	if (lsum == 0)
		return electronRepulsion_aacc_ssss();
	else if (lsum == 1)
		return electronRepulsion_aacc_psss();

	FmCoefs Cm(lsum, p().size(), q().size());
	Eigen::ArrayXXd T = (p().centerA() - q().centerA()).squaredNorm()
		* widthsReduced();
	for (int i = 0; i < 3; i++)
		elecRepPrim1d_aacc(i, Cm);

	int m = Cm.maxM();
	Eigen::ArrayXXd F = Fm(m, T);
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd A = Cm[m] * F;
	for (--m; m >= 0; --m)
	{
		F = (expmT + 2*T*F) / (2*m+1);
		A += Cm[m] * F;
	}

	return weightsAB().transpose()
		* (KK() * A * invWidthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion_abcd() const
{
	const Eigen::Vector3i& lsA = p().f().ls();
	const Eigen::Vector3i& lsB = p().g().ls();
	const Eigen::Vector3i& lsC = q().f().ls();
	const Eigen::Vector3i& lsD = q().g().ls();
	Eigen::Vector3i ls = lsA + lsB + lsC + lsD;

	int lsum = ls.sum();

	if (lsum == 0)
		return electronRepulsion_abcd_ssss();
	else if (lsum == 1)
		return electronRepulsion_abcd_psss();

	FmCoefs Cm(lsum, p().size(), q().size());
	Eigen::ArrayXXd T = Eigen::ArrayXXd::Zero(p().size(), q().size());
	for (int i = 0; i < 3; i++)
	{
		elecRepPrim1d_abcd(i, Cm);
		T += dPQ(i).square();
	}
	T *= widthsReduced();

	int m = Cm.maxM();
	Eigen::ArrayXXd F = Fm(lsum, T);
	Eigen::ArrayXXd expmT = (-T).exp();
	Eigen::ArrayXXd A = Cm[m] * F;
	for (--m; m >= 0; --m)
	{
		F = (expmT + 2*T*F) / (2*m+1);
		A += Cm[m] * F;
	}

	return weightsAB().transpose()
		* (KK() * A * invWidthsSum().sqrt()).matrix()
		* weightsCD();
}

double CGTOQuad::electronRepulsion() const
{
	switch (_pos_sym)
	{
		case POS_SYM_AAAA:
			return electronRepulsion_aaaa();
		case POS_SYM_AACC:
			return electronRepulsion_aacc();
		default:
			return electronRepulsion_abcd();
	}
}