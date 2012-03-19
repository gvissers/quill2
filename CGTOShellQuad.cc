#include <Exception.hh>
#include "CGTOShellQuad.hh"

CGTOShellQuad::CGTOShellQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD):
	_pAB(pAB), _pCD(pCD),
	_pos_sym(symmetry(pAB, pCD)),
	_lAB(pAB.lsum()), _lCD(pCD.lsum()), _lsum(_lAB+_lCD),
	_data(_pAB.size(), _pCD.size(), 10),
	_ints((_lAB+1)*(_lAB+1)*(_lAB+1), (_lCD+1)*(_lCD+1)*(_lCD+1)),
	_have_eri(false)
{
	_data[0] = (widthsAB().replicate(1, pCD.size()).rowwise()
		+ widthsCD()).inverse();
	_data[1] = Q(0).replicate(_pAB.size(), 1).colwise() - P(0);
	_data[2] = Q(1).replicate(_pAB.size(), 1).colwise() - P(1);
	_data[3] = Q(2).replicate(_pAB.size(), 1).colwise() - P(2);
	if (_lAB > 0)
	{
		_data[4] = (dPQ(0) * invWidthsSum()).rowwise() * widthsCD();
		_data[5] = (dPQ(1) * invWidthsSum()).rowwise() * widthsCD();
		_data[6] = (dPQ(2) * invWidthsSum()).rowwise() * widthsCD();
	}
	if (_lCD > 0)
	{
		_data[7] = (dPQ(0) * invWidthsSum()).colwise() * (-widthsAB());
		_data[8] = (dPQ(1) * invWidthsSum()).colwise() * (-widthsAB());
		_data[9] = (dPQ(2) * invWidthsSum()).colwise() * (-widthsAB());
	}
}

//double CGTOShellQuad::mulWeights(const Eigen::ArrayXXd& C) const
//{
//	return ((C.colwise() * ColArray::MapAligned(_pAB.weights().data(), _pAB.size())).colwise().sum()
//		* RowArray::MapAligned(_pCD.weights().data(), _pCD.size())).sum();
//}

void CGTOShellQuad::elecRepPrim1d_abcd(int i, EriCoefs& coefs) const
{
	int idx = 0;
	// C_0,0,0,0
	coefs[idx++].setOnes();
	if (_lsum == 1)
	{
		if (_lAB == 1)
		{
			// C_1,0,0,0
			coefs[idx++].colwise() = dAP(i);
			coefs[idx] = dPW(i);
		}
		else
		{
			// C_0,0,1,0
			coefs[idx++].rowwise() = dCQ(i);
			coefs[idx] = dQW(i);
		}
		return;
	}
	if (_lsum == 2)
	{
		if (_lAB == 2)
		{
			auto dAPi = dAP(i);
			auto dPWi = dPW(i);

			// C_1,0,0,0
			coefs[idx++].colwise() = dAPi;
			coefs[idx++] = dPWi;
			// C_2,0,0,0
			coefs[idx++].colwise() = dAPi.square() + hInvWidthsAB();
			coefs[idx++] = dPWi.colwise()*dAPi*2 - rho1();
			coefs[idx] = dPWi.square();
		}
		else if (_lAB == 1)
		{
			auto dCQi = dCQ(i);
			auto dQWi = dQW(i);
			auto dAPi = dAP(i);
			auto dPWi = dPW(i);

			// C_0,0,1,0
			coefs[idx++].rowwise() = dCQi;
			coefs[idx++] = dQWi;
			// C_1,0,0,0
			coefs[idx++].colwise() = dAPi;
			coefs[idx++] = dPWi;
			// C_1,0,1,0
			coefs[idx++] = dCQi.replicate(_pAB.size(), 1).colwise()*dAPi;
			coefs[idx++] = dQWi.colwise()*dAPi + dPWi.rowwise()*dCQi
				+ 0.5*invWidthsSum();
			coefs[idx] = dPWi*dQWi;
		}
		else
		{
			auto dCQi = dCQ(i);
			auto dQWi = dQW(i);

			// C_0,0,1,0
			coefs[idx++].rowwise() = dCQi;
			coefs[idx++] = dQWi;
			// C_0,0,2,0
			coefs[idx++].rowwise() = dCQi.square() + hInvWidthsCD();
			coefs[idx++] = dQWi.rowwise()*dCQi*2 - rho2();
			coefs[idx] = dQWi.square();
		}
		return;
	}

	if (_lCD > 0)
	{
		auto dCQi = dCQ(i);
		auto inv_eta = hInvWidthsCD();
		auto dQWi = dQW(i);
		auto rho2 = this->rho2();

		// C_0,0,1,0
		coefs[idx++].rowwise() = dCQi;
		coefs[idx++] = dQWi;

		// C_0,0,c,0
		for (int iC = 1; iC < _lCD; ++iC)
		{
			coefs[idx].rowwise() = coefs[idx-iC-1].row(0)*dCQi
				+ coefs[idx-2*iC-1].row(0)*(iC*inv_eta);
			++idx;
			for (int m = 1; m < iC; ++m, ++idx)
			{
				coefs[idx] = coefs[idx-iC-1].rowwise()*dCQi
					+ coefs[idx-2*iC-1].rowwise()*(iC*inv_eta)
					+ dQWi*coefs[idx-iC-2]
					- iC*rho2*coefs[idx-2*iC-2];
			}
			coefs[idx] = coefs[idx-iC-1].rowwise()*dCQi
				+ dQWi*coefs[idx-iC-2]
				- iC*rho2*coefs[idx-2*iC-2];
			++idx;
			coefs[idx] = dQWi*coefs[idx-iC-2];
			++idx;
		}
	}

	if (_lAB > 0)
	{
		auto dAPi = dAP(i);
		auto inv_zeta = hInvWidthsAB();
		auto dPWi = dPW(i);
		auto rho1 = this->rho1();
		auto inv_sum = 0.5 * invWidthsSum();

		int n2 = (_lCD+1)*(_lCD+2)/2;
		// C_1,0,c,0
		for (int iC = 0; iC <= _lCD; ++iC, ++n2)
		{
			coefs[idx] = coefs[idx-n2].colwise()*dAPi;
			++idx;
			for (int m = 1; m <= iC; ++m, ++idx)
			{
				coefs[idx] = coefs[idx-n2].colwise()*dAPi
					+ dPWi*coefs[idx-n2-1]
					+ iC*inv_sum*coefs[idx-n2-iC-1];
			}
			coefs[idx] = dPWi*coefs[idx-n2-1];
			++idx;
		}

		// C_a,0,c,0
		for (int iA = 1; iA < _lAB; ++iA)
		{
			coefs[idx] = coefs[idx-n2].colwise()*dAPi
				+ coefs[idx-2*n2+_lCD+1].colwise()*(iA*inv_zeta);
			++idx;
			for (int m = 1; m < iA; ++m, ++idx)
			{
				coefs[idx] = coefs[idx-n2].colwise()*dAPi
					+ coefs[idx-2*n2+_lCD+1].colwise()*(iA*inv_zeta)
					+ dPWi*coefs[idx-n2-1]
					- iA*rho1*coefs[idx-2*n2+_lCD];
			}
			coefs[idx] = coefs[idx-n2].colwise()*dAPi
				+ dPWi*coefs[idx-n2-1]
				- iA*rho1*coefs[idx-2*n2+_lCD];
			++idx;
			coefs[idx] = dPWi*coefs[idx-n2-1];
			++idx;
			++n2;

			for (int iC = 1; iC <= _lCD; ++iC, ++n2)
			{
				coefs[idx] = coefs[idx-n2].colwise()*dAPi
					+ coefs[idx-2*n2+_lCD+1].colwise()*(iA*inv_zeta);
				++idx;
				for (int m = 1; m < iA+iC; ++m, ++idx)
				{
					coefs[idx] = coefs[idx-n2].colwise()*dAPi
						+ coefs[idx-2*n2+_lCD+1].colwise()*(iA*inv_zeta)
						+ dPWi*coefs[idx-n2-1]
						- iA*rho1*coefs[idx-2*n2+_lCD]
						+ iC*inv_sum*coefs[idx-n2-iC-iA-1];
				}
				coefs[idx] = coefs[idx-n2].colwise()*dAPi
					+ dPWi*coefs[idx-n2-1]
					- iA*rho1*coefs[idx-2*n2+_lCD]
					+ iC*inv_sum*coefs[idx-n2-iC-iA-1];
				++idx;
				coefs[idx] = dPWi*coefs[idx-n2-1];
				++idx;
			}
		}
	}
}

double CGTOShellQuad::eri_xx(int lx1, int lx2, const EriCoefs& Cx, const Fms& fms,
	Eigen::ArrayXXd& Ctot) const
{
	int m = lx1 + lx2;
	EriCoefs::AllMBlock Cxl = Cx.allM(lx1, lx2);
	Ctot = Cxl[m] * fms[m];
	for (--m; m >= 0; --m)
		Ctot += Cxl[m] * fms[m];

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
	Ctot = Cxl[lx] * Cyl[ly] * fms[m];
	for (--m; m > 1; --m)
	{
		Cm.setZero();
		int mxmax = std::min(m, lx);
		int mxmin = std::max(m-ly, 0);
		for (int mx = mxmax, my = m-mxmax ; mx >= mxmin; --mx, ++my)
			Cm += Cxl[mx] * Cyl[my];
		Ctot += Cm * fms[m];
	}
	Ctot += (Cxl[1] * Cyl[0] + Cxl[0] * Cyl[1]) * fms[1]
		+ Cxl[0] * Cyl[0] * fms[0];

	return mulWeights(Ctot);
}

double CGTOShellQuad::eri_xx(const AngMomPairIterator& iter,
	const EriCoefs& Cx, const EriCoefs& Cy, const EriCoefs& Cz,
	const Fms& fms, Eigen::ArrayXXd& Cm, Eigen::ArrayXXd& Ctot) const
{
	int lx = iter.lx(), ly = iter.ly(), lz = iter.lz();
	
	if (lx == 0)
		return eri_xx(iter.ly1(), iter.lz1(), iter.ly2(), iter.lz2(), Cy, Cz, fms, Cm, Ctot);
	if (ly == 0)
		return eri_xx(iter.lx1(), iter.lz1(), iter.lx2(), iter.lz2(), Cx, Cz, fms, Cm, Ctot);
	if (lz == 0)
		return eri_xx(iter.lx1(), iter.ly1(), iter.lx2(), iter.ly2(), Cx, Cy, fms, Cm, Ctot);
	
	int m = iter.ltot();
	EriCoefs::AllMBlock Cxl = Cx.allM(iter.lx1(), iter.lx2());
	EriCoefs::AllMBlock Cyl = Cy.allM(iter.ly1(), iter.ly2());
	EriCoefs::AllMBlock Czl = Cz.allM(iter.lz1(), iter.lz2());
	Ctot = fms[m] * Cxl[lx] * Cyl[ly] * Czl[lz]
		+ fms[m-1] * (Cxl[lx-1] * Cyl[ly] * Czl[lz]
			+ Cxl[lx] * (Cyl[ly-1] * Czl[lz] + Cyl[ly] * Czl[lz-1]));
	for (m-=2; m > 1; --m)
	{
		Cm.setZero();
		int mxmax = std::min(m, lx);
		for (int mx = mxmax; mx >= 0; --mx)
		{
			int mymax = std::min(m-mx, ly);
			int mymin = std::max(m-mx-lz, 0);
			for (int my = mymax, mz = m-mx-mymax; my >= mymin; --my, ++mz)
				Cm += Cxl[mx] * Cyl[my] * Czl[mz];
		}
		Ctot += fms[m] * Cm;
	}
	Ctot += fms[1] * (Cxl[1] * Cyl[0] * Czl[0]
			+ Cxl[0] * (Cyl[1] * Czl[0] + Cyl[0] * Czl[1]))
		+ fms[0] * Cxl[0] * Cyl[0] * Czl[0];

	return mulWeights(Ctot);
}

double CGTOShellQuad::eri_10(const EriCoefs::AllMBlock& Cxl,
	const Fms& fms) const
{
	return mulWeights(Cxl[1]*fms[1] + Cxl[0]*fms[0]);
}

double CGTOShellQuad::eri_20(const EriCoefs::AllMBlock& Cxl,
	const Fms& fms) const
{
	return mulWeights(Cxl[2]*fms[2] + Cxl[1]*fms[1] + Cxl[0]*fms[0]);
}

double CGTOShellQuad::eri_11(const EriCoefs::AllMBlock& Cxl,
	const EriCoefs::AllMBlock& Cyl,
	const Fms& fms, Eigen::ArrayXXd& Ctot) const
{
	Ctot = Cxl[1]*Cyl[1]*fms[2]
		+ (Cxl[1]*Cyl[0]+Cxl[0]*Cyl[1])*fms[1]
		+ Cxl[0]*Cyl[0]*fms[0];
	return mulWeights(Ctot);
}

void CGTOShellQuad::setEri() const
{
	Eigen::ArrayXXd T = widthsReduced()
		* (dPQ(0).square() + dPQ(1).square() + dPQ(2).square());
	Fms fms(_lsum, T, KKW());

	eri(0, 0, 0, 0, 0, 0) = mulWeights(fms[0]);
	if (_lsum == 0)
	{
		_have_eri = true;
		return;
	}

	int n1 = _pAB.size(), n2 = _pCD.size();
	EriCoefs Cx(_lAB, _lCD, n1, n2);
	EriCoefs Cy(_lAB, _lCD, n1, n2);
	EriCoefs Cz(_lAB, _lCD, n1, n2);

	elecRepPrim1d_abcd(0, Cx);
	elecRepPrim1d_abcd(1, Cy);
	elecRepPrim1d_abcd(2, Cz);

	if (_lAB >= 1)
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
	if (_lCD >= 1)
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
	if (_lAB >= 2)
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
	if (_lAB >= 1 && _lCD >= 1)
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
	if (_lCD >= 2)
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
		AngMomPairIterator iter(l, _lAB, _lCD);
		for ( ; !iter.end(); ++iter)
			eri(iter) = eri_xx(iter, Cx, Cy, Cz, fms, Cm, Ctot);
	}

	_have_eri = true;
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1,
	int lx2, int ly2, int lz2, int lzD) const
{
	if (lzD == 0)
		return eri(lx1, ly1, lz1, lx2, ly2, lz2);
	return eri(lx1, ly1, lz1, lx2, ly2, lz2, lzD-1)
		- _pCD.dAB(2) * eri(lx1, ly1, lz1, lx2, ly2, lz2-1, lzD-1);
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1,
	int lx2, int ly2, int lz2, int lyD, int lzD) const
{
	if (lyD == 0)
		return eri(lx1, ly1, lz1, lx2, ly2, lz2, lzD);
	return eri(lx1, ly1, lz1, lx2, ly2, lz2, lyD-1, lzD)
		- _pCD.dAB(1) * eri(lx1, ly1, lz1, lx2, ly2-1, lz2, lyD-1, lzD);
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1,
	int lx2, int ly2, int lz2, int lxD, int lyD, int lzD) const
{
	if (_pos_sym == POS_SYM_AACD || _pos_sym == POS_SYM_ABCD)
	{
		if (lxD == 0)
			return eri(lx1, ly1, lz1, lx2, ly2, lz2, lyD, lzD);
		return eri(lx1, ly1, lz1, lx2, ly2, lz2, lxD-1, lyD, lzD)
			- _pCD.dAB(0) * eri(lx1, ly1, lz1, lx2-1, ly2, lz2, lxD-1, lyD, lzD);
	}
	return eri(lx1, ly1, lz1, lx2, ly2, lz2);
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1, int lzB,
	int lx2, int ly2, int lz2, int lxD, int lyD, int lzD) const
{
	if (lzB == 0)
		return eri(lx1, ly1, lz1, lx2, ly2, lz2, lxD, lyD, lzD);
	return eri(lx1, ly1, lz1, lzB-1, lx2, ly2, lz2, lxD, lyD, lzD)
		- _pAB.dAB(2) * eri(lx1, ly1, lz1-1, lzB-1, lx2, ly2, lz2, lxD, lyD, lzD);
}

inline double CGTOShellQuad::eri(int lx1, int ly1, int lz1, int lyB, int lzB,
	int lx2, int ly2, int lz2, int lxD, int lyD, int lzD) const
{
	if (lyB == 0)
		return eri(lx1, ly1, lz1, lzB, lx2, ly2, lz2, lxD, lyD, lzD);
	return eri(lx1, ly1, lz1, lyB-1, lzB, lx2, ly2, lz2, lxD, lyD, lzD)
		- _pAB.dAB(1) * eri(lx1, ly1-1, lz1, lyB-1, lzB, lx2, ly2, lz2, lxD, lyD, lzD);
}

double CGTOShellQuad::eri(int lx1, int ly1, int lz1, int lxB, int lyB, int lzB,
	int lx2, int ly2, int lz2, int lxD, int lyD, int lzD) const
{
	if (_pos_sym == POS_SYM_ABCC || _pos_sym == POS_SYM_ABCD)
	{
		if (lxB == 0)
			return eri(lx1, ly1, lz1, lyB, lzB, lx2, ly2, lz2, lxD, lyD, lzD);
		if (lxB == 1)
			return eri(lx1, ly1, lz1, lyB, lzB, lx2, ly2, lz2, lxD, lyD, lzD)
				- _pAB.dAB(0) * eri(lx1-1, ly1, lz1, lyB, lzB, lx2, ly2, lz2, lxD, lyD, lzD);
		return eri(lx1, ly1, lz1, lxB-1, lyB, lzB, lx2, ly2, lz2, lxD, lyD, lzD)
			- _pAB.dAB(0) * eri(lx1-1, ly1, lz1, lxB-1, lyB, lzB, lx2, ly2, lz2, lxD, lyD, lzD);
	}
	return eri(lx1, ly1, lz1, lx2, ly2, lz2, lxD, lyD, lzD);
}

CGTOShellQuad::PositionSymmetry CGTOShellQuad::symmetry(const CGTOShellPair& pAB,
	const CGTOShellPair& pCD)
{
	if (pAB.samePositionId())
	{
		if (pCD.samePositionId())
			return pAB.positionIdA() == pCD.positionIdA()
				? POS_SYM_AAAA : POS_SYM_AACC;
		return POS_SYM_AACD;
	}
	else if (pCD.samePositionId())
	{
		return POS_SYM_ABCC;
	}
	return POS_SYM_ABCD;
}
