#ifndef FMCOEFS_HH
#define FMCOEFS_HH

#include <Eigen/Core>
#include "EriCoefs.hh"

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
		const Eigen::ArrayBase<Derived1>& C1);
	template <typename Derived0, typename Derived1>
	void multiplyRow(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1);
	template <typename Derived0, typename Derived1>
	void multiply(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1);
	template <typename Derived1>
	void multiply_noC0(const Eigen::ArrayBase<Derived1>& C1);
	template <typename Derived0, typename Derived1, typename Derived2>
	void multiplyCol(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2);
	template <typename Derived0, typename Derived1, typename Derived2>
	void multiplyRow(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2);
	template <typename Derived0, typename Derived1, typename Derived2>
	void multiply(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2);
	template <typename Derived1, typename Derived2>
	void multiply_noC0(const Eigen::ArrayBase<Derived1>& C1,
		const Eigen::ArrayBase<Derived2>& C2);
	void multiply(const EriCoefs& coefs, int l1, int l2);

private:
	int _p_size;
	int _q_size;
	Eigen::ArrayXXd _C;
	int _mmin;
	int _mmax;
};

template <typename Derived0, typename Derived1>
void FmCoefs::multiplyCol(const Eigen::ArrayBase<Derived0>& C0,
	const Eigen::ArrayBase<Derived1>& C1)
{
	block(_mmax+1) = block(_mmax) * C1;
	for (int m = _mmax; m > _mmin; --m)
		block(m) = block(m).colwise()*C0 + block(m-1)*C1;
	block(_mmin).colwise() *= C0;
	++_mmax;
}

template <typename Derived0, typename Derived1>
void FmCoefs::multiplyRow(const Eigen::ArrayBase<Derived0>& C0,
		const Eigen::ArrayBase<Derived1>& C1)
{
	block(_mmax+1) = block(_mmax) * C1;
	for (int m = _mmax; m > _mmin; --m)
		block(m) = block(m).rowwise()*C0 + block(m-1)*C1;
	block(_mmin).rowwise() *= C0;
	++_mmax;
}

template <typename Derived0, typename Derived1>
void FmCoefs::multiply(const Eigen::ArrayBase<Derived0>& C0,
	const Eigen::ArrayBase<Derived1>& C1)
{
	block(_mmax+1) = block(_mmax) * C1;
	for (int m = _mmax; m > _mmin; --m)
		block(m) = block(m)*C0 + block(m-1)*C1;
	block(_mmin) *= C0;
	++_mmax;
}

template <typename Derived1>
void FmCoefs::multiply_noC0(const Eigen::ArrayBase<Derived1>& C1)
{
	block(_mmin, _mmax-_mmin+1) *= C1.replicate(1, _mmax-_mmin+1);
	++_mmin;
	++_mmax;
}

template <typename Derived0, typename Derived1, typename Derived2>
void FmCoefs::multiplyCol(const Eigen::ArrayBase<Derived0>& C0,
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
void FmCoefs::multiplyRow(const Eigen::ArrayBase<Derived0>& C0,
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
void FmCoefs::multiply(const Eigen::ArrayBase<Derived0>& C0,
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
void FmCoefs::multiply_noC0(const Eigen::ArrayBase<Derived1>& C1,
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

#endif // FMCOEFS_HH