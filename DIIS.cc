#include <iostream>
#include <Eigen/Dense>
#include <Eigen/MatrixFunctions>
#include "DIIS.hh"

void DIIS::step(Eigen::MatrixXd& F, const Eigen::MatrixXd& D,
	const Eigen::MatrixXd& S, const Eigen::MatrixXd& X, double Etot)
{
	if (_err_vecs_used >= _err_vecs.cols())
	{
		if (_err_vecs_used == 0)
		{
			_size = F.rows();
			_err_vecs.resize((_size-1)*_size/2, 10);
		}
		else
		{
			Eigen::MatrixXd new_err_vecs(_err_vecs.rows(), 2*_err_vecs.cols());
			new_err_vecs.leftCols(_err_vecs_used) = _err_vecs;
			_err_vecs.swap(new_err_vecs);
		}
	}
	
        // The error matrix is X(FDS - SDF)X = XFDX^-1 - h.c, which since S, P,
        // F and X are all symmetric, is antisymmetric. Therefore we only store
        // the upper triangle in vector err.
	Eigen::MatrixXd FDS = X*F*D*S*X;
	int idx = 0;
	for (int i = 1; i < _size; ++i)
		for (int j = 0; j < i; ++j, ++idx)
			_err_vecs(idx, _err_vecs_used) = FDS(i,j) - FDS(j,i);

	_max_err = _err_vecs.col(_err_vecs_used).lpNorm<Eigen::Infinity>();
	if (_max_err > std::abs(0.5*Etot))
		// This Fock matrix will not have a meaningful contribution
		// to convergence. Ignore it.
		return;

	++_err_vecs_used;
	_values.push_back(F);
	if (!_started)
	{
		if (_max_err >= std::abs(0.1*Etot))
			return;
		std::cout << "Starting DIIS\n";
		_started = true;
	}

	auto e = _err_vecs.leftCols(_err_vecs_used);
	Eigen::MatrixXd B(_err_vecs_used+1, _err_vecs_used+1);
	B.topLeftCorner(_err_vecs_used, _err_vecs_used) = e.transpose() * e * 2;
	B.rightCols(1).fill(-1);
	B.bottomRows(1).fill(-1);

	Eigen::VectorXd y = Eigen::VectorXd::Zero(_err_vecs_used+1);
	y[_err_vecs_used] = -1;

	Eigen::VectorXd coefs = B.colPivHouseholderQr().solve(y);
	F *= coefs[_err_vecs_used-1];
	for (int i = 0; i < _err_vecs_used-1; ++i)
		F += coefs[i] * _values[i];
}

void DIIS::step(Eigen::MatrixXd& Fa, const Eigen::MatrixXd& Da,
	Eigen::MatrixXd& Fb, const Eigen::MatrixXd& Db,
	const Eigen::MatrixXd& S, const Eigen::MatrixXd& X, double Etot)
{
	if (_err_vecs_used >= _err_vecs.cols())
	{
		if (_err_vecs_used == 0)
		{
			_size = Fa.rows();
			_err_vecs.resize((_size-1)*_size, 10);
		}
		else
		{
			Eigen::MatrixXd new_err_vecs(_err_vecs.rows(), 2*_err_vecs.cols());
			new_err_vecs.leftCols(_err_vecs_used) = _err_vecs;
			_err_vecs.swap(new_err_vecs);
		}
	}

	Eigen::MatrixXd FDS = X*Fa*Da*S*X;
	int idx = 0;
	for (int i = 1; i < _size; ++i)
		for (int j = 0; j < i; ++j, ++idx)
			_err_vecs(idx, _err_vecs_used) = FDS(i,j) - FDS(j,i);
	FDS = X*Fb*Db*S*X;
	for (int i = 1; i < _size; ++i)
		for (int j = 0; j < i; ++j, ++idx)
			_err_vecs(idx, _err_vecs_used) = FDS(i,j) - FDS(j,i);

	_max_err = _err_vecs.col(_err_vecs_used).lpNorm<Eigen::Infinity>();
	if (_max_err > std::abs(0.5*Etot))
		// This Fock matrix will not have a meaningful contribution
		// to convergence. Ignore it.
		return;

	++_err_vecs_used;
	Eigen::MatrixXd F(_size, 2*_size);
	F.leftCols(_size) = Fa;
	F.rightCols(_size) = Fb;
	_values.push_back(F);
	if (!_started)
	{
		if (_max_err >= std::abs(0.1*Etot))
			return;
		std::cout << "Starting DIIS\n";
		_started = true;
	}

	auto e = _err_vecs.leftCols(_err_vecs_used);
	Eigen::MatrixXd B(_err_vecs_used+1, _err_vecs_used+1);
	B.topLeftCorner(_err_vecs_used, _err_vecs_used) = e.transpose() * e * 2;
	B.rightCols(1).fill(-1);
	B.bottomRows(1).fill(-1);

	Eigen::VectorXd y = Eigen::VectorXd::Zero(_err_vecs_used+1);
	y[_err_vecs_used] = -1;

	Eigen::VectorXd coefs = B.colPivHouseholderQr().solve(y);
	F *= coefs[_err_vecs_used-1];
	for (int i = 0; i < _err_vecs_used-1; ++i)
		F += coefs[i] * _values[i];

	Fa = F.leftCols(_size);
	Fb = F.rightCols(_size);
}