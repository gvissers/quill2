#include <iostream>
#include <Eigen/Dense>
#include "DIIS.hh"

void DIIS::step(Eigen::MatrixXd& F, const Eigen::MatrixXd& D,
	const Eigen::MatrixXd& S, double Etot)
{
	if (_err_used >= _errors.cols())
	{
		Eigen::MatrixXd new_errors(_errors.rows(), 2*_errors.cols());
		new_errors.leftCols(_err_used) = _errors;
		_errors.swap(new_errors);
	}
	
        // The error matrix is FDS - SDF, which since S, P and F are
        // all symmetric, is antisymmetric. Therefore we only store the
        // upper triangle in vector err.
	Eigen::MatrixXd FDS = F*D*S;
	int idx = 0;
	for (int i = 1; i < _size; ++i)
		for (int j = 0; j < i; ++j, ++idx)
			_errors(idx, _err_used) = FDS(i,j) - FDS(j,i);

	double maxerr = _errors.col(_err_used).lpNorm<Eigen::Infinity>();
	if (maxerr > std::abs(0.5*Etot))
		// This Fock matrix will not have a meaningful contribution
		// to convergence. Ignore it.
		return;

	++_err_used;
	_values.push_back(F);
	if (!_started)
	{
		if (maxerr >= std::abs(0.1*Etot))
			return;
		std::cout << "Starting DIIS\n";
		_started = true;
	}

	auto e = _errors.leftCols(_err_used);
	Eigen::MatrixXd B(_err_used+1, _err_used+1);
	B.topLeftCorner(_err_used, _err_used) = e.transpose() * e * 2;
	B.rightCols(1).fill(-1);
	B.bottomRows(1).fill(-1);

	Eigen::VectorXd y = Eigen::VectorXd::Zero(_err_used+1);
	y[_err_used] = -1;

	Eigen::VectorXd coefs = B.colPivHouseholderQr().solve(y);
	F *= coefs[_err_used-1];
	for (int i = 0; i < _err_used-1; ++i)
		F += coefs[i] * _values[i];
}