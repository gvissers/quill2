#include "Damper.hh"

void Damper::step(Eigen::MatrixXd& F, double factor, double Etot)
{
	if (!_started)
	{
		_last = F;
		_error = std::abs(Etot);
		_started = true;
	}
	else
	{
		_error = factor * (_last - F).lpNorm<Eigen::Infinity>();
		F = _last = (1-factor)*F + factor*_last;
	}
}

void Damper::step(Eigen::MatrixXd& Fa, Eigen::MatrixXd& Fb, double factor, double Etot)
{
	int size = Fa.rows();
	if (!_started)
	{
		_last = Eigen::MatrixXd::Zero(size, 2*size);
		_last.leftCols(size) = Fa;
		_last.rightCols(size) = Fb;
		_error = std::abs(Etot);
		_started = true;
	}
	else
	{
		_last.leftCols(size) = (1-factor)*Fa + factor*_last.leftCols(size);
		_last.rightCols(size) = (1-factor)*Fb + factor*_last.rightCols(size);
		_error = std::max((_last.leftCols(size) - Fa).lpNorm<Eigen::Infinity>(),
			(_last.rightCols(size) - Fb).lpNorm<Eigen::Infinity>());
	}
}
