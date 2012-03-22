#include "HFConverger.hh"
#include "DIIS.hh"
#include "Damper.hh"
#include "support.hh"

HFConverger::Ptr HFConverger::create(const std::string& which,
	HartreeFock& hf, const Basis& basis)
{
	std::string lw = lower(which);
	if (which == "none")
		return Ptr(new HFConverger(hf, basis));
	else if (which == "damp")
		return Ptr(new HFConvergerImpl<Damper>(hf, basis));
	else if (which == "diis")
		return Ptr(new HFConvergerImpl<DIIS>(hf, basis));
	else
		throw Li::Exception("Unknown convergence method \"" + which + "\"");
}

template <>
void HFConvergerImpl<Damper>::step(Eigen::MatrixXd& F, double Etot)
{
	_algo.step(F, 0.25, Etot);
}

template <>
void HFConvergerImpl<Damper>::step(Eigen::MatrixXd& Fa, Eigen::MatrixXd& Fb,
	double Etot)
{
	_algo.step(Fa, Fb, 0.25, Etot);
}

template <>
void HFConvergerImpl<DIIS>::step(Eigen::MatrixXd& F, double Etot)
{
	_algo.step(F, hf.density(), basis.overlap(), basis.ortho(), Etot);
}

template <>
void HFConvergerImpl<DIIS>::step(Eigen::MatrixXd& Fa, Eigen::MatrixXd& Fb,
	double Etot)
{
	int nr_func = basis.size();
	Eigen::MatrixXd::ConstColsBlockXpr Da = hf.density().leftCols(nr_func);
	Eigen::MatrixXd::ConstColsBlockXpr Db = hf.density().rightCols(nr_func);
	_algo.step(Fa, Da, Fb, Db, basis.overlap(), basis.ortho(), Etot);
}
