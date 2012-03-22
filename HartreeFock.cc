#include <iomanip>
#include <Eigen/Dense>
#include <Exception.hh>
#include "HartreeFock.hh"
#include "HFConverger.hh"
#include "exceptions.hh"

const int HartreeFock::default_max_iter = 100;
const double HartreeFock::default_tolerance = 1.e-10;

HartreeFock::HartreeFock(): _max_iter(default_max_iter),
	_tolerance(default_tolerance), _conv_method("diis"),  _status(),
	_orbitals(), _orb_ener(), _density() {}

HartreeFock::HartreeFock(const Basis& basis, const Geometry& geometry,
	int multiplicity, bool restricted): _max_iter(default_max_iter),
	_tolerance(default_tolerance), _conv_method("diis"), _status(),
	_orbitals(), _orb_ener(), _density()
{
	iterate(basis, geometry, multiplicity, restricted);
}

void HartreeFock::iterate(const Basis& basis, const Geometry& geometry,
	int multiplicity, bool restricted)
{
	_status.reset();
	_restricted = restricted;
	setMultiplicity(geometry, multiplicity);

	double nuc_rep = geometry.nuclearRepulsion();
	std::cout << "nuclear repulsion = " << nuc_rep << "\n";

	Eigen::MatrixXd H = basis.kineticEnergy()
		+ basis.nuclearAttraction(geometry.positions(), geometry.charges());
	const Eigen::MatrixXd& S = basis.overlap();

	if (restricted)
		iterateRestricted(basis, H, S, nuc_rep);
	else
		iterateUnrestricted(basis, H, S, nuc_rep);
}

void HartreeFock::iterateRestricted(const Basis& basis,
	const Eigen::MatrixXd& H, const Eigen::MatrixXd& S,
	double nuc_rep)
{
	calcOrbitals(H, S);

	_energy = nuc_rep;
	if (_nr_double > 0) _energy += 2 * _orb_ener.head(_nr_double).sum();
	if (_nr_single > 0) _energy += _orb_ener.segment(_nr_double, _nr_single).sum();
	std::cout << "Energy: " << _energy << "\n";

	Eigen::MatrixXd J, K;
	HFConverger::Ptr converger = HFConverger::create(_conv_method, *this,
		basis);
	for (int iter = 0; iter < _max_iter; ++iter)
	{
		double last_energy = _energy;
		
		const Eigen::MatrixXd& D = density();
		basis.twoElectron(D, J, K);
		_energy = (H + 0.5*(J - 0.5*K)).cwiseProduct(D).sum() + nuc_rep;
		std::cout << std::setw(3) << iter << " "
			<< std::setw(20) << std::setprecision(15) << std::fixed << _energy << " "
			<< std::setw(13) << std::setprecision(6) << std::scientific << converger->error()
			<< "\n";
		//std::cout << "orbital energies: " << _orb_ener.transpose() << "\n";

		if (std::abs(_energy - last_energy) < _tolerance)
		{
			_status.set(ENERGY_CONVERGED);
			return;
		}

		Eigen::MatrixXd F = H + J - 0.5*K;
		converger->step(F, _energy);
		calcOrbitals(F, S);
	}
	
	throw NoConvergence();
}

void HartreeFock::iterateUnrestricted(const Basis& basis,
	const Eigen::MatrixXd& H, const Eigen::MatrixXd& S,
	double nuc_rep)
{
	int nr_func = basis.size();

	calcOrbitals(H, H, S);	
	_energy = nuc_rep + _orb_ener.head(_nr_alpha).sum()
		+ _orb_ener.segment(nr_func, _nr_beta).sum();
	std::cout << "Energy: " << _energy << "\n";

	Eigen::MatrixXd Ja, Ka, Jb, Kb;
	HFConverger::Ptr converger = HFConverger::create(_conv_method, *this,
		basis);
	for (int iter = 0; iter < _max_iter; ++iter)
	{
		double last_energy = _energy;

		Eigen::MatrixXd::ConstColsBlockXpr Da = density().leftCols(nr_func);
		Eigen::MatrixXd::ConstColsBlockXpr Db = density().rightCols(nr_func);
		basis.twoElectron(Da, Ja, Ka);
		basis.twoElectron(Db, Jb, Kb);

		double energyA = (H + 0.5*(Ja+Jb-Ka)).cwiseProduct(Da).sum();
		double energyB = (H + 0.5*(Ja+Jb-Kb)).cwiseProduct(Db).sum();
		_energy = energyA + energyB + nuc_rep;

		std::cout << std::setw(3) << iter << " "
			<< std::setw(20) << std::setprecision(15) << std::fixed << _energy << " "
			<< std::setw(13) << std::setprecision(6) << std::scientific << converger->error() << "\n";
		//std::cout << "orbital energies: " << _orb_ener.transpose() << "\n";

		if (std::abs(_energy - last_energy) < _tolerance)
		{
			_status.set(ENERGY_CONVERGED);
			return;
		}

		Eigen::MatrixXd Fa = H + Ja + Jb - Ka;
		Eigen::MatrixXd Fb = H + Ja + Jb - Kb;
		converger->step(Fa, Fb, _energy);
		calcOrbitals(Fa, Fb, S);
	}

	throw NoConvergence();
}

double HartreeFock::energy()
{
	return _energy;
}

void HartreeFock::setMultiplicity(const Geometry& geometry, int multiplicity)
{
	int nr_elec = geometry.totalCharge();

	if (multiplicity < 1)
		// Use lowest possible spin multiplicity
		multiplicity = 1 + nr_elec % 2;
	else if (nr_elec % 2 == multiplicity % 2 || multiplicity > nr_elec+1)
		throw Li::Exception("Multiplicity does not match number of electrons");

	_nr_single = multiplicity - 1;
	_nr_double = (nr_elec - _nr_single) / 2;
	if (!_restricted)
	{
		_nr_alpha = _nr_double + _nr_single;
		_nr_beta = _nr_double;
	}
}

void HartreeFock::calcDensity()
{
	if (!_status.test(ORBITALS_CURRENT))
		throw Li::Exception("No orbitals to compute density matrix with");

	int nr_func = _orbitals.rows();
	if (_restricted)
	{
		_density.resize(nr_func, nr_func);
		for (int j = 0; j < nr_func; j++)
		{
			for (int i = 0; i <= j; i++)
			{
				double sum = 0.0;
				if (_nr_double)
					sum += 2 * _orbitals.row(i).head(_nr_double).dot(
						_orbitals.row(j).head(_nr_double));
				if (_nr_single)
					sum += _orbitals.row(i).segment(_nr_double, _nr_single).dot(
						_orbitals.row(j).segment(_nr_double, _nr_single));
				_density(j,i) = _density(i,j) = sum;
			}
		}
        }
        else
	{
		_density.resize(nr_func, 2*nr_func);
		if (_nr_alpha == 0)
		{
			_density.leftCols(nr_func).setZero();
		}
		else
		{
			Eigen::MatrixXd::ColsBlockXpr occ
				= _orbitals.leftCols(_nr_alpha);
			_density.leftCols(nr_func) = occ * occ.transpose();
		}
		if (_nr_beta == 0)
		{
			_density.rightCols(nr_func).setZero();
		}
		else
		{
			Eigen::Block<Eigen::MatrixXd> occ
				= _orbitals.block(0, nr_func, nr_func, _nr_beta);
			_density.rightCols(nr_func) = occ * occ.transpose();
		}
	}

	_status.set(DENSITY_CURRENT);
}

void HartreeFock::calcOrbitals(const Eigen::MatrixXd& F,
	const Eigen::MatrixXd& S)
{
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(F, S);
	_orbitals = solver.eigenvectors();
	_orb_ener = solver.eigenvalues();
	_status.set(ORBITALS_CURRENT);
	_status.reset(DENSITY_CURRENT);
}

void HartreeFock::calcOrbitals(const Eigen::MatrixXd& Fa,
	const Eigen::MatrixXd& Fb, const Eigen::MatrixXd& S)
{
	int nr_func = S.rows();
	_orbitals.resize(nr_func, 2*nr_func);
	_orb_ener.resize(2*nr_func);

	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(Fa, S);
	_orbitals.leftCols(nr_func) = solver.eigenvectors();
	_orb_ener.head(nr_func) = solver.eigenvalues();
	solver.compute(Fb, S);
	_orbitals.rightCols(nr_func) = solver.eigenvectors();
	_orb_ener.tail(nr_func) = solver.eigenvalues();

	_status.set(ORBITALS_CURRENT);
	_status.reset(DENSITY_CURRENT);
}