#include <iomanip>
#include <Eigen/Dense>
#include <Exception.hh>
#include "HartreeFock.hh"
#include "DIIS.hh"
#include "exceptions.hh"

const int HartreeFock::default_max_iter = 50;
const double HartreeFock::default_tolerance = 1.e-10;

HartreeFock::HartreeFock(): _max_iter(default_max_iter),
	_tolerance(default_tolerance), _status(),
	_orbitals(), _orb_ener(), _density() {}

HartreeFock::HartreeFock(const Basis& basis, const Geometry& geometry,
	int multiplicity): _max_iter(default_max_iter),
	_tolerance(default_tolerance), _status(),
	_orbitals(), _orb_ener(), _density()
{
	iterate(basis, geometry, multiplicity);
}

void HartreeFock::iterate(const Basis& basis, const Geometry& geometry,
	int multiplicity)
{
	_status.reset();
	setMultiplicity(geometry, multiplicity);

	double nuc_rep = geometry.nuclearRepulsion();
	std::cout << "nuclear repulsion = " << nuc_rep << "\n";

	Eigen::MatrixXd H = basis.kineticEnergy()
		+ basis.nuclearAttraction(geometry.positions(), geometry.charges());
	calcOrbitals(H, basis.overlap());

	_energy = nuc_rep;
	if (_nr_double > 0) _energy += 2 * _orb_ener.head(_nr_double).sum();
	if (_nr_single > 0) _energy += _orb_ener.segment(_nr_double, _nr_single).sum();
	std::cout << "Energy: " << _energy << "\n";
	std::cout << "orbital energies: " << _orb_ener.transpose() << "\n";

	Eigen::MatrixXd G;
	DIIS diis(basis.size());
	for (int iter = 0; iter < _max_iter; ++iter)
	{
		double last_energy = _energy;
		
		const Eigen::MatrixXd& D = density();
		basis.electronRepulsion(D, G);
		_energy = (H + 0.5*G).cwiseProduct(D).sum() + nuc_rep;
 		std::cout << "Energy: " << std::setprecision(15) << _energy << "\n";
		//std::cout << "orbital energies: " << _orb_ener.transpose() << "\n";

		if (std::abs(_energy - last_energy) < _tolerance)
		{
			_status.set(ENERGY_CONVERGED);
			return;
		}

		Eigen::MatrixXd F = G + H;
		diis.step(F, D, basis.overlap(), _energy);
		calcOrbitals(F, basis.overlap());
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
}

void HartreeFock::calcDensity()
{
	if (!_status.test(ORBITALS_CURRENT))
		throw Li::Exception("No orbitals to compute density matrix with");

	int nr_func = _orbitals.cols();
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