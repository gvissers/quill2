#ifndef HARTREEFOCK_HH
#define HARTREEFOCK_HH

/*!
 * \file HartreeFock.hh
 * \brief Definition of the HartreeFock class
 */

#include "Basis.hh"
#include "Geometry.hh"

/*!
 * \brief Class for Hartree-Fock calculations
 *
 * Class HartreeFock can be used to perform spin-restricted self-consistent
 * field calculations using the Hartree-Fock method.
 */
class HartreeFock
{
public:
	//! The default maximum number of SCF iterations
	static const int default_max_iter;
	//! The SCF convergence limit for the total energy
	static const double default_tolerance;

	//! Enumeration for the different status flags
	enum StatusFlag
	{
		//! Molecular orbitals have been calculated from the Fock matrix
		ORBITALS_CURRENT,
		//! The density matrix is computed from the current orbitals
		DENSITY_CURRENT,
		//! The total energy of the system has converged
		ENERGY_CONVERGED,
		//! The total number of status flags
		NR_FLAGS
	};

	/*!
	 * \brief Constructor
	 *
	 * Create a new Hartree-Fock iterator. The HF calculation itself can
	 * be started by calling iterate().
	 */
	HartreeFock();
	/*!
	 * \brief Constructor
	 *
	 * Create a new Hartree-Fock iterator, and start a Hartree-Fock
	 * calculation immediately.
	 * \param basis        The electronic basis in which to compute
	 * \param geometry     The geometry of the system
	 * \param multiplicity The spin multiplicity of the state to compute
	 * \param restricted   Whether to do a spin-restricted calculation
	 */
	HartreeFock(const Basis& basis, const Geometry& geometry,
		int multiplicity, bool restricted=true);

	/*!
	 * \brief Start a Hartree-Fock calculation
	 *
	 * Start a Hartree-Fock self-consistent field calculation for the system
	 * described by \a geometry in basis \a basis.
	 * \param basis        The electronic basis in which to compute
	 * \param geometry     The geometry of the system
	 * \param multiplicity The spin multiplicity of the state to compute
	 * \param restricted   Whether to do a spin-restricted calculation
	 */
	void iterate(const Basis& basis, const Geometry& geometry,
		int multiplicity, bool restricted=true);

	//! Return the total energy of the system
	double energy();

	//! Return the density matrix of the system, calculating it if necessary
	const Eigen::MatrixXd& density()
	{
		if (!_status.test(DENSITY_CURRENT))
			calcDensity();
		return _density;
	}

	//! Set the maximum number of SCF iterations allowed
	void setMaxIterations(int max_iter) { _max_iter = max_iter; }
	//! Set the convergence limit in total energy
	void setTolerance(double tolerance) { _tolerance = tolerance; }

private:
	//! The maximum number of SCF iterations
	int _max_iter;
	//! The convergence limit in total energy
	double _tolerance;

	//! The state of the calculation
	std::bitset<NR_FLAGS> _status;
	//! The number of singly occupied orbitals in a restricted calculation
	int _nr_single;
	//! The number of doubly occupied orbitals in a restricted calculation
	int _nr_double;
	//! The number of orbitals with \f$\alpha\f$ spin in an unrestricted calculation
	int _nr_alpha;
	//! The number of orbitals with \f$\beta\f$ spin in an unrestricted calculation
	int _nr_beta;
	//! Whether to do a spin-restricted or spin-unrestriced calculation
	bool _restricted;
	//! The molecular orbitals
	Eigen::MatrixXd _orbitals;
	//! The orbital energies
	Eigen::VectorXd _orb_ener;
	//! The density matrix calculated from the orbitals
	Eigen::MatrixXd _density;
	//! The total energy of the system
	double _energy;

	/*!
	 * \brief Do a spin-restricted calculation
	 *
	 * Do a spin-restricted self-consistent field calculation.
	 * \param basis   The electronic basis in which to compute
	 * \param H       The one-electron Hamiltonian
	 * \param S       The overlap matrix
	 * \param nuc_rep The nuclear repulsion energy of the system
	 */
	void iterateRestricted(const Basis& basis, const Eigen::MatrixXd& H,
		const Eigen::MatrixXd& S, double nuc_rep);
	/*!
	 * \brief Do a spin-unrestricted calculation
	 *
	 * Do a spin-unrestricted self-consistent field calculation.
	 * \param basis   The electronic basis in which to compute
	 * \param H       The one-electron Hamiltonian
	 * \param S       The overlap matrix
	 * \param nuc_rep The nuclear repulsion energy of the system
	 */
	void iterateUnrestricted(const Basis& basis, const Eigen::MatrixXd& H,
		const Eigen::MatrixXd& S, double nuc_rep);

	/*!
	 * \brief Set the spin multiplicity
	 *
	 * Set the spin multiplicity for the system described by \a geometry to
	 * \a multiplicity, and compute the number of doubly and singly occupied
	 * orbitals.
	 * \param geometry The geometry of the system
	 * \param multiplicity The spin multiplicity \f$2S+1\f$ of the system
	 */
	void setMultiplicity(const Geometry& geometry, int multiplicity);
	/*!
	 * \brief Compute the density matrix
	 *
	 * Compute the density matrix from the current orbitals, and the
	 * occupation of the orbitals.
	 */
	void calcDensity();
	/*!
	 * \brief Compute the orbitals
	 *
	 * Compute the electronic orbitals from the Fock matrix \a F and the
	 * overlap matrix \a S.
	 * \param F The current Fock matrix
	 * \param S The overlap matrix of the basis
	 */
	void calcOrbitals(const Eigen::MatrixXd& F, const Eigen::MatrixXd& S);
	/*!
	 * \brief Compute the orbitals
	 *
	 * Compute the electronic orbitals in a spin-unrestricted calculation
	 * from the Fock matrices \a Fa and \a Fb, and the overlap matrix \a S.
	 * \param Fa The current Fock matrix of the \f$\alpha\f$-spin electrons
	 * \param Fb The current Fock matrix of the \f$\beta\f$-spin electrons
	 * \param S  The overlap matrix of the basis
	 */
	void calcOrbitals(const Eigen::MatrixXd& Fa, const Eigen::MatrixXd& Fb,
		const Eigen::MatrixXd& S);
};

#endif // HARTREEFOCK_HH