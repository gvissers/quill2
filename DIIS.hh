#ifndef DIIS_HH
#define DIIS_HH

/*!
 * \file DIIS.hh
 * \brief Definition of the DIIS class
 */

#include <Eigen/Core>

/*!
 * \brief The DIIS convergence optimizer
 *
 * Class DIIS implements the Direct Inversion in the Iterative Subspace
 * extrapolation technique to speed up the convergence in an SCF calculation.
 * It tries to make a better guess for the Fock matrix in the next iteration
 * by creating a linear combination of Fock matrices from previous iterations,
 * where the error term \f$FDS-SDF\f$ (with \f$F\f$ the Fock matrix, \f$D\f$
 * the density matrix, and \f$S\f$ the overlap matrix) is minimized.
 */
class DIIS
{
public:
	/*!
	 * \brief Constructor
	 *
	 * Create a new DIIS object, for a basis of size \a size.
	 * \param size The size of the matrices (basis set size).
	 */
	DIIS(int size): _size(size), _err_used(0),
		_errors((_size-1)*size/2, 10), _values(),
		_started(false) {}

	/*!
	 * \brief Create a new Fock matrix
	 *
	 * Create the new Fock matrix from the current guess for the matrix and
	 * the current density matrix. If the DIIS algorithm has not been
	 * started yet (because the error term is still too large with respect
	 * to the computed energy), it simply returns \a F. Otherwise a linear
	 * combination of previous Fock matrices \f$\tilde{F}\f$ that minimizes
	 * the error term \f$\tilde{F}DS-SD\tilde{F}\f$ is computed, and stored
	 * in \a F.
	 * \param F    The current Fock matrix
	 * \param P    The density matrix from which \a F was computed
	 * \param S    The overlap matrix for the basis
	 * \param Etot The total energy
	 */
	void step(Eigen::MatrixXd& F, const Eigen::MatrixXd& D,
		const Eigen::MatrixXd& S, double Etot);

private:
	//! The basis set size
	int _size;
	//! The number of error vectors currently stored
	int _err_used;
	//! The error vectors themselves, one per column
	Eigen::MatrixXd _errors;
	//! The previous Fock matrices
	std::vector<Eigen::MatrixXd> _values;
	//! Whether DIIS was started
	bool _started;
};

#endif // DIIS_HH