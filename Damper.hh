#ifndef DAMPER_HH
#define DAMPER_HH

/*!
 * \file Damper.hh
 * \brief Definition of the Damper class
 */

#include <Eigen/Core>

/*!
 * \brief Class for damping an iteration
 *
 * Class Damper can be used to damp a matrix in an iterative process, such as
 * the Fock matrix in an SCF calculation, to prevent oscillations. Instead of
 * directly using the \f$n\$'th matrix \f$F_n\f$ in the iteration,
 * \f$cF_{n-1} + (1-c)F_{n}\f$ is used, where the damping parameter \f$c\f$ is
 * a constant between 0 and 1 that determineshow strongly the iteration is
 * damped.
 * \sa Class HFConverger
 */
class Damper
{
public:
	//! Constructor
	Damper(): _last(), _started(false) {}

	//! Return the error in the last step of the iteration
	double error() const { return _error; }

	/*!
	 * \brief Create a new matrix
	 *
	 * Create a new matrix from the current and previous iterations of the
	 * matrix. On the first call, this simply returns \a F unchanged, with
	 * an error of \a Etot. On subsequent calls, \f$cF_{-1} + (1-c)F\f$ is
	 * returned, where \f$F_{-1}\f$ is the result of the previous iteration.
	 * \param F      The current matrix matrix
	 * \param factor The damping factor \f$c\f$
	 * \param Etot   Scaling for the error
	 */
	void step(Eigen::MatrixXd& F, double factor, double Etot);
	void step(Eigen::MatrixXd& Fa, Eigen::MatrixXd& Fb, double factor, double Etot);

private:
	//! The previous matrix in the iteration
	Eigen::MatrixXd _last;
	//! The error in the iteration, \f$||F_{n} - F_{n-1}||_\infty\f$.
	double _error;
	//! Whether the iteration has been started
	bool _started;
};

#endif // DAMPER_HH