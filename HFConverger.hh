#ifndef HFCONVERGER_HH
#define HFCONVERGER_HH

/*!
 * \file HFConverger.hh
 * \brief Definition of the HFConverger class
 */

#include "HartreeFock.hh"

/*!
 * \brief Base class for Hartree-Fock convergers
 *
 * Class HFConverger is the base class for convergence accelerators of a
 * Hartree-Fock SCF calculation. It defines the interface for convergence
 * classes. A HFConverger object can itself also be used as a convergence
 * accelerator, but it will simply return the Fock matrix unchanged.
 */
class HFConverger
{
public:
	//! Typedef for a (smart) pointer to a HFConverger
	typedef std::unique_ptr<HFConverger> Ptr;

	/*!
	 * \brief Constructor
	 *
	 * Create a new HFConverger for the HF calculation \a hf using basis
	 * \a basis.
	 * \param hf    The HartreeFock object in which to accelerate convergence
	 * \param basis The basis set in which the HF calculation is done
	 */
	HFConverger(HartreeFock& hf, const Basis& basis):
		hf(hf), basis(basis) {}

	/*!
	 * \brief Create a new Fock matrix
	 *
	 * Create a new guess for the Fock matrix. On input \a F contains the
	 * current calculated iteration of the Fock matrix, on exit it is
	 * replaced by a hopefull better guess.
	 * \param F    The Fock matrix
	 * \param Etot The total energy of the system
	 */
	virtual void step(Eigen::MatrixXd& /* F */, double /* Etot */) {};
	virtual void step(Eigen::MatrixXd& /* Fa */, Eigen::MatrixXd& /* Fb */,
		double /* Etot */) {};
	//! Return the error in the last step of the calculations
	virtual double error() const { return 0; }

	/*!
	 * \brief Create a new HFConverger
	 *
	 * Pseudo-constructor for creating a HFConverger object implementing
	 * algorithm \a which. Valid values for \a which are \c "none",
	 * \c "diis", and \c "damp".
	 * \param which The type of converger to create.
	 * \param hf    The HartreeFock object in which to accelerate convergence
	 * \param basis The basis set in which the HF calculation is done
	 * \return Smart pointer to a new HFConverger object
	 */
	static Ptr create(const std::string& which, HartreeFock& hf,
		const Basis& basis);

protected:
	//! The Hartree-Fock calculation to speed up
	HartreeFock& hf;
	//! The basis set in which the calculation is done.
	const Basis& basis;
};

/*!
 * \brief Implementation of a HF convergence accelerator
 *
 * Template class HFConvergerImpl implements a HF convergence accelerator using
 * algorithm \a Algo. Specializations of this class can call the wrapped \a Algo
 * object using the necessary parameters.
 * \tparam Algo The type of the actual convergence accelerator used.
 */
template <typename Algo>
class HFConvergerImpl: public HFConverger
{
public:
	//! Constructor
	HFConvergerImpl(HartreeFock& hf, const Basis& basis):
		HFConverger(hf, basis), _algo() {}

	//! \see HFConverger::step()
	void step(Eigen::MatrixXd& F, double Etot);
	//! \see HFConverger::step()
	void step(Eigen::MatrixXd& Fa, Eigen::MatrixXd& Fb, double Etot);
	//! \see HFConverger::error()
	double error() const { return _algo.error(); }

private:
	//! The actual optimizer used
	Algo _algo;
};

#endif // HFCONVERGER_HH