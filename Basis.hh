#ifndef BASIS_HH
#define BASIS_HH

/*!
 * \file Basis.hh
 * \brief Definition of the Basis class
 */

#include <vector>
#include <bitset>
#include <tr1/memory>
#include "AbstractBFPair.hh"

/*!
 * \brief Class representing a basis
 *
 * Class Basis holds a basis for the system being computed. The basis is
 * constructed from a basis set definition, which is applied to a geometry
 * containing the positions of the atoms in the system.
 */
class Basis
{
	public:
		//! Local typedef for a (shared) pointer to a basis function
		typedef std::tr1::shared_ptr<AbstractBF> BasisFunPtr;
		//! Local typedef for a list of basis functions
		typedef std::vector<BasisFunPtr> BasisFunList;
		//! Local typedef for a (shared) pointer to a pair of basis functions
		typedef std::tr1::shared_ptr<AbstractBFPair> PairPtr;
		//! Local typedef for a list of basis function pairs
		typedef std::vector<PairPtr> PairList;

		//! Enumeration for the different status flags
		enum StatusFlag
		{
			//! Set when the list of function pairs is current
			PAIRS_CURRENT,
			//! Set when the overlap matrix is computed and up to date
			OVERLAP_CURRENT,
			//! Set when the kinetic energy matrix is computed and up to date
			KINETIC_CURRENT,
			//! Set when the nuclear attraction matrix is computed and up to date
			NUC_ATTR_CURRENT,
			//! The number of status flags
			NR_FLAGS
		};

		/*!
		 * \brief Constructor
		 *
		 * Create a new and empty basis
		 */
		Basis(): _funs(), _status(), _overlap(), _kinetic() {}

		//! Return the number of functions in this basis
		int size() const { return _funs.size(); }

		/*!
		 * \brief Add a basis function
		 *
		 * Add basis function \a bf to this basis.
		 */
		void add(const BasisFunPtr& ptr)
		{
			_funs.push_back(ptr);
			_status.reset();
		}

		/*!
		 * \brief Return the overlap matrix
		 *
		 * Return the matrix with the overlaps between the fucntions
		 * in this basis, computing it first if necessary.
		 */
		const Eigen::MatrixXd& overlap() const
		{
			if (!_status.test(OVERLAP_CURRENT))
				calcOneElectron();
			return _overlap;
		}
		/*!
		 * \brief Return the kinetic energy matrix
		 *
		 * Return the matrix with the kinetic energy integrals
		 * for the electrons in this basis, computing it first if
		 * necessary.
		 */
		const Eigen::MatrixXd& kineticEnergy() const
		{
			if (!_status.test(KINETIC_CURRENT))
				calcOneElectron();
			return _kinetic;
		}
		/*!
		 * \brief Return the nuclear attraction matrix
		 *
		 * Return the matrix with the muclear attraction integrals
		 * for the electrons in this basis, computing it first if
		 * necessary.
		 * \param nuc_pos Positions of the nuclei
		 * \param nuc_charge Charges of the nuclei
		 */
		const Eigen::MatrixXd& nuclearAttraction(
			const Eigen::MatrixXd& nuc_pos,
			const Eigen::VectorXd& nuc_charge) const
		{
			if (!_status.test(NUC_ATTR_CURRENT))
				calcNuclearAttraction(nuc_pos, nuc_charge);
			return _nuc_attr;
		}

		/*!
		 * \brief Print this basis
		 *
		 * Print a textual representation of this basis on output
		 * stream \a os.
		 * \param os The output stream to write to
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;

	private:
		//! The list of basis functions
		BasisFunList _funs;
		//! The list of function pairs
		mutable PairList _pairs;
		//! Status flags for the basis
		mutable std::bitset<NR_FLAGS> _status;
		//! The overlap matrix for this basis
		mutable Eigen::MatrixXd _overlap;
		//! The kinetic energy matrix for this basis
		mutable Eigen::MatrixXd _kinetic;
		//! The nuclear attraction matrix for this basis
		mutable Eigen::MatrixXd _nuc_attr;

		//! Create the list of basis function pairs
		void setPairs() const;
		//! Compute the overlap matrix in this basis
		void calcOverlap() const;
		//! Compute the kinetic energy matrix in this basis
		void calcKinetic() const;
		//! Compute the one-electron matrices in this basis
		void calcOneElectron() const;
		/*!
		 * \brief Compute the nuclear attraction
		 * 
		 * Compute the nuclear integrals due to the nuclei on positions
		 * \a nuc_pos with charges \a nuc_charge.
		 * \param nuc_pos Positions of the nuclei
		 * \param nuc_charge Charges of the nuclei
		 */
		void calcNuclearAttraction(const Eigen::MatrixXd& nuc_pos,
			const Eigen::VectorXd& nuc_charge) const;
};

namespace {

/*!
 * \brief Print a basis
 *
 * Print a textual representation of basis \a basis on output stream \a os.
 * \param os    The output stream to write to
 * \param basis The basis to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const Basis& basis)
{
	return basis.print(os);
}

} // namespace

#endif // BASIS_HH
