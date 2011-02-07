#ifndef BASIS_HH
#define BASIS_HH

/*!
 * \file Basis.hh
 * \brief Definition of the Basis class
 */

#include <vector>
#include <tr1/memory>
#include "AbstractBF.hh"

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
		typedef std::tr1::shared_ptr<AbstractBF> BasisFunPtr;
		//! Local typedef for a list of basis functions
		typedef std::vector<BasisFunPtr> BasisFunList;

		/*!
		 * \brief Constructor
		 *
		 * Create a new and empty basis
		 */
		Basis(): _funs() {}

		//! Return the number of functions in this basis
		int size() const { return _funs.size(); }

		/*!
		 * \brief Add a basis function
		 *
		 * Add basis function \a bf to this basis.
		 */
		void add(const BasisFunPtr& ptr) { _funs.push_back(ptr); }

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
