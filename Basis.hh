#ifndef BASIS_HH
#define BASIS_HH

/*!
 * \file Basis.hh
 * \brief Definition of the Basis class
 */

#include <vector>
#include <bitset>
#include <memory>
#include "AbstractBFQuad.hh"
#include "BFQuadPool.hh"

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
	typedef std::unique_ptr<AbstractBF> BasisFunPtr;
	//! Local typedef for a list of basis functions
	typedef std::vector<BasisFunPtr> BasisFunList;
	//! Local typedef for a (shared) pointer to a pair of basis functions
	typedef std::unique_ptr<AbstractBFPair> PairPtr;
	//! Local typedef for a list of basis function pairs
	typedef std::vector<PairPtr> PairList;
	//! Local typedef for a pointer to a quartet of basis functions
	typedef std::unique_ptr<AbstractBFQuad, BFQuadPool::Deleter> QuadPtr;
	//! Local typedef for a list of basis function quadruples
	typedef std::vector<QuadPtr> QuadList;

	//! Enumeration for the different status flags
	enum StatusFlag
	{
		//! Set when the list of function pairs is current
		PAIRS_CURRENT,
		//! Set when the list of function quadruples is current
		QUADS_CURRENT,
		//! Set when the overlap matrix is computed and up to date
		OVERLAP_CURRENT,
		//! Set when the orthogonalization matrix is computed and up to date
		ORTHO_CURRENT,
		//! Set when the kinetic energy matrix is computed and up to date
		KINETIC_CURRENT,
		//! Set when the nuclear attraction matrix is computed and up to date
		NUC_ATTR_CURRENT,
		//! Set when the electron repulsion integrals are computed and up to date
		ELEC_REP_CURRENT,
		//! The number of status flags
		NR_FLAGS
	};

	/*!
	 * \brief Constructor
	 *
	 * Create a new and empty basis
	 */
	Basis(): _funs(), _quad_pool(), _pairs(), _quads(), _status(),
		_overlap(), _kinetic() {}

	//! Return the number of functions in this basis
	int size() const { return _funs.size(); }

	/*!
	 * \brief Add a basis function
	 *
	 * Add basis function \a bf to this basis.
	 */
	void add(AbstractBF *bf)
	{
		_funs.push_back(BasisFunPtr(bf));
		_status.reset();
	}

	/*!
	 * \brief Return the overlap matrix
	 *
	 * Return the matrix with the overlaps between the fucntions in this
	 * basis, computing it first if necessary.
	 */
	const Eigen::MatrixXd& overlap() const
	{
		if (!_status.test(OVERLAP_CURRENT))
			calcOneElectron();
		return _overlap;
	}
	/*!
	 * \brief Return the orthogonalization matrix
	 *
	 * Return the matrix with the\f$X = S^{-1/2}\f$ which orthogonalizes
	 * this basis.
	 */
	const Eigen::MatrixXd& ortho() const
	{
		if (!_status.test(ORTHO_CURRENT))
			calcOrtho();
		return _ortho;
	}
	/*!
	 * \brief Return the kinetic energy matrix
	 *
	 * Return the matrix with the kinetic energy integrals for the electrons
	 * in this basis, computing it first if necessary.
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
	const Eigen::MatrixXd& nuclearAttraction(const Eigen::MatrixXd& nuc_pos,
		const Eigen::VectorXd& nuc_charge) const
	{
		if (!_status.test(NUC_ATTR_CURRENT))
			calcNuclearAttraction(nuc_pos, nuc_charge);
		return _nuc_attr;
	}

	/*!
	 * \brief Effective electron repulsion
	 *
	 * Compute the two-electron matrix \f$G = J - \frac{1}{2}K\f$, given
	 * the electron density matrix \a D.
	 * \param D The electron density matrix
	 * \param G Place to store the ERI matrix
	 */
	void twoElectron(const Eigen::MatrixXd& D, Eigen::MatrixXd& G) const;
	/*!
	 * \brief Effective electron repulsion
	 *
	 * Compute the Coulomb and exchange matrices, given the electron density
	 * matrix \a D.
	 * \param D The electron density matrix
	 * \param J Place to store the Coulomb repulsion matrix
	 * \param K Place to store the exchange matrix
	 */
	void twoElectron(const Eigen::MatrixXd& D,
		Eigen::MatrixXd& J, Eigen::MatrixXd& K) const;

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
	//! Memory pool for storing basis function quartets
	mutable BFQuadPool _quad_pool;
	//! The list of function pairs
	mutable PairList _pairs;
	//! The list of function quadruples
	mutable QuadList _quads;
	//! Status flags for the basis
	mutable std::bitset<NR_FLAGS> _status;
	//! The overlap matrix for this basis
	mutable Eigen::MatrixXd _overlap;
	//! Orhtogonalization matrix \f$S^{-1/2}\f$
	mutable Eigen::MatrixXd _ortho;
	//! The kinetic energy matrix for this basis
	mutable Eigen::MatrixXd _kinetic;
	//! The nuclear attraction matrix for this basis
	mutable Eigen::MatrixXd _nuc_attr;
	//! The electron repulsion integrals for this basis
	mutable Eigen::ArrayXd _elec_rep;

	//! Create the list of basis function pairs
	void setPairs() const;
	//! Create the list of basis function quadruples
	void setQuads() const;
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
	/*!
	 * \brief Compute the electron repulsion integrals
	 * 
	 * Compute the electron repulsion integrals for this basis
	 */
	void calcElectronRepulsion() const;
	//! Compute the orthogonalization matrix for this basis
	void calcOrtho() const;
	
	double eri(size_t i, size_t j, size_t k, size_t l) const
	{
		return _elec_rep[eriIndex(i, j, k, l)];
	}
	/*!
	 * \brief Return the index of the electron repulsion intergal
	 *    (ij,kl) in the list of ERI's,
	 */
	static size_t eriIndex(size_t i, size_t j, size_t k, size_t l)
	{
		if (i < j)
			std::swap(i, j);
		if (k < l)
			std::swap(k, l);
		size_t ij = (i*(i+1)) / 2 + j;
		size_t kl = (k*(k+1)) / 2 + l;
		if (ij < kl)
			std::swap(ij, kl);
		return (ij*(ij+1))/2 + kl;
	}
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
