#ifndef CGTO_HH
#define CGTO_HH

/*!
 * \file CGTO.hh
 * \brief Definition of the CGTO class
 */

#include <Eigen/Dense>
#include "AbstractBF.hh"
#include "exceptions.hh"

/*!
 * \brief Contracted Gaussian Type Orbital
 *
 * Class CGTO defines a contraction of primitive (cartesian) Gaussian
 * type orbitals, each of the form
 * \f[
 * G(x,y,z) = (x-x_c)^{l_x} (y-y_c)^{l_y} (z-z_c)^{l_z}
 *            \exp(-\alpha * |{\bf r}-{\bf r}_c|^2)
 * \f]
 */
class CGTO: public AbstractBF
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new contraction of primitive Gaussian type orbitals
		 * with angular momentum quantum numbers \a ls on position
		 * \a center. The weights of the primitives are given in
		 * weights, the widths in \a widths.
		 */
		CGTO(size_t cid, const Eigen::Vector3i& ls,
			const Eigen::VectorXd& weights,
			const Eigen::VectorXd& widths,
			const Eigen::Vector3d& center):
			AbstractBF(cid), _ls(ls),
			_weights(weights), _widths(widths), _center(center) {}

		//! Return the number of primitives in this contraction
		int size() const { return _weights.size(); }
		//! Return the angular momentum quantum numbers
		const Eigen::Vector3i& ls() const { return _ls; }
		//! Return the angular momentum in the \f$x\f$ direction
		int lx() const { return _ls.x(); }
		//! Return the angular momentum in the \f$y\f$ direction
		int ly() const { return _ls.y(); }
		//! Return the angular momentum in the \f$z\f$ direction
		int lz() const { return _ls.z(); }
		//! Return the weight of the \a i'th primitive
		double weight(int i) const
		{
			checkIndex(i);
			return _weights(i);
		}
		//! Return the width of the \a i'th primitive
		double width(int i) const
		{
			checkIndex(i);
			return _widths(i);
		}
		//! Return the weights of all primitives in this orbital
		const Eigen::VectorXd& weights() const { return _weights; }
		//! Return the widths of all primitives in this orbital
		const Eigen::VectorXd& widths() const { return _widths; }
		//! Return the center position of this orbital
		const Eigen::Vector3d& center() const { return _center; }

		/*!
		 * \brief Print this CGTO
		 *
		 * Write a textual representation of this orbital to output
		 * stream \a os.
		 * \param os The output stream to write to
		 * \return The updated output stream
		 */
		virtual std::ostream& print(std::ostream& os) const;

		/*!
		 * \brief Evaluate this CGTO
		 *
		 * Compute the value of this orbital at position \a pos.
		 * \param pos The point in which to evaluate the orbital
		 * \return The value of this CGTO in \a pos
		 */
		double eval(const Eigen::Vector3d& pos) const;

	private:
		//! The angular momentum quantum numbers
		Eigen::Vector3i _ls;
		//! The weights of the primitives in this contraction
		Eigen::VectorXd _weights;
		//! The widths of the primitive Gaussians
		Eigen::VectorXd _widths;
		//! The center position of this function
		Eigen::Vector3d _center;

		/*!
		 * \brief Check if a primitive index is valid
		 *
		 * Check if \a idx is a valid primitive index for this
		 * contraction. The check is only performed if the program
		 * is compiled with the DEBUG symbol defined.
		 * \param idx The index to check
		 * \exception InvalidIndex thrown when the index is out
		 *    of bounds.
		 */
		void checkIndex(int idx) const
		{
#ifdef DEBUG
			if (idx < 0 || idx >= size())
				throw InvalidIndex(idx);
#endif
		}
};

#endif // CGTO_HH
