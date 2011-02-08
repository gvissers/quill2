#ifndef CGTO_HH
#define CGTO_HH

/*!
 * \file CGTO.hh
 * \brief Definition of the CGTO class
 */

#include <Eigen/Dense>
#include "AbstractBF.hh"
#include "io/IndentingOStream.hh"

/*!
 * \brief Contracted Gaussian Type Orbital
 *
 * Class CGTO defines a contraction of primitive (cartesian) Gaussian
 * type orbitals, each of the form
 * \f[
 * G(x,y,z) = (x-x_c)^{l_x} (y-y_c)^{l_y} (z-z_c)^{l_z}
 *            \exp(-\alpha * |{\bf r}-{\bf r}_c|^2)
 * \f]
 * \tparam lx The angular momentum quantum number in the \f$x\$ coordinate
 * \tparam ly The angular momentum quantum number in the \f$y\$ coordinate
 * \tparam lz The angular momentum quantum number in the \f$z\$ coordinate
 */
template <int lx, int ly, int lz>
class CGTO: public AbstractBF
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new contraction of primitve Gaussian type orbitals
		 * on position \a center. The weights of the primitives are
		 * given in weights, the widths in \a widths.
		 */
		CGTO(const Eigen::VectorXd& weights,
			const Eigen::VectorXd& widths,
			const Eigen::Vector3d& center):
			_weights(weights), _widths(widths), _center(center) {}

		//! Return the number of primitives in this contraction
		int size() const { return _weights.size(); }
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
		std::ostream& print(std::ostream& os) const;

		/*!
		 * \brief Evaluate this CGTO
		 *
		 * Compute the value of this orbital at position \a pos.
		 * \param pos The point in which to evaluate the orbital
		 * \return The value of this CGTO in \a pos
		 */
		double eval(const Eigen::Vector3d& pos) const;

	private:
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
		 * \param idx The idnex to check
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

template <int lx, int ly, int lz>
std::ostream& CGTO<lx, ly, lz>::print(std::ostream& os) const
{
	os << "CGTO<" << lx << ", " << ly << ", " << lz << "> (\n" << indent;
	os << "center: " << _center.transpose() << "\n";
	os << "primitives:\n" << indent;
	for (int i = 0; i < size(); i++)
		os << weight(i) << " " << width(i) << "\n";
	os << dedent << dedent << ")";
	return os;
}

template <int lx, int ly, int lz>
double CGTO<lx, ly, lz>::eval(const Eigen::Vector3d& pos) const
{
	Eigen::Vector3d dr = pos-center();
	double rn = dr.squaredNorm();
	double res = (-rn * _widths).array().exp().sum();
	for (int p = 0; p < lx; p++) res *= dr.x();
	for (int p = 0; p < ly; p++) res *= dr.y();
	for (int p = 0; p < lz; p++) res *= dr.z();
	return res;
}

#endif // CGTO_HH
