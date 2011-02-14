#ifndef CGTOSPEC_HH
#define CGTOSPEC_HH

/*!
 * \file CGTOSpec.hh
 * \brief Definition of the CGTOSpec class
 */

#include "CGTO.hh"
#include "io/manipulators.hh"

/*!
 * Template class CGTOSpec defines a contraction of primitive (cartesian)
 * Gaussian type orbitals, each of the form
 * \f[
 * G(x,y,z) = (x-x_c)^{l_x} (y-y_c)^{l_y} (z-z_c)^{l_z}
 *            \exp(-\alpha * |{\bf r}-{\bf r}_c|^2)
 * \f]
 * where the angular momentum quantum numbers \f$l_x\f$, \f$l_y\f$ and
 * \f$l_z\f$ are template parameters to the class type. Most of the
 * functionality of this class is provided by its base class CGTO; class
 * CGTOSpec is only used to differentiate between orbitals with different
 * angular momentum quantum numbers on a type level.
 * \tparam tlx The angular momentum quantum number in the \f$x\f$ coordinate
 * \tparam tly The angular momentum quantum number in the \f$y\f$ coordinate
 * \tparam tlz The angular momentum quantum number in the \f$z\f$ coordinate
 */
template <int tlx, int tly, int tlz>
struct CGTOSpec: public CGTO
{
	//! Unique class ID, used in looking up integral calculation functions
	static const size_t cid;

	/*!
	 * \brief Constructor
	 *
	 * Create a new contraction of primitive Gaussian type orbitals
	 * on position \a center. The weights of the primitives are
	 * given in weights, the widths in \a widths.
	 */
	CGTOSpec(const Eigen::VectorXd& weights, const Eigen::VectorXd& widths,
		const Eigen::Vector3d& center):
		CGTO(cid, Eigen::Vector3i(tlx, tly, tlz), weights, widths, center) {}

	/*!
	 * \brief Print this orbital
	 *
	 * Write a textual representation of this orbital to output
	 *  stream \a os.
	 * \param os The output stream to write to
	 * \return The updated output stream
	 */
	std::ostream& print(std::ostream& os) const;
};

#include "Dispatcher.hh"

template <int tlx, int tly, int tlz>
const size_t CGTOSpec<tlx, tly, tlz>::cid
	= Dispatcher::singleton().classID< CGTOSpec<tlx, tly, tlz> >();

template <int tlx, int tly, int tlz>
std::ostream& CGTOSpec<tlx, tly, tlz>::print(std::ostream& os) const
{
	os << "CGTOSpec<" << tlx << ", " << tly << ", " << tlz << "> (\n" << indent;
	os << "center: " << center().transpose() << "\n";
	os << "primitives:\n" << indent;
	for (int i = 0; i < size(); i++)
		os << weight(i) << " " << width(i) << "\n";
	os << dedent << dedent << ")";
	return os;
}

#endif // CGTOSPEC_HH
