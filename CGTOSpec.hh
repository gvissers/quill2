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
	 * on position \a center.
	 * \param weights The weights of the primitives in the contraction
	 * \param ishell  Index of this orbital's shell in the CGTOShellList
	 */
	CGTOSpec(int ishell):
		CGTO(cid, Eigen::Vector3i(tlx, tly, tlz), ishell) {}

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
	os << "norm: " << norm() << "\n";
	os << "shell:\n" << indent << shell() << "\n" << dedent;
	os << dedent << ")";
	return os;
}

#endif // CGTOSPEC_HH
