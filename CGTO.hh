#ifndef CGTO_HH
#define CGTO_HH

/*!
 * \file CGTO.hh
 * \brief Definition of the CGTO class
 */

#include <Eigen/Dense>
#include "AbstractBF.hh"
#include "CGTOShellList.hh"
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
	//! Unique class ID, used in looking up integral calculation functions
	static const size_t cid;

	/*!
	 * \brief Constructor
	 *
	 * Create a new contraction of primitive Gaussian type orbitals with
	 * angular momentum quantum numbers \a ls on position \a center.
	 * \param ls     Angular momentum of the orbital
 	 * \param shell  This orbital's shell
	 */
	CGTO(const Eigen::Array3i& ls, const CGTOShell& shell):
		AbstractBF(cid), _ls(ls), _shell(shell) {}

	//! Return the number of primitives in this contraction
	int size() const { return shell().size(); }
	//! Return the angular momentum quantum numbers
	const Eigen::Array3i& ls() const { return _ls; }
	//! Return the angular momentum in the \f$x\f$ direction
	int lx() const { return _ls.x(); }
	//! Return the angular momentum in the \f$y\f$ direction
	int ly() const { return _ls.y(); }
	//! Return the angular momentum in the \f$z\f$ direction
	int lz() const { return _ls.z(); }
	//! Return the angular momentum in the \a i direction
	int l(int i) const { return _ls[i]; }
	//! Return the total angular momentum
	int lsum() const { return _ls.sum(); }
	//! Return the shell of this orbital
	const CGTOShell& shell() const { return _shell; }
	//! Return the width of the \a i'th primitive
	double width(int i) const { return shell().width(i); }
	//! Return the widths of all primitives in this orbital
	const Eigen::ArrayXd& widths() const { return shell().widths(); }
	//! Return the position ID of this orbital's center
	int positionId() const { return shell().positionId(); }
	//! Return the center position of this orbital
	const Eigen::Vector3d& center() const { return shell().center(); }
	//! Return the \a i coordinate of the center of this orbital
	double center(int i) const { return shell().center(i); }

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
	//! The angular momentum quantum numbers
	Eigen::Array3i _ls;
	//! The shell parameters: widths, weights, and position
	const CGTOShell& _shell;

	/*!
	 * \brief Check if a primitive index is valid
	 *
	 * Check if \a idx is a valid primitive index for this contraction. The
	 * check is only performed if the program is compiled with the DEBUG
	 * symbol defined.
	 * \param idx The index to check
	 * \exception InvalidIndex thrown when the index is out of bounds.
	 */
#ifdef DEBUG
	void checkIndex(int idx) const
	{
		if (idx < 0 || idx >= size())
			throw InvalidIndex(idx);
	}
#else
	void checkIndex(int) const {}
#endif
};

#endif // CGTO_HH
