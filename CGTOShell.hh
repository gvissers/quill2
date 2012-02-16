#ifndef SHELL_HH
#define SHELL_HH

/*!
 * \file CGTOShell.hh
 * \brief Definition of the CGTOShell class
 */

#include <Eigen/Core>

class CGTOShell
{
public:
	//! Constructor
	CGTOShell(const Eigen::ArrayXd& widths,
		int ipos, const Eigen::Vector3d& center):
		_widths(widths), _ipos(ipos), _center(center) {}

	//! Return the number of primitive widths in this shell
	int size() const { return _widths.size(); }
		
	//! Return the width of the \a i'th primitive
	double width(int i) const { return _widths[i]; }
	//! Return the widths of all primitives in this shell
	const Eigen::ArrayXd& widths() const { return _widths; }
	//! Return the position ID of this shell's center
	int positionId() const { return _ipos; }
	//! Return the center position of this shell
	const Eigen::Vector3d& center() const { return _center; }
	//! Return the \a i coordinate of the center of this shell
	double center(int i) const { return _center[i]; }

	/*!
	 * \brief Print this CGTOShell
	 *
	 * Write a textual representation of this shell to output stream \a os.
	 * \param os The output stream to write to
	 * \return The updated output stream
	 */
	std::ostream& print(std::ostream& os) const;

private:
	//! The widths of the primitive Gaussians
	Eigen::ArrayXd _widths;
	//! Position identifier of the center
	int _ipos;
	//! The center position of this function
	Eigen::Vector3d _center;
};

namespace {

/*!
 * \brief Print a shell
 *
 * Write a textual representation of CGTO shell \a shell to output stream \a os.
 * \param os    The output stream to write to
 * \param shell The shell to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const CGTOShell& shell)
{
	return shell.print(os);
}

} // namespace

#endif // SHELL_HH