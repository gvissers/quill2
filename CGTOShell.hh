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
	/*!
	 * \brief Constructor
	 *
	 * Create a new shell for contracted Gaussian type orbitals with total
	 * angular momentum \a lsum, consisting of primitive Gaussians defined
	 * by \a weights and \a widths, located on center \a center.
	 * \param index   The index of this shell in the CGTOShellList
	 * \param lsum    The total angular momentum of this shell
	 * \param weights The weights of the primitives in the contraction
	 * \param widths  The widths of the primitives
	 * \param ipos    Position ID of the orbital's center
	 * \param center  The orbital center
	 */
	CGTOShell(int index, int lsum, const Eigen::ArrayXd& weights,
		const Eigen::ArrayXd& widths, int ipos,
		const Eigen::Vector3d& center):
		_index(index), _lsum(lsum), _weights(weights), _widths(widths),
		_ipos(ipos), _center(center) { scaleWeights(); }

	//! Return the number of primitive widths in this shell
	int size() const { return _widths.size(); }

	//! Return this shell's index in the CGTOShellList
	int index() const { return _index; }
	//! Return the total angular momentum of the shell
	int lsum() const { return _lsum; }
	//! Return the weight of the \a i'th primitive
	double weight(int i) const { return _weights[i]; }
	//! Return the weights of all primitives in this shell
	const Eigen::ArrayXd& weights() const { return _weights; }
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
	//! Index of this shell in the CGTOShellList
	int _index;
	//! The total angular momentum of the shell
	int _lsum;
	//! The widths of the primitive Gaussians
	Eigen::ArrayXd _weights;
	//! The widths of the primitive Gaussians
	Eigen::ArrayXd _widths;
	//! Position identifier of the center
	int _ipos;
	//! The center position of this function
	Eigen::Vector3d _center;

	//! Scale the weights of the primitives as a first step in normalization
	void scaleWeights();
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