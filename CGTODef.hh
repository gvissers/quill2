#ifndef CGTODEF_HH
#define CGTODEF_HH

/*!
 * \file CGTODef.hh
 * \brief Definition of the CGTO class
 */

#include "Eigen/Dense"
#include "AbstractBFDef.hh"
#include "IndentingOStream.hh"

/*!
 * \brief Contracted Gaussian Type Orbital
 *
 * Class CGTODef defines a contraction of primitive (cartesian) Gaussian
 * type orbitals, each of the form
 * \f[
 * G(x,y,z) = (x-x_c)^i (y-y_c)^j (z-z_c)^k
 *            \exp(-\alpha * |{\bf r}-{\bf r}_c|^2)
 * \f]
 * Here, the angular momentum quantum numbers \f$i\f$, \f$j\f$, and \f$k\f$
 * add up to the template parameter \a l, the total angular momentum of the
 * orbital. The definition only includes the total angular momentum and
 * the widths and weights of the primitives. The center of the function
 * and the angular momentum quantum numbers of the individual directions
 * are applied when the definition is expanded on an atom in the system.
 * \tparam l The total angular momentum of the orbital
 */
template <int l>
class CGTODef: public AbstractBFDef
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new contracted GTO using the weights and widths
		 * in \a ww. In each element of \a ww, the first number is
		 * the weight of the primitive, and  the second number is
		 * the width.
		 * \param ww The weights and widths of the primitives
		 */
		CGTODef(const std::vector< std::pair<double, double> >& ww);

		//! Return the number of functions in this contraction
		int size() const { return _weights.size(); }
		//! Return the weight of the \a i'th primitive
		double weight(int i) const { return _weights[i]; }
		//! Return the width of the \a i'th primitive
		double width(int i) const { return _widths[i]; }

		/*!
		 * \brief Print a basis function
		 *
		 * Print a textual representation of this CGTODef on output
		 * stream \a os.
		 * \param os The output stream to write to
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;

	private:
		//! Weights of the primitive GTOs in this contraction
		Eigen::VectorXd _weights;
		//! Widths of the primitive GTOs in this contraction
		Eigen::VectorXd _widths;
};

template <int l>
CGTODef<l>::CGTODef(const std::vector< std::pair<double, double> >& ww):
	_weights(ww.size()), _widths(ww.size())
{
	for (size_t i = 0; i < ww.size(); i++)
	{
		_weights[i] = ww[i].first;
		_widths[i] = ww[i].second;
	}
}

template <int l>
std::ostream& CGTODef<l>::print(std::ostream& os) const
{
	os << "CGTODef<" << l << "> (\n" << indent;
	for (int i = 0; i < size(); i++)
		os << weight(i) << " " << width(i) << "\n";
	return os << dedent << ")";
}

#endif // CGTODEF_HH
