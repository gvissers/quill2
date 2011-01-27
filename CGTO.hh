#ifndef CGTO_HH
#define CGTO_HH

/*!
 * \file CGTO.hh
 * \brief Definition of the CGTO class
 */

#include "Eigen/Dense"
#include "AbstractBF.hh"
#include "IndentingOStream.hh"

/*!
 * \brief Contracted Gaussian Type Orbital
 *
 * Class CGTO represents a contraction of primitive (cartesian) Gaussian
 * type orbitals, each of the form
 * \f[
 * G(x,y,z) = (x-x_c)^i (y-y_c)^j (z-z_c)^k
 *            \exp(-\alpha * |{\bf r}-{\bf r}_c|^2)
 * \f]
 * Here, the angular momentum quantum numbers \f$i\f$, \f$j\f$, and \f$k\f$
 * add up to the template parameter \a l, the total angular momentum of the
 * orbital.
 * \tparam l The total angular momentum of the orbital
 */
template <int l>
class CGTO: public AbstractBF
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
		CGTO(const std::vector< std::pair<double, double> >& ww);

		//! Return the number of functions in this contraction
		int size() const { return _weights.size(); }
		//! Return the weight of the \a i'th primitive
		double weight(int i) const { return _weights[i]; }
		//! Return the width of the \a i'th primitive
		double width(int i) const { return _widths[i]; }

		/*!
		 * \brief Print a basis function
		 *
		 * Print a textual representation of this CGTO on output
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
CGTO<l>::CGTO(const std::vector< std::pair<double, double> >& ww):
	_weights(ww.size()), _widths(ww.size())
{
	for (size_t i = 0; i < ww.size(); i++)
	{
		_weights[i] = ww[i].first;
		_widths[i] = ww[i].second;
	}
}

template <int l>
std::ostream& CGTO<l>::print(std::ostream& os) const
{
	os << "CGTO<" << l << "> (\n" << indent;
	for (int i = 0; i < size(); i++)
		os << weight(i) << " " << width(i) << "\n";
	return os << dedent << ")";
}

#endif // CGTO_HH
