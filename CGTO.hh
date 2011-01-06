#ifndef CGTO_HH
#define CGTO_HH

#include "Eigen/Dense"
#include "AbstractBF.hh"

/*!
 * \brief Contracted Gaussian Type Orbital
 *
 * Class CGTO represents a contraction of primitive (cartesian) Gaussian
 * type orbitals, each of the form
 * * \f[
 * G(x,y,z) = (x-x_c)^i (y-y_c)^j (z-z_c)^k
 *            \exp(-\alpha * |{\bf r}-{\bf r}_c|^2)
 * \f]
 * Here, the angular momentum quantum numbers $i$, $j$, and $k$ add up to
 * the template parameter \a l, the total angular momentum of the orbital.
 */
template <int l>
class CGTO: public AbstractBF
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new, empty, contracted GTO.
		 */
		CGTO(): _cws() {}

		/*!
		 * \brief Add a function
		 *
		 * Add a primitive GTO with weight \a weight and width \a width
		 * to this contraction.	
		 * \param weight The weight of the GTO in this contraction
		 * \param width  The width ($\alpha$) of the primitive GTO
		 */
		void add(double weight, double width)
		{
			int nf = size();
			_cws.conservativeResize(Eigen::NoChange, nf+1);
			_cws(0, nf) = weight;
			_cws(1, nf) = width;
		}

		//! Return the number of functions in this contraction
		int size() const { return _cws.cols(); }
		//! Return the weight of the \a i'th primitive
		int weight(int i) const { return _cws(0, i); }
		//! Return the width of the \a i'th primitive
		int width(int i) const { return _cws(1, i); }

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
		//! Weights and widths of the primitive GTOs in this contraction
		Eigen::Matrix<double, 2, Eigen::Dynamic> _cws;
};

template <int l>
std::ostream& CGTO<l>::print(std::ostream& os) const
{
	os << "CGTO<" << l << "> (\n";
	for (int i = 0; i < size(); i++)
		os << weight(i) << " " << width(i) << "\n";
	return os << ")";
}

#endif // CGTO_HH
