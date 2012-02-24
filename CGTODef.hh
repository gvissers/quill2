#ifndef CGTODEF_HH
#define CGTODEF_HH

/*!
 * \file CGTODef.hh
 * \brief Definition of the CGTO class
 */

#include "Eigen/Dense"
#include "AbstractBFDef.hh"
#include "CGTOShellList.hh"
#include "CGTODefExpander.hh"

/*!
 * \brief Contracted Gaussian Type Orbital
 *
 * Class CGTODef defines a contraction of primitive (cartesian) Gaussian
 * type orbitals, each of the form
 * \f[
 * G(x,y,z) = (x-x_c)^{l_x} (y-y_c)^{l_y} (z-z_c)^{l_z}
 *            \exp(-\alpha * |{\bf r}-{\bf r}_c|^2)
 * \f]
 * Here, the angular momentum quantum numbers \f$l_x\f$, \f$l_y\f$, and
 * \f$l_z\f$ add up to the template parameter \a l, the total angular momentum
 * of the orbital. The definition only includes the total angular momentum and
 * the widths and weights of the primitives. The center of the function
 * and the angular momentum quantum numbers of the individual directions
 * are applied when the definition is expanded on an atom in the system.
 * \tparam l The total angular momentum of the orbital
 */
template <unsigned int l>
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
		 * \brief Expand this basis function definition
		 *
		 * Expand this definition, adding contracted GTO functions
		 * centered on position \a pos to basis \a basis.
		 * \param ipos  Position identifier
		 * \param pos   Position on which the functions are centered
		 * \param basis The basis to which the functions are added
		 */
		void expand(int ipos, const Eigen::Vector3d& pos,
			Basis *basis) const
		{
			int ishell = CGTOShellList::singleton().addShell(
				l, _widths, ipos, pos);
			CGTODefExpander<l, 0, 0>::exec(_weights, ishell, basis);
		}

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
		Eigen::ArrayXd _widths;
};

template <unsigned int l>
CGTODef<l>::CGTODef(const std::vector< std::pair<double, double> >& ww):
	_weights(ww.size()), _widths(ww.size())
{
	for (size_t i = 0; i < ww.size(); i++)
	{
		_weights[i] = ww[i].first;
		_widths[i] = ww[i].second;
	}
}

template <>
void CGTODef<0>::expand(int ipos, const Eigen::Vector3d& pos,
	Basis* basis) const
{
	int ishell = CGTOShellList::singleton().addShell(0, _widths, ipos, pos);
#if LMAX_SPECIALIZED >= 0
	basis->add(new CGTOSpec<0, 0, 0>(_weights, ishell));
#else
	basis->add(new CGTO(Eigen::Vector3i(0, 0, 0), _weights, ishell));
#endif
}

template <>
void CGTODef<1>::expand(int ipos, const Eigen::Vector3d& pos,
	Basis* basis) const
{
	int ishell = CGTOShellList::singleton().addShell(1, _widths, ipos, pos);
#if LMAX_SPECIALIZED >= 1
	basis->add(new CGTOSpec<1, 0, 0>(_weights, ishell));
	basis->add(new CGTOSpec<0, 1, 0>(_weights, ishell));
	basis->add(new CGTOSpec<0, 0, 1>(_weights, ishell));
#else
	basis->add(new CGTO(Eigen::Vector3i(1, 0, 0), _weights, ishell));
	basis->add(new CGTO(Eigen::Vector3i(0, 1, 0), _weights, ishell));
	basis->add(new CGTO(Eigen::Vector3i(0, 0, 1), _weights, ishell));
#endif
}

template <unsigned int l>
std::ostream& CGTODef<l>::print(std::ostream& os) const
{
	os << "CGTODef<" << l << "> (\n" << indent;
	for (int i = 0; i < size(); i++)
		os << weight(i) << " " << width(i) << "\n";
	return os << dedent << ")";
}

#endif // CGTODEF_HH
