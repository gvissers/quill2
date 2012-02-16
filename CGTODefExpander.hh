#ifndef CGTODEFEXPANDER_HH
#define CGTODEFEXPANDER_HH

/*!
 * \file CGTODefExpander.hh
 * \brief Definition of the CGTODefExpander class
 */

#include "CGTOSpec.hh"
#include "limits.hh"

/*!
 * \brief Struct for expanding a CGTO definition
 *
 * Struct CGTODefExpander and its specializations are used to generate all
 * different combinations of angular momentum quantum numbers with a
 * sum equal to the total angular momentum of the CGTO definition, and add
 * basis functions with these quantum numbers to the basis.
 */
template <unsigned int lx, unsigned int ly, unsigned int lz>
struct CGTODefExpander
{
	/*!
	 * \brief Add a basis function
	 *
	 * Add a contracted Gaussian type orbital with angular momentum
	 * quantum numbers \a lx, \a ly, and \a lz to the basis, and continue
	 * with the next function.
	 * \param weights Weights of the primitives in the contraction
	 * \param ishell  Index of this orbital's shell in the CGTOShellList
	 * \param basis   The basis to be filled
	 */
	static void exec(const Eigen::VectorXd& weights, int ishell,
		Basis *basis)
	{
		if (Limits::lmax_specialized >= int(lx+ly+lz))
			basis->add(new CGTOSpec<lx, ly, lz>(weights, ishell));
		else
			basis->add(new CGTO(Eigen::Vector3i(lx, ly, lz),
				weights, ishell));
		CGTODefExpander<lx, ly-1, lz+1>::exec(weights, ishell, basis);
	}
};

template <unsigned int lx, unsigned int lz>
struct CGTODefExpander<lx, 0, lz>
{
	static void exec(const Eigen::VectorXd& weights, int ishell,
		Basis *basis)
	{
		if (Limits::lmax_specialized >= int(lx+lz))
			basis->add(new CGTOSpec<lx, 0, lz>(weights, ishell));
		else
			basis->add(new CGTO(Eigen::Vector3i(lx, 0, lz), weights,
				ishell));
		CGTODefExpander<lx-1, lz+1, 0>::exec(weights, ishell, basis);
	}
};

template <unsigned int lz>
struct CGTODefExpander<0, 0, lz>
{
	static void exec(const Eigen::VectorXd& weights, int ishell,
		Basis *basis)
	{
		if (Limits::lmax_specialized >= int(lz))
			basis->add(new CGTOSpec<0, 0, lz>(weights, ishell));
		else
			basis->add(new CGTO(Eigen::Vector3i(0, 0, lz), weights,
				ishell));
	}
};

#endif // CGTODEFEXPANDER_HH
