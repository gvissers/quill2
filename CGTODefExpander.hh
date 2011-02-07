#ifndef CGTODEFEXPANDER_HH
#define CGTODEFEXPANDER_HH

#include "CGTO.hh"

/*!
 * \brief Struct for expanding a CGTO definition
 *
 * Struct CGTODevExpander and its specializations are used to generate all
 * different combinations of angular momentum quantum numbers with a
 * sum equal to the total angular momentum of the CGTO definition, and add
 * basis functions with these quantum numbers to the basis.
 */
template <unsigned int lx, unsigned int ly, unsigned int lz>
struct CGTODefExpander
{
	static void exec(const Eigen::VectorXd& weights,
		const Eigen::VectorXd& widths, const Eigen::Vector3d& pos,
		Basis *basis)
	{
		basis->add(std::tr1::shared_ptr<AbstractBF>(
			new CGTO<lx, ly, lz>(weights, widths, pos)));
		CGTODefExpander<lx, ly+1, lz-1>::exec(weights, widths,
			pos, basis);
	}
};

template <unsigned int lx, unsigned int ly>
struct CGTODefExpander<lx, ly, 0>
{
	static void exec(const Eigen::VectorXd& weights,
		const Eigen::VectorXd& widths, const Eigen::Vector3d& pos,
		Basis *basis)
	{
		basis->add(std::tr1::shared_ptr<AbstractBF>(
			new CGTO<lx, ly, 0>(weights, widths, pos)));
		CGTODefExpander<lx+1, 0, ly-1>::exec(weights, widths,
			pos, basis);
	}
};

template <unsigned int lx>
struct CGTODefExpander<lx, 0, 0>
{
	static void exec(const Eigen::VectorXd& weights,
		const Eigen::VectorXd& widths, const Eigen::Vector3d& pos,
		Basis *basis)
	{
		basis->add(std::tr1::shared_ptr<AbstractBF>(
			new CGTO<lx, 0, 0>(weights, widths, pos)));
	}
};

#endif // CGTODEFEXPANDER_HH
