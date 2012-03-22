#include <cmath>
#include "CGTOShell.hh"
#include "Basis.hh"
#include "AngMomIterator.hh"
#include "io/manipulators.hh"

void CGTOShell::expand(Basis& basis) const
{
	for (AngMomIterator it(lsum()); !it.end(); ++it)
		basis.add(new CGTO(it.ls(), *this));
}

std::ostream& CGTOShell::print(std::ostream& os) const
{
	os << "CGTOShell (\n" << indent;
	os << "angular momentum: " << lsum() << "\n";
	os << "center: " << center().transpose() << "\n";
	os << "primitives:\n" << indent;
	for (int i = 0; i < widths().size(); i++)
 		os << weight(i) << " " << width(i) << "\n";
	os << dedent << dedent << ")";
	return os;
}

void CGTOShell::scaleWeights()
{
	Eigen::ArrayXd alpha = 2 * widths() / M_PI;
	alpha *= alpha.sqrt();
	for (int l = 0; l < lsum(); ++l)
		alpha *= 4 * widths() / (2*l+1);
	_weights *= alpha.sqrt();
}
