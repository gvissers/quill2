#include <cmath>
#include "CGTOShell.hh"
#include "io/manipulators.hh"

std::ostream& CGTOShell::print(std::ostream& os) const
{
	os << "CGTOShell (\n" << indent;
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
		alpha *= 4 * widths();
	_weights *= alpha.sqrt();
}