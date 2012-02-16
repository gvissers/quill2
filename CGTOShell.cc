#include <cmath>
#include "CGTOShell.hh"
#include "io/manipulators.hh"

std::ostream& CGTOShell::print(std::ostream& os) const
{
	os << "CGTOShell (\n" << indent;
	os << "center: " << center().transpose() << "\n";
	os << "widths:\n" << indent;
	for (int i = 0; i < widths().size(); i++)
		os << width(i) << "\n";
	os << dedent << dedent << ")";
	return os;
}
