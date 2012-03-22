#include "CGTO.hh"
#include "Dispatcher.hh"
#include "quillmath.hh"
#include "io/manipulators.hh"

const size_t CGTO::cid = Dispatcher::singleton().classID<CGTO>();

std::ostream& CGTO::print(std::ostream& os) const
{
	os << "CGTO (\n" << indent;
	os << "angular momentum: " << lx() << ", " << ly() << ", " << lz() << "\n";
	os << "shell:\n" << indent << shell() << "\n" << dedent;
	os << dedent << ")";
	return os;
}

double CGTO::eval(const Eigen::Vector3d& pos) const
{
	Eigen::Vector3d dr = pos-center();
	double rn = dr.squaredNorm();
	double res = (-rn * widths()).qexp().sum();
	for (int p = 0; p < lx(); p++) res *= dr.x();
	for (int p = 0; p < ly(); p++) res *= dr.y();
	for (int p = 0; p < lz(); p++) res *= dr.z();
	return res;
}
