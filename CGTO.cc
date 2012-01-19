#include "CGTO.hh"
#include "Dispatcher.hh"
#include "io/manipulators.hh"

const size_t CGTO::cid = Dispatcher::singleton().classID<CGTO>();

std::ostream& CGTO::print(std::ostream& os) const
{
	os << "CGTO (\n" << indent;
	os << "angular momentum: " << lx() << ", " << ly() << ", " << lz() << "\n";
	os << "center: " << center().transpose() << "\n";
	os << "primitives:\n" << indent;
	for (int i = 0; i < size(); i++)
		os << weight(i) << " " << width(i) << "\n";
	os << dedent << dedent << ")";
	return os;
}

double CGTO::eval(const Eigen::Vector3d& pos) const
{
	Eigen::Vector3d dr = pos-center();
	double rn = dr.squaredNorm();
	double res = (-rn * _widths).array().exp().sum();
	for (int p = 0; p < lx(); p++) res *= dr.x();
	for (int p = 0; p < ly(); p++) res *= dr.y();
	for (int p = 0; p < lz(); p++) res *= dr.z();
	return res;
}

