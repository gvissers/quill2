#include <Eigen/Dense>
#include "Geometry.hh"
#include "PeriodicTable.hh"
#include "XYZMatrix.hh"
#include "ZMatrix.hh"
#include "io/manipulators.hh"
#include "exceptions.hh"

void Geometry::setAtom(int idx, const std::string& symbol,
	double x, double y, double z)
{
	checkIndex(idx);

	const Element& elem = PeriodicTable::singleton().findBySymbol(symbol);

	_positions.col(idx) << x, y, z;
	_masses(idx) = elem.mass();
	_charges(idx) = elem.number();
	_symbols[idx] = symbol;
}

std::ostream& Geometry::print(std::ostream& os) const
{
	os << "Geometry (\n" << indent;
	for (int i = 0; i < size(); i++)
		os << _charges(i) << "\t" << _masses(i) << "\t"
			<< symbol(i) << "\t"
			<< position(i).transpose() << "\n";
	os << dedent << ")";
	return os;
}

JobIStream& Geometry::scan(JobIStream& is)
{
	is >> getline;
	if (is.eof())
		throw UnexpectedEOF();

	std::string elem;
	double x, y, z;
	is >> element(elem, false) >> x >> y >> z;
	if (!is.fail())
	{
		is.ungetLastLine();

		XYZMatrix mat;
		is >> mat;
		mat.fillGeometry(this);
	}
	else
	{
		is.ungetLastLine();

		ZMatrix mat;
		is >> mat;
		mat.fillGeometry(this);
	}

	return is;
}

void Geometry::toPrincipalAxes()
{
	// Move center of mass to origin
	Eigen::Vector3d cm = (positions() * masses()) / masses().sum();
	_positions.colwise() -= cm;

	// Compute principal axes
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(inertia());
	const Eigen::Matrix3d& axes = solver.eigenvectors();
	Eigen::Vector3d c = axes.col(0), b = axes.col(1), a = axes.col(0);
	if (a.cross(b).dot(c) < 0)
		// axis system is left-handed, flip one
		a = -a;

	// First two angles put the z axis in place
	double theta = std::acos(c.z());
	double psi = theta > 2 * std::numeric_limits<double>::epsilon()
		? std::atan2(c.y(), c.x()) : 0;
	
	// Rotate the a axis into the xy-plane
	Eigen::Matrix3d R = (Eigen::AngleAxisd(-theta, Eigen::Vector3d::UnitY())
		* Eigen::AngleAxisd(-psi, Eigen::Vector3d::UnitZ()))
		.toRotationMatrix();
	a = R * a;
	
	// Last angle aligns a axis with x (and b with y)
	double phi = std::atan2(a.y(), a.x());
	R = Eigen::AngleAxisd(-phi, Eigen::Vector3d::UnitZ()) * R;
	_positions = R * positions();
}

Eigen::Matrix3d Geometry::inertia() const
{
	Eigen::Matrix3d I;
        for (int j = 0; j < 3; j++)
        {
                for (int i = 0; i <= j; i++)
                {
			I(j,i) = I(i,j) =
				-_positions.row(i).cwiseProduct(_positions.row(j))
					.dot(masses());
		}
        }
        
        double trace = I.trace();
        for (int i = 0; i < 3; i++)
		I(i,i) += trace;
	
	return I;
}