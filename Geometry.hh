#ifndef GEOMETRY_HH
#define GEOMETRY_HH

/*!
 * \file Geometry.hh
 * \brief Definition of the Geometry class
 */

#include <string>
#include <vector>
#include <Eigen/Core>
#include "JobIStream.hh"

/*!
 * \brief Class for a system's geometry
 *
 * Class geometry holds a number of nuclei that each hold their own positions.
 * It describes the current geometry of the molecule or other system that
 * is being computed.
 */
class Geometry
{
	public:
		//! Struct describing a nucleus
		struct Nucleus
		{
			//! Element symbol of the nucleus
			std::string symbol;
			//! Position of the nucleus
			Eigen::Vector3d position;

			//! Constructor
			Nucleus(const std::string& symbol, double x,
				double y, double z):
				symbol(symbol), position(x, y, z) {}
		};
		//! Local typedef for the list of nuclei
		typedef std::vector<Nucleus> NucleusList;

		/*!
		 * \brief Constructor
		 *
		 * Create a new Geometry without any nuclei.
		 */
		Geometry(): _nuclei() {}

		/*!
		 * \brief Print this Geometry
		 *
		 * Write a textual representation of this Geometry to output
		 * stream \a os.
		 * \param os The output stream to write to
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;
		/*!
		 * \brief Read this Geometry
		 *
		 * Read the nuclei and their positions from input stream \a is.
		 * The input can be in XYZ or Z-matrix format.
		 * \param is The input stream to read from
		 * \return The updated input stream
		 * \exception UnexpectedEOF thrown when no data is available
		 */
		JobIStream& scan(JobIStream& is);

	private:
		//! The list of nuclei in the system
		NucleusList _nuclei;

		/*!
		 * \brief Read a geometry in XYZ format
		 *
		 * Read a geometry in XYZ format from input stream \a is.
		 * When not a single element can be read, an exception is
		 * thrown.
		 * \param is The input stream to read from
		 * \return The updated input stream
		 * \exception UnexpectedEOF  thrown when end of file was
		 *    reached before a line of input could be read
		 * \exception ParseError     thrown when the first line read
		 *    does not describe an element.
		 * \exception UnknownElement thrown when the element symbol
		 *    is not recognised by Quill.
		 */
		JobIStream& scanXYZ(JobIStream& is);
		/*!
		 * \brief Read a single XYZ line
		 *
		 * Read a single XYZ line from input stream \a is. When failing
		 * to parse a line, the result depends on the value of
		 * \a except: if \a except is \c true, an exception is thrown,
		 * otherwise \c false is returned.
		 * \param is The input stream to read from
		 * \param except Whether to throw an exception when failing
		 *    to read an element.
		 * \return \c true when an element was successfully read,
		 *    \c false otherwise.
		 */
		bool scanXYZLine(JobIStream& is, bool except=false);
		JobIStream& scanZMatrix(JobIStream& is);
};

namespace {

/*!
 * \brief Print a Geometry
 *
 * Write a textual representation of Geometry \a geom to output stream \a os.
 * \param os   The output stream to write to
 * \param geom The Geometry object to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const Geometry& geom)
{
	return geom.print(os);
}
/*!
 * \brief Read a Geometry
 *
 * Read the nuclei and their positions from input stream \a is and store them
 * in Geometry \a geom. The input can be in XYZ or Z-matrix format.
 * \param is   The input stream to read from
 * \param geom Place to store the geometry
 * \return The updated input stream
 */
inline JobIStream& operator>>(JobIStream& is, Geometry& geom)
{
	return geom.scan(is);
}

} // namespace

#endif // GEOMETRY_HH
