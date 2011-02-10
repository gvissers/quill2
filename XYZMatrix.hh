#ifndef XYZMATRIX_HH
#define XYZMATRIX_HH

/*!
 * \file XYZMatrix.hh
 * \brief Definition of the XYZMatrix class
 */

#include "Geometry.hh"

/*!
 * \brief Class for XYZ coordinates
 *
 * Class XYZMatrix represents the geometry of a system in Cartesian
 * XYZ coordinates.
 */
class XYZMatrix
{
	public:
		//! Internal representation of the position of a single atom
		struct AtomPos
		{
			//! Element symbol of the atom
			std::string symbol;
			//! \f$x\f$-coordinates of the atom position
			double x;
			//! \f$y\f$-coordinates of the atom position
			double y;
			//! \f$z\f$-coordinates of the atom position
			double z;

			//! Constructor
			AtomPos(const std::string& symbol, double x, double y,
				double z): symbol(symbol), x(x), y(y), z(z) {}
		};
		//! Local typedef for a list of atom positions
		typedef std::vector<AtomPos> AtomList;

		//! Return the number of atoms in the matrix
		int size() const { return _positions.size(); }

		/*!
		 * \brief Print the coordinates
		 *
		 * Write a textual representation of the cartesian atom
		 * positions to output stream \a os
		 * \param os The output stream to write to
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;
		/*!
		 * \brief Read the coordinates
		 *
		 * Read the set of cartesian coordinates from input stream
		 * \a is.
		 * \param is The input stream to read from
		 * \return The updated input stream
		 */
		JobIStream& scan(JobIStream& is);

		/*!
		 * \brief Fill a geometry
		 *
		 * Fill system geometry \a geom with the atom positions in
		 * this matrix.
		 * \param geom The geometry to fill
		 */
		void fillGeometry(Geometry *geom);

	private:
		//! The list of atom positions
		AtomList _positions;
};

namespace {

/*!
 * \brief Print an XYZ matrix
 *
 * Write a textual representation of the cartesian atom positions in \a mat
 * to output stream \a os.
 * \param os  The output stream to write to
 * \param mat The atom coordinates
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const XYZMatrix& mat)
{
	return mat.print(os);
}

/*!
 * \brief Read an XYZMatrix
 *
 * Read a set of cartesian coordinates from input stream \a is and store them
 * in \a mat.
 * \param is  The input stream to read from
 * \param mat The matrix to read into
 * \return The updated input stream
 */
inline JobIStream& operator>>(JobIStream& is, XYZMatrix& mat)
{
	return mat.scan(is);
}

} // namespace

#endif // XYZMATRIX_HH
