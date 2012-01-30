#ifndef GEOMETRY_HH
#define GEOMETRY_HH

/*!
 * \file Geometry.hh
 * \brief Definition of the Geometry class
 */

#include <vector>
#include <Eigen/Core>
#include "io/JobIStream.hh"

/*!
 * \brief The geometry of the system
 *
 * Class Geometry hold the positions and atom types of the nuclei in the
 * system.
 */
class Geometry
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new, empty Geometry
		 */
		Geometry(): _positions(), _masses(), _charges(), _symbols() {}

		//! Return the number of atoms in the Geometry
		int size() const { return _positions.cols(); }

		//! Return the element symbol of atom \a idx
		const std::string& symbol(int idx) const
		{
			checkIndex(idx);
			return _symbols[idx];
		}

		//! Return the nuclear charges of the atoms
		const Eigen::VectorXd& charges() const
		{
			return _charges;
		}
		//! Return the total charge of the the nuclei
		double totalCharge() const
		{
			return _charges.sum();
		}
		//! Return the energy due to mutual repulsion of the nuclei
		double nuclearRepulsion() const;

		//! Return the masses of the atoms
		const Eigen::VectorXd& masses() const
		{
			return _masses;
		}

		//! Return the position of atom \a idx
		Eigen::Vector3d position(int idx) const
		{
			checkIndex(idx);
			return _positions.col(idx);
		}
		//! Return tthe positions of the atoms
		const Eigen::MatrixXd& positions() const
		{
			return _positions;
		}
		
		//! Compute the moment of inertia tensor for this molecule
		Eigen::Matrix3d inertia() const;

		/*!
		 * \brief Change the size of the Geometry
		 *
		 * Change the number of atoms in the Geometry to \a n. Any
		 * information currently in the Geomtery will be lost.
		 */
		void resize(int n)
		{
			_positions.resize(3, n);
			_masses.resize(n);
			_charges.resize(n);
			_symbols.resize(n);
		}
		/*!
		 * \brief Set an atom
		 *
		 * On position \a idx in this Geometry, place an atom with
		 * element symbol \a symbol at position (\a x, \a y, \a z).
		 * Any previous atom with the same index will be overwritten.
		 * \param idx    The index of the atom in this Geometry
		 * \param symbol The element symbol of the atom
		 * \param x      The x-coordinate of the atom
		 * \param y      The y-coordinate of the atom
		 * \param z      The z-coordinate of the atom
		 * \exception UnknownElement throw when the element symbol
		 *    cannot be found in the periodic table.
		 */
		void setAtom(int idx, const std::string& symbol, double x,
			double y, double z);
		/*!
		 * \brief Update the position of an atom
		 *
		 * Set the position of the atom at index \a idx to \a pos.
		 * \param idx The index of the atom in this Geometry
		 * \param pos The new position of the atom
		 */
		void setPosition(int idx, const Eigen::Vector3d& pos)
		{
			checkIndex(idx);
			_positions.col(idx) = pos;
		}
		/*!
		 * \brief Rotate to a principal axis frame
		 *
		 * Shift and rotate the system so that its center of mass is
		 * located at the origin, and its principal axes align with the
		 * \f$x\f$, \f$y\f$ and \f$z\f$ axes.
		 */
		void toPrincipalAxes();

		/*!
		 * \brief print this Geometry
		 *
		 * Write a textual representation of this Geometry to output
		 * stream \a os.
		 * \param os The output stream to write to
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;
		/*!
		 * \brief Read a Geometry
		 *
		 * Read a description of this geometry from input stream \a is.
		 * The geometry can be given either in XYZ coordinates, or
		 * in Z-matrix format.
		 * \param is The input stream to read from
		 * \return The updated input stream
		 */
		JobIStream& scan(JobIStream& is);

	private:
		//! The positions of the atoms
		Eigen::MatrixXd _positions;
		//! The mass of the atoms
		Eigen::VectorXd _masses;
		//! The charges of the atoms
		Eigen::VectorXd _charges;
		//! The atom symbols
		std::vector<std::string> _symbols;

		/*!
		 * \brief Check if an atom index is valid
		 *
		 * Check if \a idx is a valid atom index in this Geometry. The
		 * check is only performed if the program is compiled with the
		 * DEBUG symbol defined.
		 * \param idx The atom index to check
		 * \exception InvalidIndex thrown when the index is out
		 *    of bounds
		 */
		void checkIndex(int idx) const
		{
#ifdef DEBUG
			if (idx < 0 || idx >= size())
				throw InvalidIndex(idx);
#endif
		}
};

namespace {

/*!
 * \brief print a Geometry
 *
 * Write a textual representation of Geometry \a geom to output stream \a os.
 * \param os   The output stream to write to
 * \param geom The Geometry to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const Geometry& geom)
{
	return geom.print(os);
}

/*!
 * \brief Read a Geometry
 *
 * Read a description of a geometry from input stream \a is and store it
 * in \a geom. The geometry can be given either in XYZ coordinates, or in
 * Z-matrix format.
 * \param is   The input stream to read from
 * \param geom The Geomtery to read
 * \return The updated input stream
 */
inline JobIStream& operator>>(JobIStream& is, Geometry& geom)
{
	return geom.scan(is);
}

} // namespace

#endif // GEOMETRY_HH
