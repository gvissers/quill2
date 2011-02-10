#ifndef ZMATRIX_HH
#define ZMATRIX_HH

/*!
 * \file ZMatrix.hh
 * \brief Definition of the ZMatrix class
 */

#include <vector>
#include <map>
#include "io/JobIStream.hh"
#include "Geometry.hh"

/*!
 * \brief Class for Z-matrices
 *
 * Class ZMatrix represents the geometry of a system in Z-matrix format,
 * i.e. in terms of bond lengths and internal angles.
 */
class ZMatrix
{
	public:
		//! Internal representation of the position of a single atom
		struct AtomPos
		{
			//! Element symbol of the atom
			std::string symbol;
			//! Index of the atom to which this atom is connected
			int idx1;
			//! Atom index for defining the bond angle \a idx2 - \a idx1 - this atom
			int idx2;
			//! Atom index defining the reference plane for the torsional angle
			int idx3;
			//! Bond length between this atom and atom \a idx1
			double r;
			//! Angle between the bond and that between atoms \a idx1 and \a idx2
			double theta;
			//! Angle between the bond and the plane defined by atoms \a idx1, \a idx2, and \a idx3.
			double phi;

			//! Constructor
			AtomPos(const std::string& symbol,
				int idx1, double r,
				int idx2, double theta,
				int idx3, double phi):
				symbol(symbol),
				idx1(idx1), idx2(idx2), idx3(idx3),
				r(r), theta(theta), phi(phi)
				{}
		};

		struct InvalidLabel;

		//! Local typedef for a map from element symbol to index
		typedef std::map<std::string, int> IndexMap;
		//! The list of atoms in the matrix
		typedef std::vector<AtomPos> AtomList;

		/*!
		 * \brief Constructor
		 *
		 * Create a new, empty Z-matrix.
		 */
		ZMatrix(): _idx_map(), _positions() {}

		//! Return the number of nuclei defined by this Z-matrix
		int size() const { return _positions.size(); }

		/*!
		 * \brief Print this Z-matrix
		 *
		 * Write a textual representation of this Z-matrix to output
		 * stream \a os.
		 * \param os The output stream to write to
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;
		/*!
		 * \brief Read this Z-matrix
		 *
		 * Read this Z-matrix from input stream \a is, and stop at
		 * the first line that cannot be parsed.
		 * \param is The input stream to read from
		 * \return The updated input stream
		 * \exception InvalidLabel thrown when a line has the correct
		 *    format for a Z-matrix, but refers to a previous atom using
		 *    an undefined label or invalid index.
		 */
		JobIStream& scan(JobIStream& is);

		/*!
		 * \brief Clear the Z-matrix
		 *
		 * Clear this ZMatrix, removing all nuclei and labels.
		 */
		void clear()
		{
			_idx_map.clear();
			_positions.clear();
		}

		/*!
		 * \brief Add a new atom label
		 *
		 * Add a new label for the atom with index \a idx in the
		 * position list to the label map. If the label is already
		 * already in the map, nothing is doen.
		 * \param lbl The label of the atom
		 * \param idx The index of the atom in the position list
		 */
		void addLabel(const std::string& lbl, int idx)
		{
			_idx_map.insert(std::make_pair(lbl, idx));
		}
		/*!
		 * \brief Look up an index by label
		 *
		 * Find the atom index corresponding to atom label \a lbl
		 * in the label map and return it. When the label is not found,
		 * -1 is returned.
		 * \param lbl The atom label to look up
		 * \return The index of the atom in the position list, or -1
		 *     if the label is not found.
		 */
		int findLabel(const std::string& lbl) const
		{
			IndexMap::const_iterator it = _idx_map.find(lbl);
			if (it != _idx_map.end())
				return it->second;
			return -1;
		}

		/*!
		 * \brief Fill a geometry
		 *
		 * Fill system geometry \a geom with the atom positions in
		 * this Z-matrix. The conversion to Cartesian coordinates is
		 * based on Parsons et al., J. Comp. Chem. 26, 1063 (2005).
		 * \param geom The geometry to fill
		 */
		void fillGeometry(Geometry *geom);

	private:
		//! Map from atom label to index in the position list
		IndexMap _idx_map;
		//! The list of nuclei and their coordinates
		AtomList _positions;

		/*!
		 * \brief Check if an atom index is valid
		 *
		 * Check if the atom index \a idx, defined by label \a lbl
		 * corresponds to a known nucleus in the matrix. If not,
		 * an InvalidLabel exception is thrown.
		 * \param idx The index to check
		 * \param lbl The label defining \a idx, or empty when
		 *    the input was numeric
		 * \exception InvalidLabel thrown when the check fails.
		 */
		void checkIndex(int idx, const std::string& lbl) const;
};

//! %Exception thrown when encountering an undefined atom label
struct ZMatrix::InvalidLabel: public Li::Exception
{
	//! Constructor
	InvalidLabel(const std::string& lbl):
		Exception("Invalid or unknown atom label \"" + lbl + "\"") {}
	//! Constructor
	InvalidLabel(int idx): Exception()
	{
		std::ostringstream os;
		os << "Invalid or unknown atom index " << idx;
		setMsg(os.str());
	}
};

namespace {

/*!
 * \brief Print a Z-matrix
 *
 * Write a textual representation of Z-matrix \a zmat to output stream \a os.
 * \param os   The output stream to write to
 * \param zmat The Z-matrix to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const ZMatrix& zmat)
{
	return zmat.print(os);
}

/*!
 * \brief Read a Z-matrix
 *
 * Read Z-matrix \a zmat from input stream \a is, and stop at the first line
 * that cannot be parsed.
 * \param is   The input stream to read from
 * \param zmat The Z-matrix object to read
 * \return The updated input stream
 * \exception InvalidLabel thrown when a line has the correct format for a
 *    Z-matrix, but refers to a previous atom using an undefined label or
 *    invalid index.
 */
inline JobIStream& operator>>(JobIStream& is, ZMatrix& zmat)
{
	return zmat.scan(is);
}

} // namespace

#endif // ZMATRIX_HH
