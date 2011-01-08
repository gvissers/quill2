#ifndef BASISSET_HH
#define BASISSET_HH

/*!
 * \file BasisSet.hh
 * \brief definition of the BasisSet class
 */

#include <vector>
#include <map>
#include <Exception.hh>
#include "AbstractBF.hh"
#include "LineGetter.hh"

/*!
 * \brief Class representing a basis set
 *
 * Class BasisSet holds the basis functions to be used for the different
 * elements.
 */
class BasisSet
{
	public:
		//! Typedef for a list of basis functions
		typedef std::vector<AbstractBF*> BFList;
		//! Typedef for the basis function map
		typedef std::map<std::string, BFList> BFMap;
		//! Enumeration type for recognised basis set formats
		enum Format
		{
			//! Try to detect format automatically
			Auto,
			//! Turbomole format
			Turbomole,
			//! Molpro format
			Molpro
		};

		struct UnknownFormat;
		struct ParseError;
		struct UnknownShell;

		/*!
		 * \brief Constructor
		 *
		 * Create a new and empty basis set
		 */
		BasisSet(): _elements() {}

		/*!
		 * \brief Print a basis set
		 *
		 * Print a textual representation of this basis set on
		 * output stream \a os.
		 * \param os The output stream to print on
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;
		/*!
		 * \brief Read a basis set
		 *
		 * Try to read the basis set definition from input stream \a is,
		 * asuming format \a format. If \a format is \c Auto, try to
		 * determine the format from the file contents.
		 * \param format The input format to use
		 * \param is     The input stream to read from
		 */
		template <Format format>
		void read(std::istream& is);
		/*!
		 * \brief Read a basis set
		 *
		 * Try to read the basis set definition from input stream \a is,
		 * trying to determine the format from the file contents.
		 * \param is     The input stream to read from
		 */
		void read(std::istream& is);

		//! Clear this basis set, removing all basis function definitions
		void clear();

	private:
		//! Map from element name to basis function
		BFMap _elements;

		//! Read an element definition in format \a format
		template<Format format>
		void readElement(LineGetter& getter);
		//! Read a single basis function in Turbomole format
		AbstractBF* readTurbomoleBF(LineGetter& getter);

		/*!
		 * \brief Create a contracted Gaussian
		 *
		 * Create a new contracted Gaussian orbital for shell \a shell,
		 * using the weights and widths in \a ww for the primitives.
		 * \param shell Angular momentum of the orbital
		 * \param ww    Widths and weights of the primitives in the
		 *    contraction
		 * \return a new contracted Gaussian orbital
		 * \exception UknownShell thrown when \a shell is not
		 *    recognised.
		 */
		static AbstractBF* contractedGaussian(char shell,
			const std::vector< std::pair<double, double> >& ww);
};

//! %Exception thrown when a file format is not known
struct BasisSet::UnknownFormat: public Li::Exception
{
	//! Constructor
	UnknownFormat(int format): Exception("Unknown absis set format"),
		format(format) {}

	//! The format value that was unknown
	int format;
};

//! %Exception thrown when we fail to parse an input file
struct BasisSet::ParseError: public Li::Exception
{
	//! Constructor
	ParseError(int line, Format format, const std::string& msg):
		Exception(msg), line(line), format(format) {}

	//! The line number where the error occurred
	int line;
	//! The format in which we were trying to read the input
	Format format;
};

//! Exception thrown when a shell code is not known
struct BasisSet::UnknownShell: public Li::Exception
{
	//! Constructor
	UnknownShell(char shell):
		Exception(std::string("Unknown shell '") + shell + "'"),
		shell(shell) {}

	//! The shell code that was not recognised
	char shell;
};

namespace {

inline std::ostream& operator<<(std::ostream& os, const BasisSet& set)
{
	return set.print(os);
}

} // namespace

#endif // BASISSET_HH
