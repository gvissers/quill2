#ifndef BASISSET_HH
#define BASISSET_HH

#include <map>
#include <Exception.hh>
#include "AbstractBF.hh"

/*!
 * \brief Class representing a basis set
 *
 * Class BasisSet holds the basis functions to be used for the different
 * elements.
 */
class BasisSet
{
	public:
		//! Typedef for the basis function map
		typedef std::map<std::string, AbstractBF*> BFMap;
		//! Enumeration type for recognised basis set formats
		enum Format
		{
			//! Try to detect format automatically
			Auto,
			//! Turbomole format
			Turbomole
		};

		struct UnknownFormat;

		/*!
		 * \brief Constructor
		 *
		 * Create a new and empty basis set
		 */
		BasisSet(): _functions() {}

		/*!
		 * \brief Read a basis set
		 *
		 * Try to read the basis set definition from input stream \a is,
		 * asuming format \a format. If \a format is \c Auto, try to
		 * determine the format from the file contents.
		 * \param is     The input stream to read from
		 * \param format The input format to use
		 */
		void read(std::istream& is, Format format=Auto);

		//! Clear this basis set, removing all basis function definitions
		void clear();

	private:
		//! Map from element name to basis function
		BFMap _functions;

		//! Read a basis set definition from \a is, trying to determine the format automatically
		void readAuto(std::istream& is);
		//! Read a basis set definition from \a is in Turbomole format
		void readTurbomole(std::istream& is);
};

//! Exception thrown when a file format is not known
struct BasisSet::UnknownFormat: public Li::Exception
{
	//! Constructor
	UnknownFormat(int format): Exception("Unknown absis set format"),
		format(format) {}

	//! The format value that was unknown
	int format;
};

#endif // BASISSET_HH
