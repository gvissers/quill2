#ifndef PERIODICTABLE_HH
#define PERIODICTABLE_HH

/*!
 * \file PeriodicTable.hh
 * \brief Definition of the PeriodicTable class
 */

#include <map>
#include <Singleton.hh>
#include "Element.hh"

/*!
 * \brief Class for known elements
 *
 * Class PeriodicTable holds all elements known to the program.
 */
class PeriodicTable: public Li::Singleton<PeriodicTable>
{
	public:
		//! Local typedef for the map from element name to information
		typedef std::map<std::string, Element> ElementMap;

		struct UnknownElement;

		/*!
		 * \brief Constructor
		 *
		 * Create a new, empty PeriodicTable
		 */
		PeriodicTable(): Li::Singleton<PeriodicTable>(), _elements() {}
		/*!
		 * \brief Constructor
		 *
		 * Create a new PeriodicTable, reading data from file
		 * \a fname
		 * \param fname Name of the file containing the element data
		 * \exception NoFile throw when file \a fname cannot be opened
		 */
		PeriodicTable(const std::string& fname);

		/*!
		 * \brief Read element data
		 *
		 * Read element data form input stream \a is.
		 * \param is The input stream to read the periodic table from
		 */
		void read(std::istream& is);
		/*!
		 * \brief Print this PeriodicTable
		 *
		 * Print the collection of known elements on output stream
		 * \a os.
		 * \param os    The output stream to write to
		 * \return The updated output stream
		 */
		std::ostream& print(std::ostream& os) const;

		//! Return a read-only reference to an element
		const Element& operator[](const std::string& elem) const;
		//! Return a read-write reference to an element
		Element& operator[](const std::string& elem);

		/*!
		 * \brief Insert a new element
		 *
		 * Insert a new element, overwriting any old element with
		 * the same symbol.
		 * \param elem The element to insert into the table
		 */
		void insert(const Element& elem)
		{
			_elements.insert(std::make_pair(elem.symbol(), elem));
		}

	private:
		//! The elements in the periodic table
		ElementMap _elements;
};

/*!
 * \brief %Exception thrown when trying to look up an element not in the table
 */
struct PeriodicTable::UnknownElement: public Li::Exception
{
	//! Constructor
	UnknownElement(const std::string& elem):
		Exception("Unknown element '" + elem + "'") {}
};

namespace {

/*!
 * \brief Print a PeriodicTable
 *
 * Print the collection of known elements on output stream \a os.
 * \param os    The output stream to write to
 * \param table The PeriodicTable to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const PeriodicTable& table)
{
	return table.print(os);
}

} // namespace

#endif // PERIODICTABLE_HH
