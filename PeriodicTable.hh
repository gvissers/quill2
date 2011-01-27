#ifndef PERIODICTABLE_HH
#define PERIODICTABLE_HH

/*!
 * \file PeriodicTable.hh
 * \brief Definition of the PeriodicTable class
 */

#include <vector>
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
		//! Local typedef for the list of elements
		typedef std::vector<Element*> ElementList;
		//! Local typedef for the map from element symbol to information
		typedef std::map<std::string, Element*> SymbolElementMap;
		//! Local typedef for the map from element number to information
		typedef std::map<int, Element*> NumberElementMap;
		//! Local typedef for the map from element name to information
		typedef std::map<std::string, Element*> NameElementMap;

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
		//! Destructor
		~PeriodicTable();

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

		//! Check if an element with symbol \a elem is known
		bool exists(const std::string& elem) const
		{
			return _by_symbol.find(elem) != _by_symbol.end();
		}

		//! Find an element by symbol
		const Element& operator[](const std::string& elem) const
		{
			return findBySymbol(elem);
		}
		//! Find an element by symbol
		const Element& findBySymbol(const std::string& elem) const;
		//! Find an element by atomic number
		const Element& findByNumber(int number) const;
		//! Find an element by name
		const Element& findByName(const std::string& name) const;

		/*!
		 * \brief Insert a new element
		 *
		 * Insert a new element. If an element with the same symbol
		 * already exists, it will we replaced.
		 * \param elem The element to insert into the table
		 */
		void insert(const Element& elem);

	private:
		//! The list of known elements
		ElementList _elements;
		//! Map from symbol to element
		SymbolElementMap _by_symbol;
		//! Map from number to element
		NumberElementMap _by_number;
		//! Map from name to element
		NameElementMap _by_name;
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
