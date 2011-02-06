#ifndef MANIPULATORS_HH
#define MANIPULATORS_HH

#include <iostream>

//! Proxy struct for reference to element string
struct ElementRef
{
	//! Reference to the element symbol to read
	std::string& elem;
	//! Whether t check if the element is known in Quill's period table
	bool check_if_exists;

	//! Constructor
	ElementRef(std::string& str, bool check):
		elem(str), check_if_exists(check) {}
};
/*!
 * \brief Manipulator for reading element symbols
 *
 * Extracting \a eref from a \c std::istream reads the next string in the
 * stream into the string referenced by \a eref. If
 * ElementRef::check_if_exists of \a eref is true, and the string read
 * is not a known element symbol, \c failbit is set in \a is.
 * \param is   The input stream to read from
 * \param eref Element reference to the string to read in
 */
std::istream& operator>>(std::istream& is, ElementRef eref);

namespace {

/*!
 * \brief Create a new element reference
 *
 * Create a proxy for the string reference \a str, to be read in as an element
 * symbol. If the optinal parameter \a check is \c true, the symbol read will
 * be checked in the periodic table to see if it corresponds to a known
 * element.
 * \param str String in which to read the symbol
 * \param check Whether to check if the element is known
 */
inline ElementRef element(std::string& str, bool check=true)
{
	return ElementRef(str, check);
}

} // namespace

#endif // MANIPULATORS_HH
