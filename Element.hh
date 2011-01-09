#ifndef ELEMENT_HH
#define ELEMENT_HH

/*!
 * \file Element.hh
 * \brief Definition of the Element class
 */

#include <iostream>

/*!
 * \brief A class for chemical elements
 *
 * Class Element holds information about an element, such as its symbol, name,
 * nuclear charge, and various radii.
 */
class Element
{
	public:
		/*!
		 * \brief Constructor
		 *
		 * Create a new dummy element
		 */
		Element(): _symbol(""), _name(""), _number(0), _mass(0.0),
			_covrad(0.0), _vdwrad(0.0), _H_bonding(false) {}
		/*!
		 * \brief Constructor
		 *
		 * Create a new element
		 * \param symbol    The element symbol
		 * \param name      Full name of the element
		 * \param number    The element number (nuclear charge)
		 * \param mass      Mass of the element in amu
		 * \param H_bonding Whether the element forms hydrogen bonds
		 * \param covrad    Covalent radius of the element in bohr
		 * \param vdwrad    Van der Waals radius of the element in bohr
		 */
		Element(const std::string& symbol, const std::string& name,
			unsigned int number, double mass, bool H_bonding,
			double covrad, double vdwrad=0.0): _symbol(symbol),
			_name(name), _number(number), _mass(mass),
			_covrad(covrad), _vdwrad(vdwrad),
			_H_bonding(H_bonding) {}

		//! Return the symbol of this element
		const std::string& symbol() const { return _symbol; }
		//! Return the full name of this element
		const std::string& name() const { return _name; }
		//! Return the element number of this element
		unsigned int number() const { return _number; }
		//! Return The mass of this element
		double mass() const { return _mass; }
		//! Return the covalent radius of this element
		double covalentRadius() const { return _covrad; }
		//! Return the Van der Waals radius of this element
		double vanderwaalsRadius() const { return _vdwrad; }
		//! Return whether the element forms hydorgen bonds
		bool formsHBonds() const { return _H_bonding; }

	private:
		//! The symbol of the element
		std::string _symbol;
		//! The full name of the element
		std::string _name;
		//! The element number
		unsigned int _number;
		//! The mass of this element
		double _mass;
		//! The covalent radius of this element
		double _covrad;
		//! The Van der Waals radius of this element
		double _vdwrad;
		//! Whether the element forms hydrogen bonds
		bool _H_bonding;
};

/*!
 * \brief Print an element
 *
 * Write element \a elem to output stream \a os
 * \param os   The output stream to write to
 * \param elem The element to print
 * \return The updated output stream
 */
std::ostream& operator<<(std::ostream& os, const Element& elem);

#endif // ELEMENT_HH
