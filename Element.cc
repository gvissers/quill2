#include "Element.hh"

std::ostream& operator<<(std::ostream& os, const Element& elem)
{
	return os << "Element " << elem.symbol() << " (\n"
		<< "name: " << elem.name() << "\n"
		<< "number: " << elem.number() << "\n"
		<< "mass: " << elem.mass() << "\n"
		<< "covalent radius: " << elem.covalentRadius() << "\n"
		<< "Van der Waals radius: " << elem.vanderwaalsRadius() << "\n"
		<< "hydrogen bonding: " << elem.formsHBonds() << "\n"
		<< ")";
}
