#include "io/manipulators.hh"
#include "PeriodicTable.hh"
#include "support.hh"
#include "exceptions.hh"

std::istream& operator>>(std::istream& is, ElementRef ref)
{
	is >> ref.elem;
	if (is.fail())
		return is;

	toUCFirst(ref.elem);
	if (ref.check_if_exists && !PeriodicTable::getSingleton().exists(ref.elem))
		is.setstate(std::ios::failbit);

	return is;
}
