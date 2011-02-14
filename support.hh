#ifndef SUPPORT_HH
#define SUPPORT_HH

/*!
 * \file support.hh
 * \brief Various support functions and macros for Quill
 */

#include <string>
#include <algorithm>
#include <tr1/functional>
#include <cmath>

using namespace std::tr1::placeholders;

#define NON_COPYABLE(class_name) \
	private:\
		class_name(const class_name&);\
		class_name& operator=(const class_name&);\

// NON_COPYABLE


/*!
 * \brief Remove characters from a string
 *
 * Remove all characters matching condition \a cond from string \a str, and
 * return the updated string.
 * \tparam Cond The type of the condition
 * \param str  The string from which to remove characters
 * \param cond The condition to match
 * \return The string with matching characters removed
 */
template <typename Cond>
std::string& remove(std::string& str, Cond cond)
{
	std::string::iterator new_end = std::remove_if(str.begin(), str.end(),
		cond);
	str.erase(new_end, str.end());
	return str;
}
/*!
 * \brief Remove characters from the end of a string
 *
 * Remove characters at the end of string \a str matching condition \a cond.
 * \tparam Cond The type of the condition
 * \param str  The string from which to remove characters
 * \param cond The condition to match
 * \return The trimmed string
 */
template <typename Cond>
std::string& rtrim(std::string& str, Cond cond)
{
	// Okay, the tr1::bind christmas tree below is rather sick. The problem
	// is that we need to negate the condition in order to find the
	// last character *not* matching it. However, std::not1 will only
	// work on objects that define a type argument_type, which Cond
	// does not necessarily do (most notably, when Cond itself is the
	// result of a tr1::bind).
	std::string::reverse_iterator it = std::find_if(str.rbegin(), str.rend(),
		std::tr1::bind(std::logical_not<bool>(), std::tr1::bind(cond, _1)));
	str.erase(it.base(), str.end());
	return str;
}

/*!
 * \brief Create a string with only the first character upper case
 *
 * Create a copy of \a str, where the first character is capitalized, and
 * all following characters are lower case.
 * \param str The string to transform
 * \return string with only the first character in upper case
 */
std::string ucFirst(const std::string& str);
/*!
 * \brief Convert to string with only the first character upper case
 *
 * Change the case of the characters in \a str such that only the first
 * character is in upper case.
 * \param str The string to transform
 * \return reference to the updated string
 */
std::string& toUCFirst(std::string& str);

/*!
 * \brief Left rotate
 *
 * Rotate number \a x to the left by \a n bits, without carry.
 * \tparam T The type of the number
 * \param x  The number to rotate
 * \param n  The size of the rotation
 * \return The rotated number
 */
template <typename T>
inline T lrot(T x, unsigned int n)
{
	// For unsigned types, the mask can be omitted
	T mask = (1 << n) - 1;
	return (x << n) | ((x >> (sizeof(T)-n)) & mask);
}
template <>
inline size_t lrot<size_t>(size_t x, unsigned int n)
{
	return (x << n) | (x >> (sizeof(size_t)-n));
}

namespace std
{
namespace tr1
{

/*!
 * \brief Hash function for a pair of size_t
 *
 * No standard hash function for pairs exists, so we roll our own by
 * \c xor'ing the values. Since we want to differentiate betweem \c pair(a, b)
 * and \c pair(b, a), we rotate the second value by half the bit length first.
 */
template <>
struct hash< std::pair<size_t, size_t> >
{
	//! Evaluate the hash function
	size_t operator()(const std::pair<size_t, size_t>& cids) const
	{
		return cids.first ^ lrot(cids.second, 4*sizeof(size_t));
	}
};

}
}

namespace
{

//! Convert an angle in degrees to radians
inline double degToRad(double a)
{
	return a * M_PI / 180;
}

//! Convert an angle in radians to degrees
inline double radToDeg(double a)
{
	return a * 180 / M_PI;
}

} // namespace

#endif // SUPPORT_HH
