#ifndef DELETER_HH
#define DELETER_HH

/*!
 * \file Deleter.hh
 * \brief Definitions of the Deleter and ArrayDeleter classes
 */

#include <functional>

/*!
 * \brief Class for deleting objects
 *
 * Deleter objects can be used in combination with STL iterator functions
 * to free all objects in an iterator range. Example use:
 * \code
 * std::vector<Object*> obj_vec;
 * // Fill the vector
 * obj_vec.push_back(new Object());
 * obj_vec.push_back(new Object());
 * // Delete all Object instances in the vector
 * std::for_each (obj_vec.begin(), obj_vec.end(), Deleter<Object>());
 * // Note: the vector still contains the pointers to the Objects (which
 * // are now invalid). When obj_vec is still to be used, always clear it.
 * obj_vec.clear();
 * \endcode
 * \note For freeing data associated with an array (allocated with
 * \c new[]) use an ArrayDeleter instead of a Deleter.
 * \sa ArrayDeleter
 */
template <class T>
struct Deleter: public std::unary_function<T*, void>
{
	//! Delete object \a t
	void operator()(T *t) const { delete t; }
};

/*!
 * \brief Class for deleting arrays
 *
 * Struct ArrayDeleter can be used in combination with STL iterator
 * functions to delete a collection of arrays. Example use:
 * \code
 * std::vector<double*> matrix;
 * // Fill the vector
 * matrix.push_back(new double[100]);
 * matrix.push_back(new double[100]);
 * // Delete all arrays in the vector
 * std::for_each (matrix.begin(), matrix.end(), ArrayDeleter<double>());
 * // Remove the (now invalid) pointers from the vector
 * matrix.clear();
 * \endcode
 * \note The template argument to the ArrayDeleter class is the type
 * of an element in the array, not that of the array itself.
 * \note For deleting "normal" objects (i.e. not allocated with \c new[]),
 * use a Deleter object, not an ArrayDeleter.
 * \sa Deleter
 */
template <class T>
struct ArrayDeleter: public std::unary_function<T*, void>
{
	//! Delete array \a at
	void operator()(T *at) const { delete[] at; }
};

#endif // DELETER_HH
