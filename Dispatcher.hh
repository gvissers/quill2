#ifndef DISPATCHER_HH
#define DISPATCHER_HH

/*!
 * \file Dispatcher.hh
 * \brief Definition of the Dispatcher class
 */

#include <typeinfo>
#include <tr1/unordered_map>
#include <Singleton.hh>
#include "support.hh"

// Forward declaration
class AbstractBF;

/*!
 * \brief Class for dispatching integral calculations
 *
 * Class Dispatcher is used to dispatch to dispatch calls to functions
 * taking multiple abstract basis functions as arguments (notably integral
 * calculations) to different functions for the underlying concrete basis
 * function. The type of the underlying basis functions is identified using
 * a unique class ID for each basis function class. These type IDs are used
 * as an index a lookup table to find the requested function.
 */
class Dispatcher: public Li::Singleton<Dispatcher, true>
{
	public:
		//! Local typedef for an overlap function
		typedef double (*OverlapFunctionPtr)(const AbstractBF&, const AbstractBF&);
		//! Local typedef for a lokkup table of overlap functions
		typedef std::tr1::unordered_map<std::pair<size_t, size_t>, OverlapFunctionPtr> OverlapMap;

		//! Constructor
		Dispatcher();

		/*!
		 * \brief Generate a class ID
		 *
		 * Compute a class ID for class \a T. The class ID is computed
		 * as a hash of the class name.
		 * \tparam T The type for which to compute an ID
		 */
		template <typename T>
		size_t classID() const { return _hasher(typeid(T).name()); }

		/*!
		 * \brief Compute the overlap between two basis functions
		 *
		 * Compute the overlap between basis functions \a f and \a g,
		 * by finding the overlap function to use, and calling that.
		 * \param f The first basis function
		 * \param g The second basis function
		 * \return The overlap between \a f and \a g
		 */
		double overlap(const AbstractBF& f, const AbstractBF& g) const;

		//! Return the number of type pairs in this dispatcher
		int nrPairs() { return _S_funs.size(); }

	private:
		//! Hash function object for creating class IDs
		std::tr1::hash<std::string> _hasher;
		//! The map of overlap calculation functions
		OverlapMap _S_funs;

		template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2, int lsum>
		friend struct PairFunctionAdder;

		//! Set the overlap function between orbitals with type ids \a id1 and \a id2 to \a fun
		void setOverlapFunction(size_t id1, size_t id2, OverlapFunctionPtr fun)
		{
			_S_funs.insert(std::make_pair(std::make_pair(id1, id2), fun));
		}
};

#endif // DISPATCHER_HH
