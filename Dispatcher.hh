#ifndef DISPATCHER_HH
#define DISPATCHER_HH

/*!
 * \file Dispatcher.hh
 * \brief Definition of the Dispatcher class
 */

#include <typeinfo>
#include <tr1/memory>
#include <tr1/unordered_map>
#include <Singleton.hh>
#include "support.hh"

// Forward declarations
class AbstractBF;
class AbstractBFPair;

/*!
 * \brief Class for looking up basis function pairs
 *
 * Class Dispatcher is used to store function to create pairs of basis
 * functions, retaining the type information of the concrete functions in the
 * pairs. The resulting pairs stored all derive from AbstractBFPair, which
 * specifies the functionality common to all basis function pairs. The type
 * of the underlying basis functions is identified using a unique class ID
 * for each basis function class. These type IDs are used as an index a lookup
 * table to find the requested pair creation function.
 */
class Dispatcher: public Li::Singleton<Dispatcher, true>
{
	public:
		//! Local typedef for a pair creation function
		typedef AbstractBFPair* (*PairCreatorPtr)(const AbstractBF& f, const AbstractBF& g);
		//! Local typedef for a lookup table of pair creation functions
		typedef std::tr1::unordered_map<std::pair<size_t, size_t>, PairCreatorPtr> PairMap;

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
		 * \brief Create a pair of basis functions
		 *
		 * Create the basis function pair from functions \a f and \a g,
		 * by looking up the types of \a f and \a g in the lookup
		 * table, and calling the apprpriate constructor function.
		 * \param f The first function in the pair
		 * \param g The second function in the pair
		 * \return Pointer to the basis function pair (f,g)
		 */
		std::tr1::shared_ptr<AbstractBFPair> pair(const AbstractBF& f,
			const AbstractBF& g) const;

		//! Return the number of type pairs in this dispatcher
		int nrPairs() { return _pair_funs.size(); }

	private:
		//! Hash function object for creating class IDs
		std::tr1::hash<std::string> _hasher;
		//! The map of pair creation functions
		PairMap _pair_funs;

		template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2, int lsum>
		friend struct PairFunctionAdder;

		//! Set the pair creation function for orbitals with type ids \a id1 and \a id2 to \a fun
		void setPairCreator(size_t id1, size_t id2, PairCreatorPtr pc)
		{
			_pair_funs.insert(std::make_pair(std::make_pair(id1, id2), pc));
		}
};

#endif // DISPATCHER_HH
