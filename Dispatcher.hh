#ifndef DISPATCHER_HH
#define DISPATCHER_HH

/*!
 * \file Dispatcher.hh
 * \brief Definition of the Dispatcher class
 */

#include <typeinfo>
#include <memory>
#include <unordered_map>
#include <Singleton.hh>
#include "support.hh"

// Forward declarations
class AbstractBF;
class AbstractBFPair;
class AbstractBFQuad;

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
	typedef std::unordered_map<std::pair<size_t, size_t>, PairCreatorPtr> PairMap;
	//! Local typedef for a quartet creation function
	typedef AbstractBFQuad* (*QuadCreatorPtr)(const AbstractBFPair& p, const AbstractBFPair& q);
	//! Local typedef for a lookup table of pair creation functions
	typedef std::unordered_map<std::pair<size_t, size_t>, QuadCreatorPtr> QuadMap;

	//! Constructor
	Dispatcher();

	/*!
	 * \brief Generate a class ID
	 *
	 * Compute a class ID for class \a T. The class ID is computed as a hash
	 * of the class name.
	 * \tparam T The type for which to compute an ID
	 */
	template <typename T>
	size_t classID() const { return _hasher(typeid(T).name()); }

	/*!
	 * \brief Create a pair of basis functions
	 *
	 * Create the basis function pair from functions \a f and \a g, by
	 * looking up the types of \a f and \a g in the lookup table, and
	 * calling the appropriate constructor function.
	 * \param f The first function in the pair
	 * \param g The second function in the pair
	 * \return Pointer to the basis function pair (f,g)
	 */
	std::unique_ptr<AbstractBFPair> pair(const AbstractBF& f,
		const AbstractBF& g) const;
	/*!
	 * \brief Create a quartet of basis functions
	 *
	 * Create the basis function quarter from function pairs \a p and \a q,
	 * by looking up the types of \a p and \a q in the lookup table, and
	 * calling the appropriate constructor function.
	 * \param p The first function pair in the quartet
	 * \param q The second function pair in the quartet
	 * \return Pointer to the basis function quartet (p,q)
	 */
	std::unique_ptr<AbstractBFQuad> quad(const AbstractBFPair& p,
		const AbstractBFPair& q) const;

	//! Return the number of type pairs in this dispatcher
	int nrPairs() { return _pair_funs.size(); }
	//! Return the number of type quartets in this dispatcher
	int nrQuads() { return _quad_funs.size(); }

private:
	//! Hash function object for creating class IDs
	std::hash<std::string> _hasher;
	//! The map of pair creation functions
	PairMap _pair_funs;
	//! The map of quartet creation functions
	QuadMap _quad_funs;

	template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2, int lsum>
	friend struct SpecSpecAdder;
	template <int lx1, int ly1, int lz1>
	friend struct SpecGenericAdder;

	/*!
	 * \brief Set the pair creation function for orbitals with type ids
	 *    \a id1 and \a id2 to \a fun.
	 */
	void setPairCreator(size_t id1, size_t id2, PairCreatorPtr pc)
	{
		_pair_funs.insert(std::make_pair(std::make_pair(id1, id2), pc));
	}
	/*!
	 * \brief Set the quartet creation function for orbital pairs with type
	 *    ids \a id1 and \a id2 to \a fun.
	 */
	void setQuadCreator(size_t id1, size_t id2, QuadCreatorPtr qc)
	{
		_quad_funs.insert(std::make_pair(std::make_pair(id1, id2), qc));
	}
};

#endif // DISPATCHER_HH
