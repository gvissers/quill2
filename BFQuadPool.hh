#ifndef BFQUADPOOL_HH
#define BFQUADPOOL_HH

/*!
 * \file BFQuadPool.hh
 * \brief Definition of the BFQuadPool class
 */

#include <list>
#include "AbstractBFQuad.hh"
#include "MaxQuadSize.hh"
#include "limits.hh"

/*!
 * \brief Memory pool for basis function quartets
 *
 * Class BFQuadPool implements a very simple memory pool for quartets of basis
 * functions. Elements are all of fixed size (which obviously should be at least
 * as large as the size of the largest basis function type), and are never
 * deallocated unless the pool itself is deleted.
 */ 
class BFQuadPool
{
	//! The size of a single slot
	static const size_t elem_size = 16 * ((MaxQuadSize<Limits::lmax_specialized,
		Limits::lmax_specialized>::size + 15) / 16);
	//! The number of quads per chunk
	static const size_t elems_per_chunk = 1023;

	/*!
	 * \brief A single chunk of memory
	 *
	 * The memory pool is implemented as a list of fixed size chunks of
	 * memory, in which the quads are placed. Struct Chunk implements such
	 * a chunk.
	 */
	struct Chunk
	{
		//! The memory chunk itself
		unsigned char _buf[elem_size*elems_per_chunk];
		//! First address beyond the chunk
		unsigned char* _end;
		//! Place where the next quad is stored
		unsigned char* _cur;

		//! Constructor
		Chunk(): _end(_buf+elems_per_chunk*elem_size), _cur(_buf) {}

		//! Return the place in which to store the next quad
		void* alloc();
	};

public:
	/*!
	 * \brief Deleter object for this pool
	 *
	 * Memory allocated from this pool cannot be deallocated using the
	 * normal \c delete operator. Structure like \a std::unique_ptr can
	 * therefore take an optional deleter argument, that implements the
	 * action to be taken when object can be deleted. Class Deleter
	 * does this for objects allocated from the BFQuadPool.
	 */
	struct Deleter
	{
		//! The pool from which to delete objects
		BFQuadPool* _pool;
		//! Constructor
		Deleter(BFQuadPool* pool): _pool(pool) {}
		//! Remove object starting at \a p from the memory pool
		void operator()(void* p) const { _pool->dealloc(p); }
	};

	//! Constructor
	BFQuadPool(): _chunks(), _deleter(this)
	{
		_chunks.push_front(new Chunk());
	}
	//! Destructor
	~BFQuadPool();

	//! Return a pointer where a new basis fucntion quad can be stored
	void* alloc();
	//! Remove a quad from the pool. Currently a noop.
	void dealloc(void *) {}

	//! Return the Deleter object for this pool
	const Deleter& deleter() const { return _deleter; }

private:
	//! The chunks of memory in which the quads are allocated
	std::list<Chunk*> _chunks;
	//! Deleter object that removes quads from this pool
	Deleter _deleter;
};

//! Allocate a new basis function quad from BFQuadPool \a pool
inline void* operator new(size_t /*nbytes*/, BFQuadPool& pool)
{
	return pool.alloc();
}
//! Free the memory allocated at \a p in memory pool \a pool
inline void operator delete(void* p, BFQuadPool& pool)
{
	pool.dealloc(p);
}

#endif // CGTOQUADPOOL_HH