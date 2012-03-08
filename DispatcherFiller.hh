#ifndef DISPATCHERFILLER_HH
#define DISPATCHERFILLER_HH

/*!
 * \file DispatcherFiller.hh
 * \brief Functions used to initialize the Dispatcher
 */

#include "Dispatcher.hh"
#include "CGTOSpecPair.hh"
#include "CGTOQuad.hh"

/*!
 * \brief Set pair creation function
 *
 * Add the pair creation function for pairs of contracted Gaussians
 * to the dispatcher. This template inserts a pseudo-constructor for
 * a pair of generic CGTOs, specializations of this template may insert
 * pair creators for specialized templates for specific combinations of
 * angular momentum quantum numbers.
 */
template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2, int lsum>
struct SpecSpecAdder
{
	//! Add pair creation function to the dispatcher
	static void add()
	{
		Dispatcher& dispatcher = Dispatcher::singleton();
		size_t id1 = dispatcher.classID< CGTOSpec<lx1, ly1, lz1> >();
		size_t id2 = dispatcher.classID< CGTOSpec<lx2, ly2, lz2> >();
		dispatcher.setPairCreator(id1, id2, CGTOPair::create);
	}
};
//! Specialization for pairs of two s-type functions
template <>
struct SpecSpecAdder<0, 0, 0, 0, 0, 0, 0>
{
	//! Add pair creation function to the dispatcher
	static void add()
	{
		Dispatcher& dispatcher = Dispatcher::singleton();
		size_t id = dispatcher.classID< CGTOSpec<0, 0, 0> >();
		dispatcher.setPairCreator(id, id,
			CGTOSpecPair<0, 0, 0, 0, 0, 0>::create);
	}
};
//! Specialization for (s,p) and (p,s) pairs
template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
struct SpecSpecAdder<lx1, ly1, lz1, lx2, ly2, lz2, 1>
{
	//! Add pair creation function to the dispatcher
	static void add()
	{
		Dispatcher& dispatcher = Dispatcher::singleton();
		size_t id1 = dispatcher.classID< CGTOSpec<lx1, ly1, lz1> >();
		size_t id2 = dispatcher.classID< CGTOSpec<lx2, ly2, lz2> >();
		dispatcher.setPairCreator(id1, id2,
			CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::create);
	}
};
//! Specialization for (s,d), (p,p) and (d,s) pairs
template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
struct SpecSpecAdder<lx1, ly1, lz1, lx2, ly2, lz2, 2>
{
	//! Add pair creation function to the dispatcher
	static void add()
	{
		Dispatcher& dispatcher = Dispatcher::singleton();
		size_t id1 = dispatcher.classID< CGTOSpec<lx1, ly1, lz1> >();
		size_t id2 = dispatcher.classID< CGTOSpec<lx2, ly2, lz2> >();
		dispatcher.setPairCreator(id1, id2,
			CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::create);
	}
};

template <int lx, int ly, int lz>
struct SpecGenericAdder
{
	static void add()
	{
		Dispatcher& dispatcher = Dispatcher::singleton();
		size_t id1 = dispatcher.classID< CGTOSpec<lx, ly, lz> >();
		size_t id2 = dispatcher.classID< CGTO >();
		dispatcher.setPairCreator(id1, id2, CGTOPair::create);
		dispatcher.setPairCreator(id2, id1, CGTOPair::create);
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all allowed angular momenta in the \f$y\f$ direction of the second
 * function, with all other angular momentum quantum numbers fixed.
 */
template <int lx1, int ly1, int lz1, int l2, int lx2, int ly2>
struct PairMapFillerLy2
{
	//! Add pair functions, looping over \a ly2
	static void fill()
	{
		SpecSpecAdder<lx1, ly1, lz1, lx2, ly2, l2-lx2-ly2, lx1+ly1+lz1+l2>::add();
		PairMapFillerLy2<lx1, ly1, lz1, l2, lx2, ly2-1>::fill();
	}
};
//! Specialization for ly2 = 0, base case for recursion
template <int lx1, int ly1, int lz1, int l2, int lx2>
struct PairMapFillerLy2<lx1, ly1, lz1, l2, lx2, 0>
{
	//! Add pair functions for \a ly2 = 0
	static void fill()
	{
		SpecSpecAdder<lx1, ly1, lz1, lx2, 0, l2-lx2, lx1+ly1+lz1+l2>::add();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all allowed angular momenta in the \f$x\f$ direction of the second
 * function, for fixed angular momentum on the first function, and fixed total
 * angular momentum on the second function.
 */
template <int lx1, int ly1, int lz1, int l2, int lx2>
struct PairMapFillerLx2
{
	//! Add pair functions, looping over \a lx2
	static void fill()
	{
		PairMapFillerLy2<lx1, ly1, lz1, l2, lx2, l2-lx2>::fill();
		PairMapFillerLx2<lx1, ly1, lz1, l2, lx2-1>::fill();
	}
};
//! Specialization for lx2 = 0, base case for recursion
template <int lx1, int ly1, int lz1, int l2>
struct PairMapFillerLx2<lx1, ly1, lz1, l2, 0>
{
	//! Add pair functions, for \a lx2 = 0
	static void fill()
	{
		PairMapFillerLy2<lx1, ly1, lz1, l2, 0, l2>::fill();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all supported total angular momenta of the second function, holding
 * the angular momentum of the first function fixed.
 */
template <int lx1, int ly1, int lz1, int l2>
struct PairMapFillerL2
{
	//! Add pair functions, looping over \a l2
	static void fill()
	{
		PairMapFillerLx2<lx1, ly1, lz1, l2, l2>::fill();
		PairMapFillerL2<lx1, ly1, lz1, l2-1>::fill();
	}
};
//! Base case for recursion: total angular momentum 0
template <int lx1, int ly1, int lz1>
struct PairMapFillerL2<lx1, ly1, lz1, 0>
{
	//! Add pair functions, for l2 = 0
	static void fill()
	{
		PairMapFillerLx2<lx1, ly1, lz1, 0, 0>::fill();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all allowed angular momenta in the \f$y\f$ direction of the first
 * function, for fixed total angular momentum and angular momentum in the
 * \f$x\f$ direction on the first function.
 */
template <int l1, int lx1, int ly1, int lmax2>
struct PairMapFillerLy1
{
	//! Add pair functions, looping over \a ly1
	static void fill()
	{
		PairMapFillerL2<lx1, ly1, l1-lx1-ly1, lmax2>::fill();
		SpecGenericAdder<lx1, ly1, l1-lx1-ly1>::add();
		PairMapFillerLy1<l1, lx1, ly1-1, lmax2>::fill();
	}
};
//! Specialization for ly1 = 0, base case for recursion
template <int l1, int lx1, int lmax2>
struct PairMapFillerLy1<l1, lx1, 0, lmax2>
{
	//! Add pair functions, for \a ly1 = 0
	static void fill()
	{
		PairMapFillerL2<lx1, 0, l1-lx1, lmax2>::fill();
		SpecGenericAdder<lx1, 0, l1-lx1>::add();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all allowed angular momenta in the \f$x\f$ direction of the first
 * function, holding the total angular momenta fixed.
 */
template <int l1, int lx1, int lmax2>
struct PairMapFillerLx1
{
	//! Add pair functions, looping over \a lx1
	static void fill()
	{
		PairMapFillerLy1<l1, lx1, l1-lx1, lmax2>::fill();
		PairMapFillerLx1<l1, lx1-1, lmax2>::fill();
	}
};
//! Specialization for lx1 = 0, base case for recursion
template <int l1, int lmax2>
struct PairMapFillerLx1<l1, 0, lmax2>
{
	//! Add pair functions, for \a lx1  = 0
	static void fill()
	{
		PairMapFillerLy1<l1, 0, l1, lmax2>::fill();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all supported total angular momenta of the first function.
 */
template <int l1, int lmax2>
struct PairMapFiller
{
	//! Add pair functions, looping over \a l1
	static void fill()
	{
		PairMapFillerLx1<l1, l1, lmax2>::fill();
		PairMapFiller<l1-1, lmax2>::fill();
	}
};
//! Base case for recursion: total angular momentum 0
template <int lmax2>
struct PairMapFiller<0, lmax2>
{
	//! Add pair functions, for \a l1 = 0
	static void fill()
	{
		PairMapFillerLx1<0, 0, lmax2>::fill();
	}
};

/*!
 * \brief Fill the quad functions in the dispatcher
 *
 * Fill the quad functions in the dispatcher. This template is used to loop
 * over all supported total angular momenta of the second pair, for fixed total
 * angular momentum of the first pair.
 */
template <int l1, int l2>
struct QuadMapFillerL2
{
	//! Add pair functions, looping over \a l2
	static void fill()
	{
		Dispatcher::singleton().setQuadCreator(l1, l2, CGTOQuad::create);
		QuadMapFillerL2<l1, l2-1>::fill();
	}
};
//! Base case for recursion: total angular momentum 0
template <int l1>
struct QuadMapFillerL2<l1, 0>
{
	static void fill()
	{
		Dispatcher& dispatcher = Dispatcher::singleton();
		dispatcher.setQuadCreator(l1, 0, CGTOQuad::create);
	}
};
template <int l2>
struct QuadMapFillerL2<0, l2>
{
	static void fill()
	{
		Dispatcher::singleton().setQuadCreator(0, l2, CGTOQuad::create);
		Dispatcher::singleton().setQuadCreator(size_t(-1), l2, CGTOQuad::create);
		QuadMapFillerL2<0, l2-1>::fill();
	}
};
template <>
struct QuadMapFillerL2<0, 0>
{
	static void fill()
	{
		Dispatcher& dispatcher = Dispatcher::singleton();
		dispatcher.setQuadCreator(0, 0, CGTOQuad::create);
		dispatcher.setQuadCreator(0, size_t(-1), CGTOQuad::create);
		dispatcher.setQuadCreator(size_t(-1), 0, CGTOQuad::create);
		dispatcher.setQuadCreator(size_t(-1), size_t(-1), CGTOQuad::create);
	}
};

/*!
 * \brief Fill the quad functions in the dispatcher
 *
 * Fill the quad functions in the dispatcher. This template is used to loop
 * over all supported total angular momenta of the first pair.
 */
template <int l1, int l2>
struct QuadMapFiller
{
	//! Add pair functions, looping over \a l1
	static void fill()
	{
		QuadMapFillerL2<l1, l2>::fill();
		Dispatcher::singleton().setQuadCreator(l1, size_t(-1),
			CGTOQuad::create);
		QuadMapFiller<l1-1, l2>::fill();
	}
};
//! Base case for recursion: total angular momentum 0
template <int l2>
struct QuadMapFiller<0, l2>
{
	static void fill()
	{
		QuadMapFillerL2<0, l2>::fill();
	}
};

#endif // DISPATCHERFILLER_HH
