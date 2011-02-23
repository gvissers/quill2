#ifndef DISPATCHERFILLER_HH
#define DISPATCHERFILLER_HH

/*!
 * \file DispatcherFiller.hh
 * \brief Functions used to initialize the Dispatcher
 */

#include "Dispatcher.hh"
#include "CGTOSpecPair.hh"
#include "gaussint/gto_kinetic.hh"
#include "gaussint/gto_one_elec.hh"

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
struct PairFunctionAdder
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
struct PairFunctionAdder<0, 0, 0, 0, 0, 0, 0>
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
struct PairFunctionAdder<lx1, ly1, lz1, lx2, ly2, lz2, 1>
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
struct PairFunctionAdder<lx1, ly1, lz1, lx2, ly2, lz2, 2>
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

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all allowed angular momenta in the \f$y\f$ direction of the second
 * function, with all other angular momentum quantum numbers fixed.
 */
template <int l1, int lx1, int ly1, int l2, int lx2, int ly2>
struct PairMapFillerLy2
{
	//! Add pair functions, looping over \a ly2
	static void fill()
	{
		PairFunctionAdder<lx1, ly1, l1-lx1-ly1, lx2, ly2, l2-lx2-ly2, l1+l2>::add();
		PairMapFillerLy2<l1, lx1, ly1, l2, lx2, ly2-1>::fill();
	}
};
//! Specialization for ly2 = 2
template <int l1, int lx1, int ly1, int l2, int lx2>
struct PairMapFillerLy2<l1, lx1, ly1, l2, lx2, 2>
{
	//! Add pair functions, looping over \a ly2 = 2, 1, 0
	static void fill()
	{
		PairFunctionAdder<lx1, ly1, l1-lx1-ly1, lx2, 2, l2-lx2-2, l1+l2>::add();
		PairFunctionAdder<lx1, ly1, l1-lx1-ly1, lx2, 1, l2-lx2-1, l1+l2>::add();
		PairFunctionAdder<lx1, ly1, l1-lx1-ly1, lx2, 0, l2-lx2, l1+l2>::add();
	}
};
//! Specialization for ly2 = 1
template <int l1, int lx1, int ly1, int l2, int lx2>
struct PairMapFillerLy2<l1, lx1, ly1, l2, lx2, 1>
{
	//! Add pair functions, looping over \a ly2 = 1, 0
	static void fill()
	{
		PairFunctionAdder<lx1, ly1, l1-lx1-ly1, lx2, 1, l2-lx2-1, l1+l2>::add();
		PairFunctionAdder<lx1, ly1, l1-lx1-ly1, lx2, 0, l2-lx2, l1+l2>::add();
	}
};
//! Specialization for ly2 = 0, base case for recursion
template <int l1, int lx1, int ly1, int l2, int lx2>
struct PairMapFillerLy2<l1, lx1, ly1, l2, lx2, 0>
{
	//! Add pair functions for \a ly2 = 0
	static void fill()
	{
		PairFunctionAdder<lx1, ly1, l1-lx1-ly1, lx2, 0, l2-lx2, l1+l2>::add();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all allowed angular momenta in the \f$x\f$ direction of the second
 * function, for fixed angular momentum on the first function, and fix total
 * angular momentum on the second function.
 */
template <int l1, int lx1, int ly1, int l2, int lx2>
struct PairMapFillerLx2
{
	//! Add pair functions, looping over \a lx2
	static void fill()
	{
		PairMapFillerLy2<l1, lx1, ly1, l2, lx2, l2-lx2>::fill();
		PairMapFillerLx2<l1, lx1, ly1, l2, lx2-1>::fill();
	}
};
//! Specialization for lx2 = 2
template <int l1, int lx1, int ly1, int l2>
struct PairMapFillerLx2<l1, lx1, ly1, l2, 2>
{
	//! Add pair functions, looping over \a lx2 = 2, 1, 0
	static void fill()
	{
		PairMapFillerLy2<l1, lx1, ly1, l2, 2, l2-2>::fill();
		PairMapFillerLy2<l1, lx1, ly1, l2, 1, l2-1>::fill();
		PairMapFillerLy2<l1, lx1, ly1, l2, 0, l2>::fill();
	}
};
//! Specialization for lx2 = 1
template <int l1, int lx1, int ly1, int l2>
struct PairMapFillerLx2<l1, lx1, ly1, l2, 1>
{
	//! Add pair functions, looping over \a lx2 = 1, 0
	static void fill()
	{
		PairMapFillerLy2<l1, lx1, ly1, l2, 1, l2-1>::fill();
		PairMapFillerLy2<l1, lx1, ly1, l2, 0, l2>::fill();
	}
};
//! Specialization for lx2 = 0, base case for recursion
template <int l1, int lx1, int ly1, int l2>
struct PairMapFillerLx2<l1, lx1, ly1, l2, 0>
{
	//! Add pair functions, for \a lx2 = 0
	static void fill()
	{
		PairMapFillerLy2<l1, lx1, ly1, l2, 0, l2>::fill();
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
template <int l1, int lx1, int ly1, int l2>
struct PairMapFillerLy1
{
	//! Add pair functions, looping over \a ly1
	static void fill()
	{
		PairMapFillerLx2<l1, lx1, ly1, l2, l2>::fill();
		PairMapFillerLy1<l1, lx1, ly1-1, l2>::fill();
	}
};
//! Specialization for ly1 = 2
template <int l1, int lx1, int l2>
struct PairMapFillerLy1<l1, lx1, 2, l2>
{
	//! Add pair functions, looping over \a ly1 = 2, 1, 0
	static void fill()
	{
		PairMapFillerLx2<l1, lx1, 2, l2, l2>::fill();
		PairMapFillerLx2<l1, lx1, 1, l2, l2>::fill();
		PairMapFillerLx2<l1, lx1, 0, l2, l2>::fill();
	}
};
//! Specialization for ly1 = 1
template <int l1, int lx1, int l2>
struct PairMapFillerLy1<l1, lx1, 1, l2>
{
	//! Add pair functions, looping over \a ly1 = 1, 0
	static void fill()
	{
		PairMapFillerLx2<l1, lx1, 1, l2, l2>::fill();
		PairMapFillerLx2<l1, lx1, 0, l2, l2>::fill();
	}
};
//! Specialization for ly1 = 0, base case for recursion
template <int l1, int lx1, int l2>
struct PairMapFillerLy1<l1, lx1, 0, l2>
{
	//! Add pair functions, for \a ly1 = 0
	static void fill()
	{
		PairMapFillerLx2<l1, lx1, 0, l2, l2>::fill();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all allowed angular momenta in the \f$x\f$ direction of the first
 * function, holding the total angular momenta fixed.
 */
template <int l1, int lx1, int l2>
struct PairMapFillerLx1
{
	//! Add pair functions, looping over \a lx1
	static void fill()
	{
		PairMapFillerLy1<l1, lx1, l1-lx1, l2>::fill();
		PairMapFillerLx1<l1, lx1-1, l2>::fill();
	}
};
//! Specialization for lx1 = 2
template <int l1, int l2>
struct PairMapFillerLx1<l1, 2, l2>
{
	//! Add pair functions, looping over \a lx1 = 2, 1, 0
	static void fill()
	{
		PairMapFillerLy1<l1, 2, l1-2, l2>::fill();
		PairMapFillerLy1<l1, 1, l1-1, l2>::fill();
		PairMapFillerLy1<l1, 0, l1, l2>::fill();
	}
};
//! Specialization for lx1 = 1
template <int l1, int l2>
struct PairMapFillerLx1<l1, 1, l2>
{
	//! Add pair functions, looping over \a lx1 = 1, 0
	static void fill()
	{
		PairMapFillerLy1<l1, 1, l1-1, l2>::fill();
		PairMapFillerLy1<l1, 0, l1, l2>::fill();
	}
};
//! Specialization for lx1 = 0, base case for recursion
template <int l1, int l2>
struct PairMapFillerLx1<l1, 0, l2>
{
	//! Add pair functions, for \a lx1  = 0
	static void fill()
	{
		PairMapFillerLy1<l1, 0, l1, l2>::fill();
	}
};

/*!
 * \brief Fill the pair functions in the dispatcher
 *
 * Fill the pair functions in the dispatcher. This template is used to loop
 * over all supported total angular momenta of the second function, holding
 * the total angular momentum of the first function fixed.
 */
template <int l1, int l2>
struct PairMapFillerL2
{
	//! Add pair functions, looping over \a l2
	static void fill()
	{
		PairMapFillerLx1<l1, l1, l2>::fill();
		PairMapFillerL2<l1, l2-1>::fill();
	}
};
//! Specialization for highest total supported angular momentum 6
template <int l1>
struct PairMapFillerL2<l1, 6>
{
	//! Add pair functions, looping over \a l2 = 6, 5, ..., 0
	static void fill()
	{
		PairMapFillerLx1<l1, l1, 6>::fill();
		PairMapFillerLx1<l1, l1, 5>::fill();
		PairMapFillerLx1<l1, l1, 4>::fill();
		PairMapFillerLx1<l1, l1, 3>::fill();
		PairMapFillerLx1<l1, l1, 2>::fill();
		PairMapFillerLx1<l1, l1, 1>::fill();
		PairMapFillerLx1<l1, l1, 0>::fill();
	}
};
//! Base case for recursion: total angular momentum 0
template <int l1>
struct PairMapFillerL2<l1, 0>
{
	//! Add pair functions, for l2 = 0
	static void fill()
	{
		PairMapFillerLx1<l1, l1, 0>::fill();
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
		PairMapFillerL2<l1, lmax2>::fill();
		PairMapFiller<l1-1, lmax2>::fill();
	}
};
//! Specialization for highest total supported angular momentum 6
template <int lmax2>
struct PairMapFiller<6, lmax2>
{
	//! Add pair functions, looping over \a l1 = 6, 5, ..., 0
	static void fill()
	{
		PairMapFillerL2<6, lmax2>::fill();
		PairMapFillerL2<5, lmax2>::fill();
		PairMapFillerL2<4, lmax2>::fill();
		PairMapFillerL2<3, lmax2>::fill();
		PairMapFillerL2<2, lmax2>::fill();
		PairMapFillerL2<1, lmax2>::fill();
		PairMapFillerL2<0, lmax2>::fill();
	}
};
//! Base case for recursion: total angular momentum 0
template <int lmax2>
struct PairMapFiller<0, lmax2>
{
	//! Add pair functions, for \a l1 = 0
	static void fill()
	{
		PairMapFillerL2<0, lmax2>::fill();
	}
};

#endif // DISPATCHERFILLER_HH
