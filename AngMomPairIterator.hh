#ifndef ANGMOMPAIRITERATOR_HH
#define ANGMOMPAIRITERATOR_HH

/*!
 * \file AngMomPairIterator.hh
 * \brief Definition of the AngMomPairIterator class
 */

#include "AngMomIterator.hh"

/*!
 * \brief Iterator over two angular momenta
 *
 * Class AngMomPairIterator can be used to iterate over all possible
 * distributions of a given total angular momentum over two sets of angular
 * momentum quantum numbers.
 */
class AngMomPairIterator
{
public:
	/*!
	 * \brief Constructor
	 *
	 * Create a new AngMomPairIterator object, which distributes total
	 * angular momentum \a ltot over the two pairs of quantum numbers,
	 * where the total angular momentum of the first set will not exceed
	 * \a lmax1 and that of the second set is at most \a lmax2.
	 */
	AngMomPairIterator(int ltot, int lmax1, int lmax2):
		_ltot(ltot), _lmax1(lmax1), _lmax2(lmax2),
		_it1(std::min(ltot, lmax1)), _it2(ltot-l1()) {}

	/*!
	 * \brief Test for end of iteration
	 *
	 * Check if we have iterated over all possible distributions, returning
	 * \c true if we did and \c false otherwise.
	 */
	bool end() const { return l1() < 0 || l2() > _lmax2; }

	//! Return the total angular momentum iterated over
	int ltot() const { return _ltot; }
	
	//! Return the current total angular momentum of the first set
	int l1() const { return _it1.ltot(); }
	//! Return the current total angular momentum of the second set
	int l2() const { return _it2.ltot(); }

	//! Return the angular momentum in the \f$x\f$ dimension of the first set
	int lx1() const { return _it1.lx(); }
	//! Return the angular momentum in the \f$y\f$ dimension of the first set
	int ly1() const { return _it1.ly(); }
	//! Return the angular momentum in the \f$z\f$ dimension of the first set
	int lz1() const { return _it1.lz(); }
	//! Return the angular momentum in the \f$x\f$ dimension of the second set
	int lx2() const { return _it2.lx(); }
	//! Return the angular momentum in the \f$y\f$ dimension of the second set
	int ly2() const { return _it2.ly(); }
	//! Return the angular momentum in the \f$z\f$ dimension of the second set
	int lz2() const { return _it2.lz(); }

	//! Return the total angular momentum in the \f$x\f$ dimension
	int lx() const { return lx1() + lx2(); }
	//! Return the total angular momentum in the \f$y\f$ dimension
	int ly() const { return ly1() + ly2(); }
	//! Return the total angular momentum in the \f$z\f$ dimension
	int lz() const { return lz1() + lz2(); }

	/*!
	 * \brief Advance the iterator
	 *
	 * Advance the iterator to the next distribution.
	 * \return The updated iterator.
	 */
	AngMomPairIterator& operator++()
	{
		if ((++_it2).end())
		{
			if ((++_it1).end())
			{
				_it1.reset(l1()-1);
				_it2.reset(l2()+1);
			}
			else
			{
				_it2.reset();
			}
		}
		return *this;
	}

private:
	//! The total angular momentum iterated over
	int _ltot;
	//! The maximum total angular momentum of the first set
	int _lmax1;
	//! The maximum total angular momentum of the second set
	int _lmax2;
	//! Iterator over the first set
	AngMomIterator _it1;
	//! Iterator over the second set
	AngMomIterator _it2;
};

#endif // ANGMOMPAIRITERATOR_HH