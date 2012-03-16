#ifndef ANGMOMITERATOR_HH
#define ANGMOMITERATOR_HH

/*!
 * \file AngMomIterator.hh
 * \brief Definition of the AngMomIterator class
 */

/*!
 * \brief Class for iterating of possible angular momenta
 *
 * Class AngMomIterator can be used to iterate over all possible distribution
 * of a given total angular momentum over the three dimensions.
 */
class AngMomIterator
{
public:
	/*!
	 * \brief Constructor
	 *
	 * Create a new AngMomIterator object for total angular momentum \a ltot.
	 */
	AngMomIterator(int ltot) { reset(ltot); }

	/*!
	 * \brief Test for end of iteration
	 *
	 * Check if we have iterated over all possible distributions, returning
	 * \c true if we did and \c false otherwise.
	 */
	bool end() const { return lx() < 0; }

	//! Return the current angular momentum in the \f$x\f$ dimension
	int lx() const { return _lx; }
	//! Return the current angular momentum in the \f$y\f$ dimension
	int ly() const { return _ly; }
	//! Return the current angular momentum in the \f$x\f$ dimension
	int lz() const { return _lz; }
	//! Return the total angular momentum we are iterating over.
	int ltot() const { return _ltot; }

	//! Reset the iterator, starting anew with the current total angular momentum
	void reset() { reset(_ltot); }
	//! Reset the iterator, starting anew with the total angular momentum \a ltot.
	void reset(int ltot)
	{
		_ltot = ltot;
		_lx = ltot;
		_ly = 0;
		_lz = 0;
	}

	/*!
	 * \brief Advance the iterator
	 *
	 * Advance the iterator to the next distribution.
	 * \return The updated iterator.
	 */
	AngMomIterator& operator++()
	{
		if (_ly > 0)
		{
			--_ly;
			++_lz;
		}
		else
		{
			--_lx;
			_ly = _lz+1;
			_lz = 0;
		}
		return *this;
	}

private:
	//! The total angular momentum iterated over
	int _ltot;
	//! The angular momentum in the \f$x\f$ dimension
	int _lx;
	//! The angular momentum in the \f$y\f$ dimension
	int _ly;
	//! The angular momentum in the \f$z\f$ dimension
	int _lz;
};

#endif // ANGMOMITERATOR_HH
