#ifndef QUILL_FUNCTORS_HH
#define QUILL_FUNCTORS_HH

/*!
 * \brief Functor for computing the Boys functions
 *
 * Functor to compute the Boys function
 * \f{eqnarray*}
 *    F_m(t) &=& \int_0^1 \exp(-t s^2) s^{2m} ds\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) P(m+1/2, t)\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) [1-Q(m+1/2, t)]
 * \f}
 * for the values \f$t\f$ in the calling array. This functor is used as a binary
 * operation, with the second argument being \f$e^-t\f$, which can save
 * recomputing the exponential if it is already known.
 */
template <typename Scalar>
struct scalar_boys_with_exp_op
{
	scalar_boys_with_exp_op(int m): _m(m) {}
	const Scalar operator()(const Scalar& a, const Scalar& expma) const
	{
		return Fm(_m, a, expma);
	}
	int _m;
};
template <typename Scalar>
struct functor_traits< scalar_boys_with_exp_op<Scalar> >
{
	enum
	{
		Cost = 50 * NumTraits<Scalar>::MulCost,
		PacketAccess = false
	};
};

/*!
 * \brief Functor for computing the Boys functions
 *
 * Functor to compute the Boys function
 * \f{eqnarray*}
 *    F_m(t) &=& \int_0^1 \exp(-t s^2) s^{2m} ds\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) P(m+1/2, t)\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) [1-Q(m+1/2, t)]
 * \f}
 * for the values \f$t\f$ in the calling array.
 */
template <typename Scalar>
struct scalar_boys_op
{
	scalar_boys_op(int m): _m(m) {}
	const Scalar operator()(const Scalar& a) const { return Fm(_m, a); }
	int _m;
};
template <typename Scalar>
struct functor_traits< scalar_boys_op<Scalar> >
{
	enum
	{
		Cost = functor_traits < scalar_boys_with_exp_op<Scalar> >::Cost
			+ functor_traits< scalar_exp_op<Scalar> >::Cost,
		PacketAccess = false
	};
};

#endif // QUILL_FUNCTORS_HH