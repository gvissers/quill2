#ifndef QUILL_FUNCTORS_HH
#define QUILL_FUNCTORS_HH

template <typename Scalar>
struct scalar_qexp_op
{
	scalar_qexp_op() {}
	const Scalar operator()(const Scalar& a) const
	{
		return qexp(a);
	}
#ifdef __SSE2__
	template <typename Packet>
	const Packet packetOp(const Packet& a) const
	{
		return qexp(a);
	}
#endif
};
template <typename Scalar>
struct functor_traits< scalar_qexp_op<Scalar> >
{
	enum
	{
		Cost = 15 * NumTraits<Scalar>::MulCost,
#ifdef __SSE2__
		PacketAccess = true
#else
		PacketAccess = false
#endif
	};
};

template <typename Scalar>
struct scalar_qerf_op
{
	scalar_qerf_op() {}
	const Scalar operator()(const Scalar& a) const
	{
		return qerf(a);
	}
#ifdef __SSE2__
	template <typename Packet>
	const Packet packetOp(const Packet& a) const
	{
		return qerf(a);
	}
#endif
};
template <typename Scalar>
struct functor_traits< scalar_qerf_op<Scalar> >
{
	enum
	{
		Cost = 50 * NumTraits<Scalar>::MulCost,
#ifdef __SSE2__
		PacketAccess = true
#else
		PacketAccess = false
#endif
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
#ifdef __SSE2__
	template <typename Packet>
	const Packet packetOp(const Packet& a, const Packet& expma) const
	{
		return Fm(_m, a, expma);
	}
#endif
	int _m;
};
template <typename Scalar>
struct functor_traits< scalar_boys_with_exp_op<Scalar> >
{
	enum
	{
		Cost = 50 * NumTraits<Scalar>::MulCost,
#ifdef __SSE2__
		PacketAccess = true
#else
		PacketAccess = false
#endif
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
	const Scalar operator()(const Scalar& a) const
	{
		return Fm(_m, a);
	}
#ifdef __SSE2__
	template <typename Packet>
	const Packet packetOp(const Packet& a) const
	{
		return Fm(_m, a);
	}
#endif
	int _m;
};
template <typename Scalar>
struct functor_traits< scalar_boys_op<Scalar> >
{
	enum
	{
		Cost = functor_traits < scalar_boys_with_exp_op<Scalar> >::Cost
			+ functor_traits< scalar_exp_op<Scalar> >::Cost,
#ifdef __SSE2__
		PacketAccess = true
#else
		PacketAccess = false
#endif
	};
};

#endif // QUILL_FUNCTORS_HH