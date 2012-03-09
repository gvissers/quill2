#ifndef EIGEN_ADDONS_HH
#define EIGEN_ADDONS_HH

const CwiseUnaryOp<internal::scalar_qexp_op<Scalar>, const Derived> qexp() const
{
	return derived();
}

const CwiseUnaryOp<internal::scalar_boys_op<Scalar>, const Derived>
boys(int m) const
{
	return CwiseUnaryOp<internal::scalar_boys_op<Scalar>, const Derived>(
		derived(), internal::scalar_boys_op<Scalar>(m));
}
template <typename OtherDerived>
const CwiseBinaryOp<internal::scalar_boys_with_exp_op<Scalar>, const Derived, const OtherDerived>
boys(int m, const ArrayBase<OtherDerived>& other) const
{
	return CwiseBinaryOp<internal::scalar_boys_with_exp_op<Scalar>, const Derived, const OtherDerived>(
		derived(), other.derived(), internal::scalar_boys_with_exp_op<Scalar>(m));
}

#endif // EIGEN_ADDONS_HH