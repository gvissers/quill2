#ifndef EIGEN_ADDONS_HH
#define EIGEN_ADDONS_HH

const CwiseUnaryOp<internal::scalar_boys_op<Scalar>, const Derived> boys(int m) const
{
	return CwiseUnaryOp<internal::scalar_boys_op<Scalar>, const Derived>(
		derived(), internal::scalar_boys_op<Scalar>(m));
}

#endif // EIGEN_ADDONS_HH