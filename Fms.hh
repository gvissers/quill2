#ifndef FMS_HH
#define FMS_HH

#include "boys.hh"

struct Fms: public MultiArray
{
	Fms(int maxm, const Eigen::ArrayXXd& T):
		MultiArray(T.rows(), T.cols(), maxm+1)
	{
		Eigen::ArrayXXd expmT = (-T).qexp();
		(*this)[maxm] = T.boys(maxm, expmT);
		for (int m = maxm-1; m >= 0; --m)
			(*this)[m] = (expmT + 2*T*(*this)[m+1]) / (2*m+1);
	}

	Fms(int maxm, const Eigen::ArrayXXd& T, const Eigen::ArrayXXd& KKW):
		MultiArray(T.rows(), T.cols(), maxm+1)
	{
		Eigen::ArrayXXd expmT = (-T).qexp();
		(*this)[maxm] = T.boys(maxm, expmT);
		for (int m = maxm-1; m >= 0; --m)
			(*this)[m] = (expmT + 2*T*(*this)[m+1]) / (2*m+1);
		for (int m = maxm; m >= 0; --m)
			(*this)[m] *= KKW;
	}
};

#endif // FMS_HH