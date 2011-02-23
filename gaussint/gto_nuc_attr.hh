#ifndef GTO_NUC_ATTR_HH
#define GTO_NUC_ATTR_HH

double gto_nuc_attr_generic(const Eigen::Vector3i& ls1,
	const Eigen::VectorXd& weights1, const Eigen::VectorXd& widths1,
	const Eigen::Vector3d& pos1,
	const Eigen::Vector3i& ls2,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& pos2,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge);

#endif // GTO_NUC_ATTR_HH
