#include <memory>
#include "Basis.hh"
#include "Dispatcher.hh"
#include "io/manipulators.hh"

void Basis::twoElectron(const Eigen::MatrixXd& P, Eigen::MatrixXd& G) const
{
	if (!_status.test(ELEC_REP_CURRENT))
		calcElectronRepulsion();

	int n = P.rows();
#ifdef DEBUG
	if (P.cols() != n)
		throw Li::Exception("Density matrix is not square");
#endif
	G.resize(n, n);
	G.setZero();
	int idx = 0;
	double e;
	for (int i = 0; i < n; ++i)
	{
		double Pii = P(i,i);
		for (int j = 0; j < i; ++j)
		{
			double Pij = P(i,j);
			double Pjj = P(j,j);
			for (int k = 0; k < j; ++k)
			{
				double Pik = P(i,k);
				double Pjk = P(j,k);
				for (int l = 0; l < k; ++l)
				{
					e = _elec_rep(idx++);
					G(i,j) += 2 * P(k,l) * e;
					G(i,k) -= 0.5 * P(j,l) * e;
					G(i,l) -= 0.5 * Pjk * e;
					G(j,k) -= 0.5 * P(i,l) * e;
					G(j,l) -= 0.5 * Pik * e;
					G(k,l) += 2 * Pij * e;
				}
				e = _elec_rep(idx++);
				G(i,j) += P(k,k) * e;
				G(i,k) -= 0.5 * Pjk * e;
				G(j,k) -= 0.5 * Pik * e;
				G(k,k) += 2 * Pij * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				G(i,j) += 1.5 * P(j,l) * e;
				G(i,l) -= 0.5 * Pjj * e;
				G(j,j) -= P(i,l) * e;
				G(j,l) += 1.5 * Pij * e;
			}
			e = _elec_rep(idx++);
			G(i,j) += 0.5 * Pjj * e;
			G(j,j) += Pij * e;
			for (int k = j+1; k < i; ++k)
			{
				double Pik = P(i,k);
				double Pkj = P(k,j);
				for (int l = 0; l < j; ++l)
				{
					e = _elec_rep(idx++);
					G(i,j) += 2 * P(k,l) * e;
					G(i,k) -= 0.5 * P(j,l) * e;
					G(i,l) -= 0.5 * Pkj * e;
					G(j,l) -= 0.5 * Pik * e;
					G(k,j) -= 0.5 * P(i,l) * e;
					G(k,l) += 2 * Pij * e;
				}
				e = _elec_rep(idx++);
				G(i,j) += 1.5 * Pkj * e;
				G(i,k) -= 0.5 * Pjj * e;
				G(j,j) -= Pik * e;
				G(k,j) += 1.5 * Pij * e;
				for (int l = j+1; l < k; ++l)
				{
					e = _elec_rep(idx++);
					G(i,j) += 2 * P(k,l) * e;
					G(i,k) -= 0.5 * P(l,j) * e;
					G(i,l) -= 0.5 * Pkj * e;
					G(k,l) += 2 * Pij * e;
					G(k,j) -= 0.5 * P(i,l) * e;
					G(l,j) -= 0.5 * Pik * e;
				}
				e = _elec_rep(idx++);
				G(i,j) += P(k,k) * e;
				G(i,k) -= 0.5 * Pkj * e;
				G(k,j) -= 0.5 * Pik * e;
				G(k,k) += 2 * Pij * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				G(i,i) -= P(j,l) * e;
				G(i,j) += 1.5 * P(i,l) * e;
				G(i,l) += 1.5 * Pij * e;
				G(j,l) -= 0.5 * Pii * e;
			}
			e = _elec_rep(idx++);
			G(i,i) -= 0.5 * Pjj * e;
			G(i,j) += 1.5 * Pij * e;
			G(j,j) -= 0.5 * Pii * e;
		}

		for (int k = 0; k < i; ++k)
		{
			double Pik = P(i, k);
			for (int l = 0; l < k; ++l)
			{
				e = _elec_rep(idx++);
				G(i,i) += 2 * P(k,l) * e;
				G(i,k) -= 0.5 * P(i,l) * e;
				G(i,l) -= 0.5 * Pik * e;
				G(k,l) += Pii * e;
			}
			e = _elec_rep(idx++);
			G(i,i) += P(k,k) * e;
			G(i,k) -= 0.5 * Pik * e;
			G(k,k) += Pii * e;
		}
		for (int l = 0; l < i; ++l)
		{
			e = _elec_rep(idx++);
			G(i,i) += P(i,l) * e;
			G(i,l) += 0.5 * Pii * e;
		}
		G(i,i) += 0.5 * Pii * _elec_rep(idx++);
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j)
			G(j,i) = G(i,j);
	}
}

void Basis::twoElectron(const Eigen::MatrixXd& P, Eigen::MatrixXd& J,
	Eigen::MatrixXd& K) const
{
	if (!_status.test(ELEC_REP_CURRENT))
		calcElectronRepulsion();

	int n = P.rows();
#ifdef DEBUG
	if (P.cols() != n)
		throw Li::Exception("Density matrix is not square");
#endif
	J.resize(n, n);
	J.setZero();
	K.resize(n, n);
	K.setZero();
	int idx = 0;
	double e;
	for (int i = 0; i < n; ++i)
	{
		double Pii = P(i,i);
		for (int j = 0; j < i; ++j)
		{
			double Pij = P(i,j);
			double Pjj = P(j,j);
			for (int k = 0; k < j; ++k)
			{
				double Pik = P(i,k);
				double Pjk = P(j,k);
				for (int l = 0; l < k; ++l)
				{
					e = _elec_rep(idx++);
					J(i,j) += 2 * P(k,l) * e;
					J(k,l) += 2 * Pij * e;
					K(i,k) += P(j,l) * e;
					K(i,l) += Pjk * e;
					K(j,k) += P(i,l) * e;
					K(j,l) += Pik * e;
				}
				e = _elec_rep(idx++);
				J(i,j) += P(k,k) * e;
				J(k,k) += 2 * Pij * e;
				K(i,k) += Pjk * e;
				K(j,k) += Pik * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				J(i,j) += 2 * P(j,l) * e;
				J(j,l) += 2 * Pij * e;
				K(i,j) += P(j,l) * e;
				K(i,l) += Pjj * e;
				K(j,j) += 2 * P(i,l) * e;
				K(j,l) += Pij * e;
			}
			e = _elec_rep(idx++);
			J(i,j) += Pjj * e;
			J(j,j) += 2 * Pij * e;
			K(i,j) += Pjj * e;
			K(j,j) += 2 * Pij * e;
			for (int k = j+1; k < i; ++k)
			{
				double Pik = P(i,k);
				double Pkj = P(k,j);
				for (int l = 0; l < j; ++l)
				{
					e = _elec_rep(idx++);
					J(i,j) += 2 * P(k,l) * e;
					J(k,l) += 2 * Pij * e;
					K(i,k) += P(j,l) * e;
					K(i,l) += Pkj * e;
					K(j,l) += Pik * e;
					K(k,j) += P(i,l) * e;
				}
				e = _elec_rep(idx++);
				J(i,j) += 2 * Pkj * e;
				J(k,j) += 2 * Pij * e;
				K(i,k) += Pjj * e;
				K(i,j) += Pkj * e;
				K(j,j) += 2 * Pik * e;
				K(k,j) += Pij * e;
				for (int l = j+1; l < k; ++l)
				{
					e = _elec_rep(idx++);
					J(i,j) += 2 * P(k,l) * e;
					J(k,l) += 2 * Pij * e;
					K(i,k) += P(l,j) * e;
					K(i,l) += Pkj * e;
					K(k,j) += P(i,l) * e;
					K(l,j) += Pik * e;
				}
				e = _elec_rep(idx++);
				J(i,j) += P(k,k) * e;
				J(k,k) += 2 * Pij * e;
				K(i,k) += Pkj * e;
				K(k,j) += Pik * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				J(i,j) += 2 * P(i,l) * e;
				J(i,l) += 2 * Pij * e;
				K(i,i) += 2 * P(j,l) * e;
				K(i,j) += P(i,l) * e;
				K(i,l) += Pij * e;
				K(j,l) += Pii * e;
			}
			e = _elec_rep(idx++);
			J(i,j) += 2 * Pij * e;
			K(i,i) += Pjj * e;
			K(i,j) += Pij * e;
			K(j,j) += Pii * e;
		}
		for (int k = 0; k < i; ++k)
		{
			double Pik = P(i,k);
			for (int l = 0; l < k; ++l)
			{
				e = _elec_rep(idx++);
				J(i,i) += 2 * P(k,l) * e;
				J(k,l) += Pii * e;
				K(i,k) += P(i,l) * e;
				K(i,l) += Pik * e;
			}
			e = _elec_rep(idx++);
			J(i,i) += P(k,k) * e;
			J(k,k) += Pii * e;
			K(i,k) += Pik * e;
		}
		for (int l = 0; l < i; ++l)
		{
			e = _elec_rep(idx++);
			J(i,i) += 2 * P(i,l) * e;
			J(i,l) += Pii * e;
			K(i,i) += P(i,l) * e;
			K(i,l) += Pii * e;
			K(i,i) += P(i,l) * e;
		}
		e = _elec_rep(idx++);
		J(i,i) += Pii * e;
		K(i,i) += Pii * e;
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			J(j,i) = J(i,j);
			K(j,i) = K(i,j);
		}
	}
}

std::ostream& Basis::print(std::ostream& os) const
{
	os << "Basis (\n" << indent;
	for (BasisFunList::const_iterator it = _funs.begin();
		it != _funs.end(); ++it)
	{
		os << **it << "\n";
	}
	os << dedent << ")";
	return os;
}

void Basis::setPairs() const
{
	const Dispatcher& dispatcher = Dispatcher::singleton();

	_pairs.clear();
	for (BasisFunList::const_iterator iit = _funs.begin(); iit != _funs.end(); ++iit)
	{
		for (BasisFunList::const_iterator jit = _funs.begin(); jit <= iit; ++jit)
		{
			AbstractBFPair *pair = dispatcher.pair(**iit, **jit);
			_pairs.push_back(PairPtr(pair));
		}
	}

	_status.set(PAIRS_CURRENT);
}

void Basis::setQuads() const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	const Dispatcher& dispatcher = Dispatcher::singleton();

	_quads.clear();
	for (PairList::const_iterator iit = _pairs.begin(); iit != _pairs.end(); ++iit)
	{
		for (PairList::const_iterator jit = _pairs.begin(); jit <= iit; ++jit)
		{
			AbstractBFQuad *quad = dispatcher.quad(**iit, **jit, _quad_pool);
			_quads.push_back(QuadPtr(quad, _quad_pool.deleter()));
		}
	}
	
	_status.set(QUADS_CURRENT);
}

void Basis::calcOverlap() const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	PairList::const_iterator pit = _pairs.begin();
	int n = size();
	_overlap.resize(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j, ++pit)
			_overlap(i, j) = _overlap(j, i) = (*pit)->overlap();
		_overlap(i, i) = (*pit)->overlap();
		++pit;
	}

	_status.set(OVERLAP_CURRENT);
}

void Basis::calcKinetic() const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	PairList::const_iterator pit = _pairs.begin();
	int n = size();
	_kinetic.resize(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j, ++pit)
			_kinetic(i, j) = _kinetic(j, i) = (*pit)->kineticEnergy(); 
		_kinetic(i, i) = (*pit)->kineticEnergy();
		++pit;
	}

	_status.set(KINETIC_CURRENT);
}

void Basis::calcOneElectron() const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	PairList::const_iterator pit = _pairs.begin();
	int n = size();
	_overlap.resize(n, n);
	_kinetic.resize(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j, ++pit)
		{
			(*pit)->oneElectron(_overlap(i,j), _kinetic(i,j));
			_overlap(j,i) = _overlap(i,j);
			_kinetic(j,i) = _kinetic(i,j);
		}
		(*pit)->oneElectron(_overlap(i,i), _kinetic(i,i));
		++pit;
	}

	_status.set(OVERLAP_CURRENT);
	_status.set(KINETIC_CURRENT);
}

void Basis::calcNuclearAttraction(const Eigen::MatrixXd& nuc_pos,
	const Eigen::VectorXd& nuc_charge) const
{
	if (!_status.test(PAIRS_CURRENT))
		setPairs();

	PairList::const_iterator pit = _pairs.begin();
	int n = size();
	_nuc_attr.resize(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j, ++pit)
			_nuc_attr(i, j) = _nuc_attr(j, i)
				= (*pit)->nuclearAttraction(nuc_pos, nuc_charge); 
		_nuc_attr(i, i) = (*pit)->nuclearAttraction(nuc_pos, nuc_charge);
		++pit;
	}

	_status.set(NUC_ATTR_CURRENT);
}

void Basis::calcElectronRepulsion() const
{
	if (!_status.test(QUADS_CURRENT))
		setQuads();
	
	_elec_rep.resize(_quads.size());
	for (unsigned int i = 0; i < _quads.size(); ++i)
		_elec_rep[i] = _quads[i]->electronRepulsion();

	_status.set(ELEC_REP_CURRENT);
}
