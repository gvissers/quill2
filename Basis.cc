#include <memory>
#include <Eigen/MatrixFunctions>
#include "Basis.hh"
#include "Dispatcher.hh"
#include "io/manipulators.hh"

void Basis::twoElectron(const Eigen::MatrixXd& D, Eigen::MatrixXd& G) const
{
	if (!_status.test(ELEC_REP_CURRENT))
		calcElectronRepulsion();

	int n = D.rows();
#ifdef DEBUG
	if (D.cols() != n)
		throw Li::Exception("Density matrix is not square");
#endif
	G.resize(n, n);
	G.setZero();
	int idx = 0;
	double e;
	for (int i = 0; i < n; ++i)
	{
		double Dii = D(i,i);
		for (int j = 0; j < i; ++j)
		{
			double Dij = D(i,j);
			double Djj = D(j,j);
			for (int k = 0; k < j; ++k)
			{
				double Dik = D(i,k);
				double Djk = D(j,k);
				for (int l = 0; l < k; ++l)
				{
					e = _elec_rep(idx++);
					G(i,j) += 2 * D(k,l) * e;
					G(i,k) -= 0.5 * D(j,l) * e;
					G(i,l) -= 0.5 * Djk * e;
					G(j,k) -= 0.5 * D(i,l) * e;
					G(j,l) -= 0.5 * Dik * e;
					G(k,l) += 2 * Dij * e;
				}
				e = _elec_rep(idx++);
				G(i,j) += D(k,k) * e;
				G(i,k) -= 0.5 * Djk * e;
				G(j,k) -= 0.5 * Dik * e;
				G(k,k) += 2 * Dij * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				G(i,j) += 1.5 * D(j,l) * e;
				G(i,l) -= 0.5 * Djj * e;
				G(j,j) -= D(i,l) * e;
				G(j,l) += 1.5 * Dij * e;
			}
			e = _elec_rep(idx++);
			G(i,j) += 0.5 * Djj * e;
			G(j,j) += Dij * e;
			for (int k = j+1; k < i; ++k)
			{
				double Dik = D(i,k);
				double Dkj = D(k,j);
				for (int l = 0; l < j; ++l)
				{
					e = _elec_rep(idx++);
					G(i,j) += 2 * D(k,l) * e;
					G(i,k) -= 0.5 * D(j,l) * e;
					G(i,l) -= 0.5 * Dkj * e;
					G(j,l) -= 0.5 * Dik * e;
					G(k,j) -= 0.5 * D(i,l) * e;
					G(k,l) += 2 * Dij * e;
				}
				e = _elec_rep(idx++);
				G(i,j) += 1.5 * Dkj * e;
				G(i,k) -= 0.5 * Djj * e;
				G(j,j) -= Dik * e;
				G(k,j) += 1.5 * Dij * e;
				for (int l = j+1; l < k; ++l)
				{
					e = _elec_rep(idx++);
					G(i,j) += 2 * D(k,l) * e;
					G(i,k) -= 0.5 * D(l,j) * e;
					G(i,l) -= 0.5 * Dkj * e;
					G(k,l) += 2 * Dij * e;
					G(k,j) -= 0.5 * D(i,l) * e;
					G(l,j) -= 0.5 * Dik * e;
				}
				e = _elec_rep(idx++);
				G(i,j) += D(k,k) * e;
				G(i,k) -= 0.5 * Dkj * e;
				G(k,j) -= 0.5 * Dik * e;
				G(k,k) += 2 * Dij * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				G(i,i) -= D(j,l) * e;
				G(i,j) += 1.5 * D(i,l) * e;
				G(i,l) += 1.5 * Dij * e;
				G(j,l) -= 0.5 * Dii * e;
			}
			e = _elec_rep(idx++);
			G(i,i) -= 0.5 * Djj * e;
			G(i,j) += 1.5 * Dij * e;
			G(j,j) -= 0.5 * Dii * e;
		}

		for (int k = 0; k < i; ++k)
		{
			double Dik = D(i, k);
			for (int l = 0; l < k; ++l)
			{
				e = _elec_rep(idx++);
				G(i,i) += 2 * D(k,l) * e;
				G(i,k) -= 0.5 * D(i,l) * e;
				G(i,l) -= 0.5 * Dik * e;
				G(k,l) += Dii * e;
			}
			e = _elec_rep(idx++);
			G(i,i) += D(k,k) * e;
			G(i,k) -= 0.5 * Dik * e;
			G(k,k) += Dii * e;
		}
		for (int l = 0; l < i; ++l)
		{
			e = _elec_rep(idx++);
			G(i,i) += D(i,l) * e;
			G(i,l) += 0.5 * Dii * e;
		}
		G(i,i) += 0.5 * Dii * _elec_rep(idx++);
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < i; ++j)
			G(j,i) = G(i,j);
	}
}

void Basis::twoElectron(const Eigen::MatrixXd& D, Eigen::MatrixXd& J,
	Eigen::MatrixXd& K) const
{
	if (!_status.test(ELEC_REP_CURRENT))
		calcElectronRepulsion();

	int n = D.rows();
#ifdef DEBUG
	if (D.cols() != n)
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
		double Dii = D(i,i);
		for (int j = 0; j < i; ++j)
		{
			double Dij = D(i,j);
			double Djj = D(j,j);
			for (int k = 0; k < j; ++k)
			{
				double Dik = D(i,k);
				double Djk = D(j,k);
				for (int l = 0; l < k; ++l)
				{
					e = _elec_rep(idx++);
					J(i,j) += 2 * D(k,l) * e;
					J(k,l) += 2 * Dij * e;
					K(i,k) += D(j,l) * e;
					K(i,l) += Djk * e;
					K(j,k) += D(i,l) * e;
					K(j,l) += Dik * e;
				}
				e = _elec_rep(idx++);
				J(i,j) += D(k,k) * e;
				J(k,k) += 2 * Dij * e;
				K(i,k) += Djk * e;
				K(j,k) += Dik * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				J(i,j) += 2 * D(j,l) * e;
				J(j,l) += 2 * Dij * e;
				K(i,j) += D(j,l) * e;
				K(i,l) += Djj * e;
				K(j,j) += 2 * D(i,l) * e;
				K(j,l) += Dij * e;
			}
			e = _elec_rep(idx++);
			J(i,j) += Djj * e;
			J(j,j) += 2 * Dij * e;
			K(i,j) += Djj * e;
			K(j,j) += 2 * Dij * e;
			for (int k = j+1; k < i; ++k)
			{
				double Dik = D(i,k);
				double Dkj = D(k,j);
				for (int l = 0; l < j; ++l)
				{
					e = _elec_rep(idx++);
					J(i,j) += 2 * D(k,l) * e;
					J(k,l) += 2 * Dij * e;
					K(i,k) += D(j,l) * e;
					K(i,l) += Dkj * e;
					K(j,l) += Dik * e;
					K(k,j) += D(i,l) * e;
				}
				e = _elec_rep(idx++);
				J(i,j) += 2 * Dkj * e;
				J(k,j) += 2 * Dij * e;
				K(i,k) += Djj * e;
				K(i,j) += Dkj * e;
				K(j,j) += 2 * Dik * e;
				K(k,j) += Dij * e;
				for (int l = j+1; l < k; ++l)
				{
					e = _elec_rep(idx++);
					J(i,j) += 2 * D(k,l) * e;
					J(k,l) += 2 * Dij * e;
					K(i,k) += D(l,j) * e;
					K(i,l) += Dkj * e;
					K(k,j) += D(i,l) * e;
					K(l,j) += Dik * e;
				}
				e = _elec_rep(idx++);
				J(i,j) += D(k,k) * e;
				J(k,k) += 2 * Dij * e;
				K(i,k) += Dkj * e;
				K(k,j) += Dik * e;
			}
			for (int l = 0; l < j; ++l)
			{
				e = _elec_rep(idx++);
				J(i,j) += 2 * D(i,l) * e;
				J(i,l) += 2 * Dij * e;
				K(i,i) += 2 * D(j,l) * e;
				K(i,j) += D(i,l) * e;
				K(i,l) += Dij * e;
				K(j,l) += Dii * e;
			}
			e = _elec_rep(idx++);
			J(i,j) += 2 * Dij * e;
			K(i,i) += Djj * e;
			K(i,j) += Dij * e;
			K(j,j) += Dii * e;
		}
		for (int k = 0; k < i; ++k)
		{
			double Dik = D(i,k);
			for (int l = 0; l < k; ++l)
			{
				e = _elec_rep(idx++);
				J(i,i) += 2 * D(k,l) * e;
				J(k,l) += Dii * e;
				K(i,k) += D(i,l) * e;
				K(i,l) += Dik * e;
			}
			e = _elec_rep(idx++);
			J(i,i) += D(k,k) * e;
			J(k,k) += Dii * e;
			K(i,k) += Dik * e;
		}
		for (int l = 0; l < i; ++l)
		{
			e = _elec_rep(idx++);
			J(i,i) += 2 * D(i,l) * e;
			J(i,l) += Dii * e;
			K(i,i) += D(i,l) * e;
			K(i,l) += Dii * e;
			K(i,i) += D(i,l) * e;
		}
		e = _elec_rep(idx++);
		J(i,i) += Dii * e;
		K(i,i) += Dii * e;
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

void Basis::calcOrtho() const
{
	if (!_status.test(OVERLAP_CURRENT))
		calcOverlap();
	_ortho = _overlap.sqrt().inverse();
	_status.set(ORTHO_CURRENT);
}
