#include <memory>
#include "Basis.hh"
#include "Dispatcher.hh"
#include "io/manipulators.hh"

Eigen::MatrixXd Basis::electronRepulsion(const Eigen::MatrixXd& P)
{
	if (!_status.test(ELEC_REP_CURRENT))
		calcElectronRepulsion();

	int n = P.rows();
#ifdef DEBUG
	if (P.cols() != n)
		throw Li::Exception("Density matrix is not square");
#endif
	Eigen::MatrixXd G(n, n);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j <= i; ++j)
		{
			double gij = 0;
			for (int k = 0; k < n; ++k)
			{
				for (int l = 0; l < n; ++l)
				{
					gij += P(k,l)
						* (eri(i, k, j, l) - 0.5*eri(i, k, l, j));
				}
			}
			G(i,j) = G(j,i) = gij;
		}
	}

	return G;
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
			_pairs.push_back(dispatcher.pair(**iit, **jit));
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
			_quads.push_back(dispatcher.quad(**iit, **jit));
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

#include <CGTOQuad.hh>
void Basis::calcElectronRepulsion() const
{
	if (!_status.test(QUADS_CURRENT))
		setQuads();
	
	_elec_rep.resize(_quads.size());
	for (unsigned int i = 0; i < _quads.size(); ++i)
	{
		_elec_rep[i] = _quads[i]->electronRepulsion();
/*
		if (isnan(_elec_rep[i]))
		{
			CGTOQuad *qd = static_cast<CGTOQuad*>(_quads[i].get());
			std::cout << qd->p().f().ls().transpose() << ", "
				<< qd->p().g().ls().transpose() << ", "
				<< qd->q().f().ls().transpose() << ", "
				<< qd->q().g().ls().transpose() << "\n";
			std::cout << qd->p().f().center().transpose() << ", "
				<< qd->p().g().center().transpose() << ", "
				<< qd->q().f().center().transpose() << ", "
				<< qd->q().g().center().transpose() << "\n";
		}
*/
		std::cout << i << " " << _elec_rep[i] << "\n";
	}
	
	_status.set(ELEC_REP_CURRENT);
}
