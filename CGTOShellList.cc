#include "CGTOShellList.hh"

SINGLETON_OBJECT(CGTOShellList);

const CGTOShell& CGTOShellList::addShell(int l, const Eigen::ArrayXd& weights,
	const Eigen::ArrayXd& widths, int ipos,
	const Eigen::Vector3d& center)
{
	CGTOShell *sh = new CGTOShell(nrShells(), l, weights, widths,
		ipos, center);
	for (int i = 0; i < nrShells(); ++i)
		addPair(*sh, shell(i));
	addPair(*sh, *sh);
	_shells.push_back(ShellPtr(sh));
	return *sh;
}

int CGTOShellList::addPair(const CGTOShell& shA, const CGTOShell& shB)
{
	CGTOShellPair *p = new CGTOShellPair(shA, shB);
	for (int ij = 0; ij < nrPairs(); ++ij)
		addQuad(*p, pair(ij));
	addQuad(*p, *p);
	_pairs.push_back(ShellPairPtr(p));
	return nrPairs()-1;
}

int CGTOShellList::addQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD)
{
	CGTOShellQuad *q = new CGTOShellQuad(pAB, pCD);
	_quads.push_back(ShellQuadPtr(q));
	return nrQuads()-1;
}