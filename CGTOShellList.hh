#ifndef CGTOSHELLLIST_HH
#define CGTOSHELLLIST_HH

/*!
 * \file CGTOShellList.hh
 * \brief Definition of the CGTOShellList class
 */

#include <memory>
#include <Singleton.hh>
#include "CGTOShellQuad.hh"

class CGTOShellList: public Li::Singleton<CGTOShellList>
{
public:
	typedef std::unique_ptr<CGTOShell> ShellPtr;
	typedef std::vector<ShellPtr> ShellList;
	typedef std::unique_ptr<CGTOShellPair> ShellPairPtr;
	typedef std::vector<ShellPairPtr> ShellPairList;
	typedef std::unique_ptr<CGTOShellQuad> ShellQuadPtr;
	typedef std::vector<ShellQuadPtr> ShellQuadList;

	CGTOShellList(): Li::Singleton<CGTOShellList>(),
		_shells(), _pairs(), _quads() {}

	int nrShells() const { return _shells.size(); }
	int nrPairs() const { return _pairs.size(); }
	int nrQuads() const { return _quads.size(); }
	
	int addShell(int l, const Eigen::ArrayXd& weights,
		const Eigen::ArrayXd& widths, int ipos,
		const Eigen::Vector3d& center);
	const CGTOShell& shell(int i)
	{
		return *_shells[i];
	}

	const CGTOShellPair& pair(int i, int j)
	{
		return pair(pairIndex(i, j));
	}
	const CGTOShellPair& pair(int ij)
	{
		return *_pairs[ij];
	}

	const CGTOShellQuad& quad(int i, int j, int k, int l)
	{
		return quad(pairIndex(pairIndex(i, j), pairIndex(k, l)));
	}
	const CGTOShellQuad& quad(int ijkl)
	{
		return *_quads[ijkl];
	}

	static int pairIndex(int i, int j)
	{
		if (i < j)
			std::swap(i, j);
		return i*(i+1)/2 + j;
	}

private:
	ShellList _shells;
	ShellPairList _pairs;
	ShellQuadList _quads;

	int addPair(const CGTOShell& shA, const CGTOShell& shB);
	int addQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD);
};

#endif // CGTOSHELLIST_HH