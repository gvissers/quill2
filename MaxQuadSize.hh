#ifndef MAXQUADSIZE_HH
#define MAXQUADSIZE_HH

/*!
 * \file MaxQuadSize
 * \brief Set of structures to determine the maximum size of a basis function quartet
 */

#include "CGTOSpecQuad.hh"

template <size_t s1, size_t s2, bool first>
struct Select
{
	static const size_t elem = s1;
};
template <size_t s1, size_t s2>
struct Select<s1, s2, false>
{
	static const size_t elem = s2;
};

template <size_t s1, size_t s2>
struct MaxSize
{
	static const size_t size = Select<s1, s2, s1 >= s2>::elem;
};

template <CGTOShellQuad::PositionSymmetry sym, int l1, int l2,
	int l1max, int l2max, size_t cur_size>
struct MaxQuadSizeIter
{
	static const size_t size = MaxQuadSizeIter<sym, l1, l2-1, l1max, l2max,
		MaxSize<sizeof(CGTOSpecQuad<sym, l1, l2>), cur_size>::size>::size;
};
template <CGTOShellQuad::PositionSymmetry sym, int l1, int l1max, int l2max,
	size_t cur_size>
struct MaxQuadSizeIter<sym, l1, 0, l1max, l2max, cur_size>
{
	static const size_t size = MaxQuadSizeIter<sym, l1-1, l2max,
		l1max, l2max,
		MaxSize<sizeof(CGTOSpecQuad<sym, l1, 0>), cur_size>::size>::size;
};
template <CGTOShellQuad::PositionSymmetry sym, int l1max, int l2max,
	size_t cur_size>
struct MaxQuadSizeIter<sym, 0, 0, l1max, l2max, cur_size>
{
	static const size_t size = MaxQuadSizeIter<
		CGTOShellQuad::PositionSymmetry(sym+1),
		l1max, l2max, l1max, l2max,
		MaxSize<sizeof(CGTOSpecQuad<sym, 0, 0>), cur_size>::size>::size;
};
template <int l1, int l2, int l1max, int l2max, size_t cur_size>
struct MaxQuadSizeIter<CGTOShellQuad::POS_SYM_COUNT, l1, l2, l1max, l2max,
	cur_size>
{
	static const size_t size = cur_size;
};

template <int l1, int l2>
struct MaxQuadSize
{
	static const size_t size = MaxQuadSizeIter<CGTOShellQuad::POS_SYM_ABCD,
		l1, l2, l1, l2, sizeof(CGTOQuad)>::size;
};

#endif // MAXQUADSIZE_HH