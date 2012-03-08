#ifndef MAXQUADSIZE_HH
#define MAXQUADSIZE_HH

/*!
 * \file MaxQuadSize
 * \brief Set of structures to determine the maximum size of a basis function quartet
 */

#include "CGTOQuad.hh"

template <int l1, int l2>
struct MaxQuadSize
{
	static const size_t size = sizeof(CGTOQuad);
};

#endif // MAXQUADSIZE_HH