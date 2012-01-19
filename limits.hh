#ifndef LIMITS_HH
#define LIMITS_HH

/*!
 * \file limits.hh
 * \brief Various limits used in quill
 */

// defines for the limits, used for conditional compilation
#define LMAX_SPECIALIZED 2

namespace Limits
{

//! Highest supported total angular momentum
const int lmax = 6;
//! Highest total angular momentum for which some functions are specialized
const int lmax_specialized = LMAX_SPECIALIZED;

} // namespace

#endif // LIMITS_HH
