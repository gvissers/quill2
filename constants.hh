#ifndef CONSTANTS_HH
#define CONSTANTS_HH

/*!
 * \file constants.hh
 * \brief Various constants used in calculations
 */

// For most values, see NIST:
// http://physics.nist.gov/cuu/Constants/index.html
namespace Constants
{

const double Ang = 1.0 / 0.52917720859;         //!< Angstrom in a0
const double meter = 1.0e10 * Ang;              //!< meter in a0
const double second = 1.0 / 2.418884326505e-17; //!< second in atomic units
const double kg = 1.0 / 9.10938215e-31;         //!< kilogram in m_e
const double amu = 1.660538782e-27 * kg;        //!< atomic mass unit in m_e
const double Coulomb = 1.0 / 1.602176487e-19;   //!< Coulomb in e
const double c = 299792458 * meter / second;    //!< speed of light
const double Debye = 0.1 * Coulomb * Ang * Ang / (second * c); //!< Debye in e*a0
const double cm1 = 1.0 / 219474.6313705;        //!< cm^-1 in E_h

const double pi_sqrt_pi    = 5.5683279968317078;    //!< \f$\pi^{3/2}\f$
const double sqrt_2_pi_5_4 = 5.9149671727956129;    //!< \f$\sqrt{2}\pi^{5/4}\f$
const double log_pi        = 1.14472988584940017414;//!< \f$\ln\pi\f$
const double log_sqrt_2pi  = 0.91893853320467274178;//!< \f$\ln\sqrt{2\pi}\f$


} // namespace Constants

#endif
