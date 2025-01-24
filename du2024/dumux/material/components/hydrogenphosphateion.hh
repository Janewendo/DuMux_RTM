/*
 * Na.hh
 *
 *  Created on: 28.06.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the HPO4 fluid properties
 */
#ifndef DUMUX_HPO4_HH
#define DUMUX_HPO4_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the HPO4 (Hydrogenphosphate Ion) fluid properties
 */
template <class Scalar>
class HydrogenphosphateIon
: public Components::Base<Scalar, HydrogenphosphateIon<Scalar> >
, public Components::Ion<Scalar, HydrogenphosphateIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for HPO4.
    */
    static std::string name()
    { return "HPO4"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of HPO4.
    */
    static Scalar molarMass()
    { return 96.99e-3; } // kgHPO4/molHPO4

   /*!
    * \brief The charge of the HPO4 ion.
    */
    static constexpr int charge()
    {
        return -2;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

