/*
 * Na.hh
 *
 *  Created on: 28.06.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the TiOH4 fluid properties
 */
#ifndef DUMUX_TIOH4_HH
#define DUMUX_TIOH4_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the TiOH4 (Titaniumhydroxide Ion) fluid properties
 */
template <class Scalar>
class TitaniumhydroxideIon
: public Components::Base<Scalar, TitaniumhydroxideIon<Scalar> >
, public Components::Ion<Scalar, TitaniumhydroxideIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for TiOH4.
    */
    static std::string name()
    { return "TiOH4"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of TiOH4.
    */
    static Scalar molarMass()
    { return 115.89e-3; } // kgTiOH4/molTiOH4

   /*!
    * \brief The charge of the TiOH4 ion.
    */
    static constexpr int charge()
    {
        return +4;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

