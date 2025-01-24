/*
 * Na.hh
 *
 *  Created on: 28.06.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the Mn2+ fluid properties
 */
#ifndef DUMUX_MN_HH
#define DUMUX_MN_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Mn (Manganese Ion) fluid properties
 */
template <class Scalar>
class ManganeseIon
: public Components::Base<Scalar, ManganeseIon<Scalar> >
, public Components::Ion<Scalar, ManganeseIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Na.
    */
    static std::string name()
    { return "Mn2+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Mn.
    */
    static Scalar molarMass()
    { return 54.9380e-3; } // kgMn/molMn

   /*!
    * \brief The charge of the Mn ion.
    */
    static constexpr int charge()
    {
        return +2;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

