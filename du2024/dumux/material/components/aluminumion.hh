/*
 * Al.hh
 *
 *  Created on: 2024
 *      Author: Jianwen
 */

/*!
 * \file
 *
 * \brief A class for the Al3+ fluid properties
 */
#ifndef DUMUX_AL_HH
#define DUMUX_AL_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Al (Aluminum Ion) fluid properties
 */
template <class Scalar>
class AluminumIon
: public Components::Base<Scalar, AluminumIon<Scalar> >
, public Components::Ion<Scalar, AluminumIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Al.
    */
    static std::string name()
    { return "Al3+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Al.
    */
    static Scalar molarMass()
    { return 26.9815e-3; } // kgNa/molAl

   /*!
    * \brief The charge of the Al ion.
    */
    static constexpr int charge()
    {
        return +3;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

