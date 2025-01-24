/*
 * Na.hh
 *
 *  Created on: 28.06.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the K+ fluid properties
 */
#ifndef DUMUX_K_HH
#define DUMUX_K_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the K (Potassium Ion) fluid properties
 */
template <class Scalar>
class PotassiumIon
: public Components::Base<Scalar, PotassiumIon<Scalar> >
, public Components::Ion<Scalar, PotassiumIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for K.
    */
    static std::string name()
    { return "K+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of K.
    */
    static Scalar molarMass()
    { return 39.0983e-3; } // kgNa/molNa

   /*!
    * \brief The charge of the K ion.
    */
    static constexpr int charge()
    {
        return +1;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

