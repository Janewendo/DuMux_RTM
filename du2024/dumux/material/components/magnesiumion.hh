/*
 * Mg.hh
 *
 *  Created on: 2024
 *      Author: Jianwen
 */

/*!
 * \file
 *
 * \brief A class for the Mg2+ fluid properties
 */
#ifndef DUMUX_MG_HH
#define DUMUX_MG_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the Mg (Magnesium Ion) fluid properties
 */
template <class Scalar>
class MagnesiumIon
: public Components::Base<Scalar, MagnesiumIon<Scalar> >
, public Components::Ion<Scalar, MagnesiumIon<Scalar> >
{
public:
   /*!
    * \brief A human readable name for Na.
    */
    static std::string name()
    { return "Mg2+"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Na.
    */
    static Scalar molarMass()
    { return 24.305e-3; } // kgNa/molNa

   /*!
    * \brief The charge of the Na ion.
    */
    static constexpr int charge()
    {
        return +2;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

