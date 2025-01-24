/*
 * Na.hh
 *
 *  Created on: 28.06.2011
 *      Author: kissinger
 */

/*!
 * \file
 *
 * \brief A class for the SiO2 fluid properties
 */
#ifndef DUMUX_SIO2_HH
#define DUMUX_SIO2_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the SiO2 (Silicondioxide Ion) fluid properties
 */
template <class Scalar>
class Silicondioxide
: public Components::Base<Scalar, Silicondioxide<Scalar> >
, public Components::Ion<Scalar, Silicondioxide<Scalar> >
{
public:
   /*!
    * \brief A human readable name for SiO2.
    */
    static std::string name()
    { return "SiO2"; }

   /*!
    * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of SiO2.
    */
    static Scalar molarMass()
    { return 60.08e-3; } // kgSiO2/molSiO2

   /*!
    * \brief The charge of the SiO2 ion.
    */
    static constexpr int charge()
    {
        return 0;
    }

};

} // end namespace Components
} // end namespace Dumux


#endif

