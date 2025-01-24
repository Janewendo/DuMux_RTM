// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Components
 * \brief A class for the basalt glass phase properties
 */
#ifndef DUMUX_BIRNESSITE_HH
#define DUMUX_BIRNESSITE_HH

#include <dumux/material/components/base.hh>
#include <dumux/material/components/solid.hh>

#include <dumux/material/components/magnesiumion.hh>
#include <dumux/material/components/aluminumion.hh>
#include <dumux/material/components/potassiumion.hh>
#include <dumux/material/components/manganeseion.hh>
#include <dumux/material/components/titaniumhydroxideion.hh>
#include <dumux/material/components/hydrogenphosphateion.hh>
#include <dumux/material/components/silicondioxide.hh>
#include <dumux/material/components/sodiumion.hh>
#include <dumux/material/components/calciumion.hh>
#include <dumux/material/components/iron2ion.hh>
#include <dumux/material/components/hydronion.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/o2.hh>

namespace Dumux {
namespace Components {

/*!
 * \ingroup Components
 * \brief A class for the basalt glass phase properties
 */
template <class Scalar>
class Birnessite
: public Components::Base<Scalar, Birnessite<Scalar> >
, public Components::Solid<Scalar, Birnessite<Scalar> >
{

public:
	using H2O = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using Na = Components::SodiumIon<Scalar>;
    using Ca = Components::CalciumIon<Scalar>;
    using Mg = Components::MagnesiumIon<Scalar>;
    using Al = Components::AluminumIon<Scalar>;
    using SiO2 = Components::Silicondioxide<Scalar>;
    using HPO4 = Components::HydrogenphosphateIon<Scalar>;
    using K = Components::PotassiumIon<Scalar>;
    using Mn = Components::ManganeseIon<Scalar>;
	using TiOH4 = Components::TitaniumhydroxideIon<Scalar>;
    using H = Components::HydronIon<Scalar>;
    using Fe2 = Components::Iron2Ion<Scalar>;
    using O2 = Components::O2<Scalar>;
    /*!
     * \brief A human readable name for Ferrohydrite.
     */
    static std::string name()
    { return "Birnessite"; }

    /*!
     * \brief The molar mass in \f$\mathrm{[kg/mol]}\f$ of Ferrohydrite.
     */
    static constexpr Scalar molarMass()
    { return 8*Mn::molarMass()+ 19*O2::molarMass() + 10*H::molarMass(); } // kg/mol

    /*!
     * \brief Returns true if the solid phase is assumed to be compressible
     */
    static constexpr bool solidIsCompressible()
    { return false; }

    /*!
     * \brief The density in \f$\mathrm{[kg/m^3]}\f$ of the component at a given pressure in
     *          \f$\mathrm{[Pa]}\f$ and temperature in \f$\mathrm{[K]}\f$.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static constexpr Scalar solidDensity(Scalar temperature)
    { return 5e3; } //don't know

    /*!
     * \brief Thermal conductivity of the component \f$\mathrm{[W/(m*K)]}\f$ as a solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidThermalConductivity(Scalar temperature)
    { return 0.5; } //don't know

    /*!
     * \brief Specific isobaric heat capacity of the component \f$\mathrm{[J/(kg*K)]}\f$ as a solid.
     *
     * \param temperature temperature of component in \f$\mathrm{[K]}\f$
     */
    static Scalar solidHeatCapacity(Scalar temperature)
    { return 420; } //don't know
};

} // end namespace Components
} // end namespace Dumux

#endif

