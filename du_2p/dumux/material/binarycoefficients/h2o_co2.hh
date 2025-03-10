// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and oxygen.
 */
#ifndef DUMUX_BINARY_COEFF_H2O_CO2_HH
#define DUMUX_BINARY_COEFF_H2O_CO2_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/components/co2tablereader.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/h2o.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for water and oxygen.
 */
class H2O_CO2
{
public:
  /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for molecular CO2 in liquid water.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        const Scalar E = 1672.9376;
        const Scalar F = 28.1751;
        const Scalar G = -112.4619;
        const Scalar H = 85.3807;

        return henryIAPWS(E, F, G, H, temperature);
    }


    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ of water in the CO2 phase.
     *
     * According to B. Xu et al. (2002) \cite xu2003 <BR>
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
	template <class Scalar, class CO2Table >
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
		using H2O = Components::H2O<Scalar>;
        using CO2 = Components::CO2<Scalar, CO2Table>;
		
        static const bool hasGasDiffCoeff = hasParam("BinaryCoefficients.GasDiffCoeff");
        if (!hasGasDiffCoeff) //in case one might set that user-specific as e.g. in dumux-lecture/mm/convectivemixing
        {
            //Diffusion coefficient of water in the CO2 phase
            constexpr Scalar PI = 3.141593;
            constexpr Scalar k = 1.3806504e-23; // Boltzmann constant
            constexpr Scalar c = 4; // slip parameter, can vary between 4 (slip condition) and 6 (stick condition)
            constexpr Scalar R_h = 1.72e-10; // hydrodynamic radius of the solute
            const Scalar mu = CO2::gasViscosity(temperature, pressure); // CO2 viscosity
            return k / (c * PI * R_h) * (temperature / mu);
        }
        else
        {
            static const Scalar D = getParam<Scalar>("BinaryCoefficients.GasDiffCoeff");
            return D;
        }
    }
    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular CO2 in liquid water. // used the one for O2
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     *
     * The empirical equations for estimating the diffusion coefficient in
     * infinite solution which are presented in Reid, 1987 all show a
     * linear dependency on temperature. We thus simply scale the
     * experimentally obtained diffusion coefficient of Ferrell and
     * Himmelblau by the temperature.
     *
     * See:
     *
     * R. Reid et al. (1987, pp. 599) \cite reid1987 <BR>
     *
     * R. Ferrell, D. Himmelblau (1967, pp. 111-115) \cite ferrell1967
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        const Scalar Texp = 273.15 + 25; // [K]
        const Scalar Dexp = 2.0e-9; // [m^2/s] //refer the value in brine_co2.hh
        return Dexp * temperature/Texp;
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
