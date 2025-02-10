// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for nitrogen and CO2.
 */
#ifndef DUMUX_BINARY_COEFF_N2_CO2_HH
#define DUMUX_BINARY_COEFF_N2_CO2_HH

#include <dumux/material/binarycoefficients/henryiapws.hh>
#include <dumux/material/binarycoefficients/fullermethod.hh>

#include <dumux/material/components/co2.hh>
#include <dumux/material/components/n2.hh>

namespace Dumux {
namespace BinaryCoeff {

/*!
 * \ingroup Binarycoefficients
 * \brief Binary coefficients for nitrogen and CO2.
 */
template <class Scalar>
class N2_CO2 {
public:
    /*!
     * \brief Henry coefficient \f$\mathrm{[Pa]}\f$  for molecular CO2 in liquid nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     */
    template <class Scalar>
    static Scalar henry(Scalar temperature)
    {
        DUNE_THROW(Dune::NotImplemented, "henry coefficient for CO2 in liquid nitrogen");
    }

    /*!
     * \brief Binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular CO2 in liquid nitrogen.
     *
     * Uses fullerMethod to determine the diffusion of water in nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar gasDiffCoeff(Scalar temperature, Scalar pressure)
    {
        using N2 = Dumux::Components::N2<Scalar>;
        using CO2 = Dumux::Components::CO2<Scalar>;

        // atomic diffusion volumes
        const Scalar SigmaNu[2] = { 18.1 /* N2 */,  26.9 /* CO2 */ };
        // molar masses [g/mol]
        const Scalar M[2] = { N2::molarMass()*1e3, CO2::molarMass()*1e3 };
        return fullerMethod(M, SigmaNu, temperature, pressure);
    }

    /*!
     * \brief Diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for molecular carbon dioxide in liquid nitrogen.
     * \param temperature the temperature \f$\mathrm{[K]}\f$
     * \param pressure the phase pressure \f$\mathrm{[Pa]}\f$
     */
    template <class Scalar>
    static Scalar liquidDiffCoeff(Scalar temperature, Scalar pressure)
    {
        DUNE_THROW(Dune::NotImplemented, "diffusion coefficient for liquid carbon dioxide and nitrogen");
    }
};

} // end namespace BinaryCoeff
} // end namespace Dumux

#endif
