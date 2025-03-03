/*****************************************************************************
 *   Copyright (C) 2012 by Johannes Hommel                                   *
 *                                                                           *
 *   Copyright (C) 2008-2010 by Melanie Darcis                               *
 *   Copyright (C) 2009-2010 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A fluid system with water and gas as phases and brine and CO2
 *        as components.
 */
#ifndef DUMUX_LEO_MIN_SYSTEM_HH
#define DUMUX_LEO_MIN_SYSTEM_HH

#include <dumux/material/idealgas.hh>
#include <dumux/material/constants.hh>


#include <dumux/material/fluidsystems/base.hh>

#include <dumux/material/fluidstates/adapter.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/n2.hh>
#include <dumux/material/components/tabulatedcomponent.hh>

#include "components/hydronion.hh"
#include "components/hydroxideion.hh"
#include <dumux/material/components/carbonateion.hh>
#include "components/bicarbonateion.hh"

#include <dumux/material/binarycoefficients/h2o_n2.hh>
#include "binarycoefficients/h2o_co2.hh" //added by du
#include "binarycoefficients/n2_co2.hh" //added by du
#include <dumux/material/binarycoefficients/h2o_o2.hh> //added by du
#include <dumux/material/binarycoefficients/n2_o2.hh> //added by du

#include <dumux/common/valgrind.hh>
#include <dumux/common/exceptions.hh>

// #include <assert.h>

// #ifdef DUMUX_PROPERTIES_HH
// #include <dumux/common/propertysystem.hh>
// #include <dumux/common/basicproperties.hh>
// #endif

namespace Dumux
{
namespace FluidSystems
{
/*!
 * \brief A compositional fluid with brine and carbon as
 *        components in both, the liquid and the gas (supercritical) phase,
 *        additional biomineralisation components and solid phases.
 *
 * This class provides acess to the Bio fluid system when no property system is used.
 * For Dumux users, using EnzymeMinFluid<TypeTag> and the documentation therein is
 * recommended.
 *
 *  The user can provide their own material table for CO2 properties.
 *  This fluidsystem is initialized as default with the tabulated version of
 *  water of the IAPWS-formulation, and the tabularized adapter to transfer
 *  this into brine.
 *  In the non-TypeTagged version, salinity information has to be provided with
 *  the init() methods.
 */
template<bool fastButSimplifiedRelations = true>
struct H2ON2CO2DefaultPolicy
{
    static constexpr bool useH2ODensityAsLiquidMixtureDensity() { return fastButSimplifiedRelations; }
    static constexpr bool useIdealGasDensity() { return fastButSimplifiedRelations; }
    static constexpr bool useN2ViscosityAsGasMixtureViscosity() { return fastButSimplifiedRelations; }
    static constexpr bool useN2HeatConductivityAsGasMixtureHeatConductivity() { return fastButSimplifiedRelations; }
    static constexpr bool useIdealGasHeatCapacities() { return fastButSimplifiedRelations; }
};

template <class Scalar,
          class CO2Tables,
          class H2OType = Components::TabulatedComponent<Components::H2O<Scalar>>,
          class Policy = H2ON2CO2DefaultPolicy<>	>
class LeoMinFluid
: public Base<Scalar, LeoMinFluid<Scalar, CO2Tables, H2OType, Policy> >

{
    using ThisType = LeoMinFluid<Scalar, CO2Tables, H2OType, Policy>;
    using Base = Dumux::FluidSystems::Base<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;
    using Constants = Dumux::Constants<Scalar>;
	
public:
    using H2O = H2OType;
	using CO2 = Components::CO2<Scalar,CO2Tables>;	
	using N2 = Components::N2<Scalar>;
    using CO3 = Components::CarbonateIon<Scalar>;
    using HCO3 = Components::BicarbonateIon<Scalar>;
    using H = Components::HydronIon<Scalar>;
    using OH = Components::HydroxideIon<Scalar>;
	using H2O_CO2 = BinaryCoeff::H2O_CO2;
	using N2_CO2 = BinaryCoeff::N2_CO2;

    // the type of parameter cache objects. this fluid system does not
    // cache anything, so it uses Dumux::NullParameterCache
    typedef Dumux::NullParameterCache ParameterCache;

public:

    /****************************************
     * Fluid phase related static parameters
     ****************************************/
    static const int numPhases = 2; // liquid and gas phases
    static const int wPhaseIdx = 0; // index of the liquid phase
    static const int nPhaseIdx = 1; // index of the gas phase
    static const int phase0Idx = wPhaseIdx;
    static const int phase1Idx = nPhaseIdx;

    /*!
     * \brief Return the human readable name of a fluid phase
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static std::string phaseName(int phaseIdx)
    {
        static std::string name[] = {
            "w",
            "n"
        };
        assert(0 <= phaseIdx && phaseIdx < numPhases);
        return name[phaseIdx];
    }

    /*!
     * \brief Return whether a phase is liquid
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isLiquid(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx == wPhaseIdx;
    }
    /*!
     * \brief Returns whether the fluids are miscible
     */
    static constexpr bool isMiscible()
    { return true; }
    /*!
     * \brief Return whether a phase is gaseous
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static constexpr bool isGas(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return phaseIdx == nPhaseIdx;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be an ideal mixture.
     *
     * We define an ideal mixture as a fluid phase where the fugacity
     * coefficients of all components times the pressure of the phase
     * are indepent on the fluid composition. This assumtion is true
     * if Henry's law and Rault's law apply. If you are unsure what
     * this function should return, it is safe to return false. The
     * only damage done will be (slightly) increased computation times
     * in some cases.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isIdealMixture(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /*!
     * \brief Returns true if and only if a fluid phase is assumed to
     *        be compressible.
     *
     * Compressible means that the partial derivative of the density
     * to the fluid pressure is always larger than zero.
     *
     * \param phaseIdx The index of the fluid phase to consider
     */
    static bool isCompressible(int phaseIdx)
    {
        assert(0 <= phaseIdx && phaseIdx < numPhases);

        return true;
    }

    /****************************************
     * Component related static parameters
     ****************************************/

    static const int numComponents = 4; // Water, N2,  CO2, H
    // static const int numMajorComponents = 2; //Water, N2
	static const int numSecComponents = 5; // h_total, co3, oh,  hco3, co2(aq)_total
    static const int H2OIdx =  0;
	static const int N2Idx =  1; 
    static const int wCompIdx =  H2OIdx;
    static const int nCompIdx =  N2Idx;
    static const int comp0Idx = wCompIdx;


    static const int CO2aqIdx  = 2;	
    static const int HIdx =  3;

    static const int OHIdx =  numComponents;
    static const int HCO3Idx = numComponents + 1;
    static const int CO3Idx = numComponents + 2;
    static const int CO2aqtotalIdx =  numComponents + 3;
    static const int HtotalIdx =  numComponents + 4;	

    /*!
     * \brief Return the human readable name of a component
     *
     * \param compIdx The index of the component to consider
     */
    static std::string componentName(int compIdx)
    {

        switch (compIdx) {
        case H2OIdx: return H2O::name();
        case N2Idx: return N2::name();
        case HtotalIdx: return "TotalH";
        case HIdx: return H::name();		
        case OHIdx: return OH::name();
        case HCO3Idx: return HCO3::name();
        case CO3Idx: return CO3::name();
		case CO2aqIdx: return "CO2aq";
		case CO2aqtotalIdx: return "TotalC";
        default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        };
    }

    /*!
     * \brief Return the molar mass of a component in [kg/mol].
     *
     * \param compIdx The index of the component to consider
     */
    static Scalar molarMass(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return H2O::molarMass();
// actually, the molar mass of brine is only needed for diffusion
// but since chloride and sodium are accounted for seperately
// only the molar mass of water is returned.
        case N2Idx: return N2::molarMass();
        case CO2aqIdx: return CO2::molarMass();
        case CO2aqtotalIdx: return CO2::molarMass();
        case HIdx: return H::molarMass();
        case HtotalIdx: return H::molarMass();		
        case OHIdx: return OH::molarMass();
        case HCO3Idx: return HCO3::molarMass();
        case CO3Idx: return CO3::molarMass();
        default: DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        };
    }

    /*!
     * \brief Return the charge value of a component.
     */
    static Scalar charge(int compIdx)
    {
        switch (compIdx) {
        case H2OIdx: return 0;
        case N2Idx: return 0;
        case HIdx: return H::charge();
        case HtotalIdx: return H::charge();		
        case OHIdx: return OH::charge();
        case CO2aqIdx: return 0;
        case CO2aqtotalIdx: return 0;
        case HCO3Idx: return HCO3::charge();
        case CO3Idx: return CO3::charge();
        default:DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
        }
    }


public:


    /****************************************
     * thermodynamic relations
     ****************************************/

    /*!
     * \brief Initialize the fluid system's static parameters generically
     *
     * If a tabulated H2O component is used, we do our best to create
     * tables that always work.
     */
    static void init()
    {
        init(/*tempMin=*/273.15,
             /*tempMax=*/623.15,
             /*numTemp=*/100,
             /*pMin=*/0.0,
             /*pMax=*/20e6,
             /*numP=*/200);
    }

    /*!
     * \brief Initialize the fluid system's static parameters using
     *        problem specific temperature and pressure ranges
     *
     * \param tempMin The minimum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param tempMax The maximum temperature used for tabulation of water \f$\mathrm{[K]}\f$
     * \param nTemp The number of ticks on the temperature axis of the  table of water
     * \param pressMin The minimum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param pressMax The maximum pressure used for tabulation of water \f$\mathrm{[Pa]}\f$
     * \param nPress The number of ticks on the pressure axis of the  table of water
     */
    static void init(Scalar tempMin, Scalar tempMax, unsigned nTemp,
                     Scalar pressMin, Scalar pressMax, unsigned nPress)
    {
        std::cout << "The H2O-N2-CO2 fluid system was configured with the following policy:\n";
        std::cout << " - use H2O density as liquid mixture density: " << std::boolalpha << Policy::useH2ODensityAsLiquidMixtureDensity() << "\n";
        std::cout << " - use ideal gas density: " << std::boolalpha << Policy::useIdealGasDensity() << "\n";
        std::cout << " - use N2 viscosity as gas mixture viscosity: " << std::boolalpha << Policy::useN2ViscosityAsGasMixtureViscosity() << "\n";
        std::cout << " - use N2 heat conductivity as gas mixture heat conductivity: " << std::boolalpha << Policy::useN2HeatConductivityAsGasMixtureHeatConductivity() << "\n";
        std::cout << " - use ideal gas heat capacities: " << std::boolalpha << Policy::useIdealGasHeatCapacities() << std::endl;

        if (H2O::isTabulated)
        {
            H2O::init(tempMin, tempMax, nTemp,
                               pressMin, pressMax, nPress);
        }
    }

    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * If Policy::useH2ODensityAsLiquidMixtureDensity() == false, we apply Eq. (7)
     * in Class et al. (2002a) \cite A3:class:2002b <BR>
     * for the liquid density.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx)
        {
            // assume pure water
            if (Policy::useH2ODensityAsLiquidMixtureDensity())
                return H2O::liquidDensity(T, 150000);//p);

            // See: Eq. (7) in Class et al. (2002a)
            // This assumes each gas molecule displaces exactly one
            // molecule in the liquid.
            else
                return H2O::liquidMolarDensity(T, p)
                       * (fluidState.moleFraction(wPhaseIdx, H2OIdx)*H2O::molarMass()
                          + fluidState.moleFraction(wPhaseIdx, N2Idx)*N2::molarMass()
                          + fluidState.moleFraction(wPhaseIdx, CO2aqIdx)*CO2::molarMass());
        }

        // gas phase
        else if (phaseIdx == nPhaseIdx)
        {

            // for the gas phase assume an ideal gas
            using std::max;
            if (Policy::useIdealGasDensity())
                return IdealGas::molarDensity(T, p) * fluidState.averageMolarMass(nPhaseIdx);

            // assume ideal mixture: steam, nitrogen and oxygen don't "see" each other
            else
                return H2O::gasDensity(T, fluidState.partialPressure(nPhaseIdx, H2OIdx))
                       + N2::gasDensity(T, fluidState.partialPressure(nPhaseIdx, N2Idx))
                       + CO2::gasDensity(T, fluidState.partialPressure(nPhaseIdx, CO2aqIdx));
        }

        DUNE_THROW(Dune::InvalidStateException, "Unknown phase index " << phaseIdx);
    }

    using Base::molarDensity;
    //! \copydoc Base<Scalar,ThisType>::molarDensity(const FluidState&,int)
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            // assume pure water or that each gas molecule displaces exactly one
            // molecule in the liquid.
            return H2O::liquidMolarDensity(T, 150000);//p);
        }
        else
        {
            if (Policy::useIdealGasDensity())
            {   //assume ideal gas
                return IdealGas::molarDensity(T,p);
            }

            return H2O::gasMolarDensity(T, fluidState.partialPressure(nPhaseIdx, H2OIdx))
                   + N2::gasMolarDensity(T, fluidState.partialPressure(nPhaseIdx, N2Idx))
                   + CO2::gasMolarDensity(T, fluidState.partialPressure(nPhaseIdx, CO2aqIdx));
        }
    }

    using Base::viscosity;
    /*!
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * Compositional effects in the gas phase are accounted by the Wilke method.
     * See Reid et al. (1987)  \cite reid1987 <BR>
     * 4th edition, McGraw-Hill, 1987, 407-410
     * 5th edition, McGraw-Hill, 20001, p. 9.21/22
     * \note Compositional effects for a liquid mixture have to be implemented.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
            // assume pure water for the liquid phase
			// printf("The value of H2O::liquidViscosity(T, p) is: %.e\n", H2O::liquidViscosity(T, p)); 
            return H2O::liquidViscosity(T, 150000);//p);
			
        }

        // gas phase
        if (Policy::useN2ViscosityAsGasMixtureViscosity())
        {
            // assume pure nitrogen for the gas phase
            return N2::gasViscosity(T, p); //2e-5
        }
        else
        {
            // Wilke method (Reid et al.):
            Scalar muResult = 0;
            const Scalar mu[numComponents] = {
                // h2oGasViscosityInMixture(T, p), //changed, but need to figure out where went wrong
				H2O::gasViscosity(T, p),
                N2::gasViscosity(T, p),
                CO2::gasViscosity(T, p)
            };

            Scalar sumx = 0.0;
            using std::max;
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                sumx += fluidState.moleFraction(phaseIdx, compIdx);
            sumx = max(1e-10, sumx);

            for (int i = 0; i < numComponents; ++i) {
                Scalar divisor = 0;
                using std::pow;
                using std::sqrt;
                for (int j = 0; j < numComponents; ++j) {
                    Scalar phiIJ = 1 + sqrt(mu[i]/mu[j]) * pow(molarMass(j)/molarMass(i), 1/4.0);
                    phiIJ *= phiIJ;
                    phiIJ /= sqrt(8*(1 + molarMass(i)/molarMass(j)));
                    divisor += fluidState.moleFraction(phaseIdx, j)/sumx * phiIJ;
                }
                muResult += fluidState.moleFraction(phaseIdx, i)/sumx * mu[i] / divisor;
            }
            return muResult;
        }
    }

    using Base::fugacityCoefficient;
    /*!
     * \brief Returns the fugacity coefficient \f$\mathrm{[-]}\f$ of a component in a
     *        phase.
     *
     * The fugacity coefficient \f$\phi^\kappa_\alpha\f$ of
     * component \f$\kappa\f$ in phase \f$\alpha\f$ is connected to
     * the fugacity \f$f^\kappa_\alpha\f$ and the component's mole
     * fraction \f$x^\kappa_\alpha\f$ by means of the relation
     *
     * \f[
     f^\kappa_\alpha = \phi^\kappa_\alpha\;x^\kappa_\alpha\;p_\alpha
     \f]
     * where \f$p_\alpha\f$ is the pressure of the fluid phase.
     *
     * For liquids with very low miscibility this boils down to the
     * Henry constant for the solutes and the saturated vapor pressure
     * both divided by phase pressure.
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar fugacityCoefficient(const FluidState &fluidState,
                                      int phaseIdx,
                                      int compIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        assert(0 <= compIdx  && compIdx < numComponents+numSecComponents); // added by du

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        
		if (phaseIdx == wPhaseIdx)
        {
            switch(compIdx){
                case H2OIdx: return H2O::vaporPressure(T)/p;

                case N2Idx: return BinaryCoeff::H2O_N2::henry(T)/p;
				case CO2aqIdx: return BinaryCoeff::H2O_CO2::henry(T)/p;
				case HIdx: return 0.0;
				
            };
        }

        // for the gas phase, assume an ideal gas when it comes to
        // fugacity (-> fugacity == partial pressure)
        return 1.0;

    }
   using Base::diffusionCoefficient;
    //! \copydoc Base<Scalar,ThisType>::diffusionCoefficient(const FluidState&,int,int)
    template <class FluidState>
    static Scalar diffusionCoefficient(const FluidState &fluidState,
                                       int phaseIdx,
                                       int compIdx)
    {
        DUNE_THROW(Dune::NotImplemented, "Diffusion coefficients");
    }

    using Base::binaryDiffusionCoefficient;
    //! \copydoc Base<Scalar,ThisType>::binaryDiffusionCoefficient(const FluidState&,int,int,int)
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)

    {
        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

#ifndef NDEBUG
        if (compIIdx == compJIdx ||
            phaseIdx > numPhases - 1 ||
            compJIdx > numComponents - 1)
        {
            DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
#endif

        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        // reached here printf("The value of (T, p)");
		
        // gas phase
        if (phaseIdx == nPhaseIdx) {

			if (compIIdx == N2Idx && compJIdx == H2OIdx)
            {
				// printf("The value of N2_H2O::liquidDiffCoeff(T, p) is: %.e\n", BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p));
				return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);
		    }
			if (compIIdx == N2Idx && compJIdx == CO2aqIdx)
            {
				// printf("The value of N2_CO2::liquidDiffCoeff(T, p) is: %.e\n", N2_CO2::gasDiffCoeff<Scalar, CO2Table>(T, p));
				return BinaryCoeff::N2_CO2::gasDiffCoeff<Scalar, CO2Tables>(T, p); // 2e-5
				// return 4e-7;//0.0;
		    }

			else if (compJIdx <numComponents)
            {
				return 1e-12;
            }
			else 
				DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
        // liquid phase
        else if (phaseIdx == wPhaseIdx) {
            if (compIIdx == H2OIdx && compJIdx == N2Idx)
			{
				// printf("The value of H2O_N2::liquidDiffCoeff(T, p) is: %.e\n", BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p));
				return BinaryCoeff::H2O_N2::liquidDiffCoeff(T, p);
				// return 0.0;
			}
            if (compIIdx == H2OIdx && compJIdx == CO2aqIdx)
			{
				// printf("The value of H2O_CO2::liquidDiffCoeff(T, p) is: %.e\n", H2O_CO2::liquidDiffCoeff(T, p));
				return H2O_CO2::liquidDiffCoeff(T, p); //2e-09 
				// return 4e-11;
			    // return 0.0;
			}
			else if (compJIdx <numComponents)
			{
				return 1.587e-9;        //[m²/s]        //J. Phys. D: Appl. Phys. 40 (2007) 2769-2776 //old Value from Anozie 1e-9
                // return 0.0;
				// printf("The value of liquidDiffCoeff(T, p)");
			}
            else 
				DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx
                       << " in phase " << phaseIdx << " is undefined!\n");
        }
        DUNE_THROW(Dune::InvalidStateException,
                  "Binary diffusion coefficient of components "
                  << compIIdx << " and " << compJIdx
                  << " in phase " << phaseIdx << " is undefined!\n");
    }

    using Base::enthalpy;
    /*!
     * \brief Given a phase's composition, temperature, pressure and
     *        density, calculate its specific enthalpy \f$\mathrm{[J/kg]}\f$.
     *
     *  \note This fluid system neglects the contribution of
     *        gas-molecules in the liquid phase. This contribution is
     *        probably not big. Somebody would have to find out the
     *        enthalpy of solution for this system. ...
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar enthalpy(const FluidState &fluidState,
                           int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);

        // liquid phase
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidEnthalpy(T, p);
        }
        // gas phase
        else if (phaseIdx == nPhaseIdx)
        {
            // assume ideal mixture: which means
            // that the total specific enthalpy is the sum of the
            // "partial specific enthalpies" of the components.
            Scalar hH2O =
                fluidState.massFraction(nPhaseIdx, H2OIdx)
                * H2O::gasEnthalpy(T, p);
            Scalar hN2 =
                fluidState.massFraction(nPhaseIdx, N2Idx)
                * N2::gasEnthalpy(T,p);
            Scalar hCO2 =
                fluidState.massFraction(nPhaseIdx, CO2aqIdx)
                * CO2::gasEnthalpy(T,p);
            return hH2O + hN2 + hCO2;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the specific enthalpy \f$\mathrm{[J/kg]}\f$ of a component in a specific phase
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param componentIdx The index of the component to consider
     */
    template <class FluidState>
    static Scalar componentEnthalpy(const FluidState &fluidState,
                                    int phaseIdx,
                                    int componentIdx)
    {
        const Scalar T = fluidState.temperature(phaseIdx);
        const Scalar p = fluidState.pressure(phaseIdx);

        if (phaseIdx == phase0Idx)
        {
            if (componentIdx == H2OIdx)
                return H2O::liquidEnthalpy(T, p);
            else if (componentIdx == N2Idx)
                DUNE_THROW(Dune::NotImplemented, "Component enthalpy of nitrogen in liquid phase");
            else if (componentIdx == CO2aqIdx)
                DUNE_THROW(Dune::NotImplemented, "Component enthalpy of carbon dioxide in liquid phase");
            else
                DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        else if (phaseIdx == phase1Idx)
        {
            if (componentIdx == H2OIdx)
                return H2O::gasEnthalpy(T, p);
            else if (componentIdx == N2Idx)
                return N2::gasEnthalpy(T, p);
            else if (componentIdx == CO2aqIdx)
                return CO2::gasEnthalpy(T, p);
            DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << componentIdx);
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    using Base::thermalConductivity;
    /*!
     * \brief Thermal conductivity of a fluid phase \f$\mathrm{[W/(m K)]}\f$.
     *
     * Use the conductivity of air and water as a first approximation.
     *
     * http://en.wikipedia.org/wiki/List_of_thermal_conductivities
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar thermalConductivity(const FluidState &fluidState,
                                      const int phaseIdx)
    {
        assert(0 <= phaseIdx  && phaseIdx < numPhases);
        Scalar temperature  = fluidState.temperature(phaseIdx) ;
        Scalar pressure = fluidState.pressure(phaseIdx);

        if (phaseIdx == wPhaseIdx)
        {
            return H2O::liquidThermalConductivity(temperature, pressure);
        }
        else
        {
            Scalar lambdaPureN2 = N2::gasThermalConductivity(temperature, pressure);
            Scalar lambdaPureCO2 = CO2::gasThermalConductivity(temperature, pressure);
            if (!Policy::useN2HeatConductivityAsGasMixtureHeatConductivity())
            {
                Scalar xN2 = fluidState.moleFraction(phaseIdx, N2Idx);
                Scalar xCO2 = fluidState.moleFraction(phaseIdx, CO2aqIdx);
                Scalar xH2O = fluidState.moleFraction(phaseIdx, H2OIdx);
                Scalar lambdaN2 = xN2 * lambdaPureN2;
                Scalar lambdaCO2 = xCO2 * lambdaPureCO2;
                Scalar partialPressure  = pressure * xH2O;
                Scalar lambdaH2O = xH2O * H2O::gasThermalConductivity(temperature, partialPressure);
                return lambdaN2 + lambdaH2O + lambdaCO2;
            }
            else
                return lambdaPureN2;
        }
    }

    using Base::heatCapacity;
    //! \copydoc Base<Scalar,ThisType>::heatCapacity(const FluidState&,int)
    template <class FluidState>
    static Scalar heatCapacity(const FluidState &fluidState,
                               int phaseIdx)
    {
        if (phaseIdx == wPhaseIdx) {
            return H2O::liquidHeatCapacity(fluidState.temperature(phaseIdx),
                                           fluidState.pressure(phaseIdx));
        }

        Scalar c_pN2;
        Scalar c_pCO2;
        Scalar c_pH2O;
        // let the water and nitrogen components do things their own way
        if (!Policy::useIdealGasHeatCapacities()) {
            c_pN2 = N2::gasHeatCapacity(fluidState.temperature(phaseIdx),
                                        fluidState.pressure(phaseIdx)
                                        * fluidState.moleFraction(phaseIdx, N2Idx));
            c_pH2O = H2O::gasHeatCapacity(fluidState.temperature(phaseIdx),
                                          fluidState.pressure(phaseIdx)
                                          * fluidState.moleFraction(phaseIdx, H2OIdx));
            c_pCO2 = CO2::gasHeatCapacity(fluidState.temperature(phaseIdx),
                                        fluidState.pressure(phaseIdx)
                                        * fluidState.moleFraction(phaseIdx, CO2aqIdx));
        }
        else {
            // assume an ideal gas for both components. See:
            //
            //http://en.wikipedia.org/wiki/Heat_capacity
            Scalar c_vN2molar = Constants::R*2.50;
            Scalar c_pN2molar = Constants::R + c_vN2molar;

            Scalar c_vCO2molar = Constants::R*2.43;
            Scalar c_pCO2molar = Constants::R + c_vCO2molar;

            Scalar c_vH2Omolar = Constants::R*3.25; // <- correct??
            Scalar c_pH2Omolar = Constants::R + c_vH2Omolar;
			
            c_pN2 = c_pN2molar/molarMass(N2Idx);
            c_pCO2 = c_pCO2molar/molarMass(CO2aqIdx);
            c_pH2O = c_pH2Omolar/molarMass(H2OIdx);
        }

        // mangle all components together
        return
            c_pH2O*fluidState.massFraction(nPhaseIdx, H2OIdx)
            + c_pN2*fluidState.massFraction(nPhaseIdx, N2Idx)
            + c_pCO2*fluidState.massFraction(nPhaseIdx, CO2aqIdx);
    }

};

} // end namespace FluidSystems
} // end namespace Dumux

#endif
