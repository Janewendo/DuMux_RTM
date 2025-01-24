// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup TwoPICPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase induced calcium carbonate precipitation model.
 */
#ifndef DUMUX_2PICP_VOLUME_VARIABLES_HH
#define DUMUX_2PICP_VOLUME_VARIABLES_HH

// #include <iostream>
// #include <vector>

// #include <dumux/common/math.hh>
// #include <dumux/common/properties.hh>
// #include <dumux/discretization/methods.hh>
#include <dumux/porousmediumflow/2pnc/volumevariables.hh>
#include <dumux/material/fluidstates/compositionalsecondarycomponent.hh>
// #include <dumux/material/constraintsolvers/miscible2pnccomposition.hh>
// #include <dumux/material/constraintsolvers/computefromreferencephase.hh>

#include <dumux/material/solidstates/updatesolidvolumefractions.hh>

// #include "indices.hh" // for formulation

#include <dumux/material/binarycoefficients/brine_co2.hh>

namespace Dumux
{

/*!
 * \ingroup TwoPICPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase induced calcium carbonate precipitation model.
 */
template <class Traits>
class TwoPICPVolumeVariables
: public TwoPNCVolumeVariables<Traits>
, public EnergyVolumeVariables<Traits, TwoPICPVolumeVariables<Traits> >
{
    using ParentType = TwoPNCVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPICPVolumeVariables<Traits> >;
//     using Implementation = GetPropType<TypeTag, Properties::VolumeVariables>;
//     using GridView = GetPropType<TypeTag, Properties::GridView>;
//     using Problem = GetPropType<TypeTag, Properties::Problem>;
//     using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
//     using FVElementGeometry = typename GetPropType<TypeTag, Properties::FVGridGeometry>::LocalView;
//     using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using PermeabilityType = typename Traits::PermeabilityType;
    using FS =  typename Traits::FluidSystem;
    using SSY = typename Traits::SolidSystem; // added by du
    using ModelTraits = typename Traits::ModelTraits;

    using CO2Tables = typename Traits::CO2Tables;
    using Chemistry = typename Traits::Chemistry;
    using Brine_CO2 = Dumux::BinaryCoeff::Brine_CO2<Scalar, CO2Tables, true>;

    // static constexpr int numFluidComps = ParentType::numFluidComponents();
    static constexpr int numFluidComps = ModelTraits::numFluidComponents();
    enum
    {
//         numPhases = ModelTraits::numPhases(),
        numMajorComponents = ModelTraits::numPhases(),
        numComponents = ModelTraits::numFluidComponents(),
        numSecComponents = ModelTraits::numSecComponents(),

        // phase indices
        wPhaseIdx = FS::wPhaseIdx,
        nPhaseIdx = FS::nPhaseIdx,
		// gPhaseIdx = SSY::GlassIdx, // added by du
		gPhaseIdx = SSY::GlassIdx, // added by du
		// cPhaseIdx = SSY::CalciteIdx, // added by du
		fPhaseIdx = SSY::FerrohydriteIdx, // added by du
		pPhaseIdx = SSY::ProtoImogoliteIdx, // added by du
		bPhaseIdx = SSY::BirnessiteIdx, // added by du
		hPhaseIdx = SSY::HydroxyapatiteIdx, // added by du
		sPhaseIdx = SSY::SepioliteIdx, // added by du
        phase0Idx = wPhaseIdx,
        phase1Idx = nPhaseIdx,

        // component indices
        wCompIdx = FS::wCompIdx,
        nCompIdx = FS::nCompIdx,
        comp0Idx = wCompIdx,
        comp1Idx = nCompIdx,
        NaIdx = FS::NaIdx,
        CaIdx = FS::CaIdx,
        // ClIdx = FS::ClIdx,
        // CtotIdx = FS::CtotIdx,
        HIdx = FS::HIdx,
        // CO2gIdx = FS::CO2gIdx,
		CO2aqIdx = FS::CO2aqIdx,
		CO2aqonlyIdx = FS::CO2aqonlyIdx,
        // phase presence enums
        secondPhaseOnly = ModelTraits::Indices::secondPhaseOnly,
        firstPhaseOnly = ModelTraits::Indices::firstPhaseOnly,
        bothPhases = ModelTraits::Indices::bothPhases,
        wPhaseOnly = firstPhaseOnly,
        nPhaseOnly = secondPhaseOnly,

        // primary variable indices
        pressureIdx = ModelTraits::Indices::pressureIdx,
        switchIdx = ModelTraits::Indices::switchIdx
    };

    static constexpr auto formulation = ModelTraits::priVarFormulation();
    static constexpr bool setFirstPhaseMoleFractions = ModelTraits::setMoleFractionsForFirstPhase();

    using MiscibleMultiPhaseComposition = Dumux::MiscibleMultiPhaseComposition<Scalar, FS>;
    using ComputeFromReferencePhase = Dumux::ComputeFromReferencePhase<Scalar, FS>;

public:
    //! return number of secondary components considered by the model
    static constexpr int numSecFluidComponents() { return Traits::ModelTraits::numSecComponents(); }
    //! export fluid state type
    using FluidState = typename Traits::FluidState;
    //! export fluid system type
    using FluidSystem = typename Traits::FluidSystem;
    //! export type of solid state
    using SolidState = typename Traits::SolidState;
    //! export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    //! return whether moles or masses are balanced
    static constexpr bool useMoles() { return Traits::ModelTraits::useMoles(); }
    //! return the two-phase formulation used here
    static constexpr TwoPFormulation priVarFormulation() { return formulation; }
    // check for permissive specifications
    static_assert(useMoles(), "use moles has to be set true in the 2pnc model");
    static_assert(ModelTraits::numPhases() == 2, "NumPhases set in the model is not two!");
    static_assert((formulation == TwoPFormulation::p0s1 || formulation == TwoPFormulation::p1s0), "Chosen TwoPFormulation not supported!");

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        // added by du
		// const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);

        // precompute the minimum capillary pressure (entry pressure)
        // needed to make sure we don't compute unphysical capillary pressures and thus saturations
        // minPc_ = fluidMatrixInteraction.endPointPc();
        // added till here by du
		
        //we need porosity and permeability in completeFluidState for Leverett scaling of capillary pressure!
        updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_, numFluidComps);
		// printf("The numFluidComps is: %.10f\n", numFluidComps);
		// printf("The numFluidComps is: %.10f\n", numComponents);
		// printf("The solidvolumefraction glass is: %.20e\n", solidState().volumeFraction(gPhaseIdx));
		// printf("The solidvolumefraction proto is: %.20e\n", solidState().volumeFraction(pPhaseIdx));
		// printf("The solidvolumefraction ferr is: %.20e\n", solidState().volumeFraction(fPhaseIdx));
		// printf("The solidvolumefraction birn is: %.20e\n", solidState().volumeFraction(bPhaseIdx));
        
		porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
		// printf("Porosity is: %.14e\n", porosity_);	
		// if (porosity_ != 0.39)
		// {
		// 	printf("The porosity1 is: %.10e\n", porosity_);
		// }
        permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
        
        Scalar ii = 0;
        completeFluidState(elemSol, problem, element, scv, fluidState_, solidState_);

        /////////////
        // calculate the remaining quantities
        /////////////
        // Second instance of a parameter cache.
        // Could be avoided if diffusion coefficients also
        // became part of the fluid state.
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState_);

        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            // relative permeabilities
            Scalar kr;
            if (phaseIdx == wPhaseIdx)
                kr = MaterialLaw::krw(materialParams, saturation(wPhaseIdx));
            else // ATTENTION: krn requires the wetting-phase saturation as parameter!
                kr = MaterialLaw::krn(materialParams, saturation(wPhaseIdx));

            mobility_[phaseIdx] = kr / fluidState_.viscosity(phaseIdx);
        //     int compIIdx = phaseIdx;
        //     for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        //     {
        //         // binary diffusion coefficients
        //         if(compIIdx!= compJIdx)
        //         {
        //             setDiffusionCoefficient_(phaseIdx, compJIdx,
        //                                      FluidSystem::binaryDiffusionCoefficient(fluidState_,
        //                                                                              paramCache,
        //                                                                              phaseIdx,
        //                                                                              compIIdx,
        //                                                                              compJIdx));
        //         }
        //     }
        }
        // binary diffusion coefficients, added by du
        for (unsigned int compJIdx = 0; compJIdx < numComponents; ++compJIdx)
        {
            if(compJIdx != comp0Idx)
                setDiffusionCoefficient_( phase0Idx, compJIdx,
                                          FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                  paramCache,
                                                                                  phase0Idx,
                                                                                  comp0Idx,
                                                                                  compJIdx) );
            // printf("The value0 is: %.e\n", FluidSystem::binaryDiffusionCoefficient(fluidState_,
            //                                                                     paramCache,
            //                                                                     phase0Idx,
            //                                                                     comp0Idx,
            //                                                                     compJIdx)); 																				  
            if(compJIdx != comp1Idx)
                setDiffusionCoefficient_( phase1Idx, compJIdx,
                                          FluidSystem::binaryDiffusionCoefficient(fluidState_,
                                                                                  paramCache,
                                                                                  phase1Idx,
                                                                                  comp1Idx,
                                                                                  compJIdx) );
			// printf("The value1 is: %.e\n", FluidSystem::binaryDiffusionCoefficient(fluidState_,
            //                                                                    paramCache,
            //                                                                    phase1Idx,
            //                                                                    comp1Idx,
            //                                                                    compJIdx)); 																				  
    																	  
   
        }   
        // calculate the remaining quantities
        EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv, solidState_);
    }
    // Scalar init ii = 0;  
    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void completeFluidState(const ElemSol& elemSol,
                            const Problem& problem,
                            const Element& element,
                            const Scv& scv,
                            FluidState& fluidState,
                            SolidState& solidState)
    {
        EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState, solidState);

        const auto& priVars = elemSol[scv.localDofIndex()];
        const auto phasePresence = priVars.state();
        // added by du
		// const auto fluidMatrixInteraction = problem.spatialParams().fluidMatrixInteraction(element, scv, elemSol);
        // ended here 
		
        using MaterialLaw = typename Problem::SpatialParams::MaterialLaw;
        const auto& materialParams = problem.spatialParams().materialLawParams(element, scv, elemSol);
        const int wPhaseIdx = problem.spatialParams().template wettingPhase<FluidSystem>(element, scv, elemSol);
        fluidState.setWettingPhase(wPhaseIdx);

        // added by du
		// compute the capillary pressure to compute the saturation
        // make sure that we the capillary pressure is not smaller than the minimum pc
        // this would possibly return unphysical values from regularized material laws
        // using std::max;
        // const Scalar pc = max(minPc_, problem.nonwettingReferencePressure() - fluidState.pressure(wPhaseIdx));
        // const Scalar sw = fluidMatrixInteraction.sw(pc);
        // fluidState.setSaturation(phase0Idx, sw);
        // fluidState.setSaturation(phase1Idx, 1.0-sw);
		
        // set the saturations changed by du
        if (phasePresence == secondPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 0.0);
            fluidState.setSaturation(phase1Idx, 1.0);
        }
        else if (phasePresence == firstPhaseOnly)
        {
            fluidState.setSaturation(phase0Idx, 1.0);
            fluidState.setSaturation(phase1Idx, 0.0);
        }
        else if (phasePresence == bothPhases)
        {
            if (formulation == TwoPFormulation::p0s1)
            {
                fluidState.setSaturation(phase1Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase0Idx, 1 - priVars[switchIdx]);
            }
            else
            {
                fluidState.setSaturation(phase0Idx, priVars[switchIdx]);
                fluidState.setSaturation(phase1Idx, 1 - priVars[switchIdx]);
            	// printf("p1s0priVarsswitchIdx.\n");
	            // printf("The value of priVarsswitchIdx is: %.2lf\n", priVars[switchIdx]);

			}
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "phasePresence: " << phasePresence << " is invalid.");

        // set pressures of the fluid phases
        refPC_ = MaterialLaw::pc(materialParams, fluidState.saturation(wPhaseIdx));
		// printf("refPC.\n");
	    // printf("The value of refPC is: %.2lf\n", refPC_);
        refPorosity_     = getParam<Scalar>("SpatialParams.ReferencePorosity");
        refPermeability_ = getParam<Scalar>("SpatialParams.ReferencePermeability");
        pc_ = refPC_ * sqrt((refPermeability_*porosity_)/
                               (permeability_*refPorosity_));
        if (formulation == TwoPFormulation::p0s1)
        {
            fluidState.setPressure(phase0Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase1Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] + pc_
                                                                       : priVars[pressureIdx] - pc_);
        }
        else
        {
		    // printf("pc_.\n");
	        // printf("The value of PC is: %.2lf\n", pc_);
            fluidState.setPressure(phase1Idx, priVars[pressureIdx]);
            fluidState.setPressure(phase0Idx, (wPhaseIdx == phase0Idx) ? priVars[pressureIdx] - pc_
                                                                       : priVars[pressureIdx] + pc_);
        }


        // calculate the phase compositions
        typename FluidSystem::ParameterCache paramCache;

        // now comes the tricky part: calculate phase composition
        Dune::FieldVector<Scalar, numComponents + numSecComponents> moleFrac(0.0);
        // //Du moved here
		// for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
        // {
        //     fluidState.setMoleFraction(wPhaseIdx, compIdx, priVars[compIdx]);
        // }
		// // printf("The value of moleFrac[CO2aqIdx]1 is: %.8e\n", moleFrac[CO2aqIdx]); //  0.00000000e+00
        // for (int compIdx= 0; compIdx<numComponents; ++compIdx)
        // {
        //     moleFrac[compIdx] =fluidState.moleFraction(wPhaseIdx, compIdx);
		// }
		// // moleFrac[CO2aqonlyIdx] =2.38640000e-07;
        // // printf("moleFrac[CO2aqIdx]): %.8e\n", moleFrac[CO2aqIdx]);		
		// // printf("The value of moleFrac[CO2aqIdx]2 is: %.8e\n", moleFrac[CO2aqIdx]);	 // 2.65687305e-06	
		// Scalar rhoMolar = 0.0; // added by du Declare rhoMolar outside the loop and initialize it
		// rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phase0Idx);
		// // printf("The value of moleFracforchem0 is: %.5e\n", moleFrac[16]);				
	    // // printf("fluidState.moleFraction(0, CO2aqonlyIdx)0: %.8e\n", fluidState.moleFraction(0, CO2aqonlyIdx));
		// Chemistry chemistry;
		// // printf("reacehed here 1");
        // // auto variable = chemistry.calculateEquilibriumChemistry(fluidState, phasePresence, moleFrac, rhoMolar);
        // chemistry.calculateEquilibriumChemistry(fluidState, phasePresence, moleFrac, rhoMolar);
		// 
		// // printf("The value of moleFracforchem-1 is: %.5e\n", moleFrac[16]);			
        // for (int compIdx=numComponents; compIdx<numComponents + numSecComponents; ++compIdx)
        // {
        //     fluidState.setMoleFractionSecComp(phase0Idx, compIdx, moleFrac[compIdx]);
        //     fluidState.setMoleFractionSecComp(phase1Idx, compIdx, 0);
        // }

	    // printf("moleFrac[CO2aqonlyIdx]): %.8e\n", moleFrac[CO2aqonlyIdx]);
		// printf("variable[CO2aqonlyIdx]): %.8e\n", variable);//(CO2aqonlyIdx));
	    // printf("fluidState.moleFraction(0, CO2aqonlyIdx): %.8e\n", fluidState.moleFraction(0, CO2aqonlyIdx));
		// printf("The value of moleFrac[CO2aqonlyIdx1] is: %.8e\n", moleFrac[CO2aqonlyIdx]);// 2.38639366e-07
		// double ratio = 1.0;
		// printf("The value of moleFrac[CO2aqonlyIdx1] is: %.8e\n", fluidState.moleFraction(phase0Idx, CO2aqonlyIdx));
		// printf("The value of moleFrac[CO2aqIdx] is: %.8e\n", priVars[CO2aqIdx]);
		// ratio = fluidState.moleFraction(phase0Idx, CO2aqonlyIdx) / priVars[CO2aqIdx];
		// 
		// if (ratio < 0.0001)
		// {
		// 	ratio = 2.38640000e-07 / 2.65687305e-06;
		// 	printf("ratio1: %.5e\n",ratio);
		// }
		// Scalar co2aqn = fluidState.moleFraction(phase0Idx, CO2aqonlyIdx) * fluidState.fugacityCoefficient(phase0Idx, CO2aqIdx);
        // printf("The value of fluidState.moleFraction(phase0Idx, CO2aqonlyIdx) is: %.8e\n", fluidState.moleFraction(phase0Idx, CO2aqonlyIdx));
        // printf("The value of fluidState.fugacityCoefficient(phase0Idx, CO2aqIdx is: %.8e\n", fluidState.fugacityCoefficient(phase0Idx, CO2aqIdx));
        // printf("The value of co2aqn is: %.8e\n", co2aqn);
        // ii += 1;
		// printf("The value of ii is: %.8e\n", ii);
		if (phasePresence == bothPhases)
        {
            // both phases are present, phase composition results from
            // the first <-> second phase equilibrium. This is the job
            // of the "MiscibleMultiPhaseComposition" constraint solver

            // set the known mole fractions in the fluidState so that they
            // can be used by the MiscibleMultiPhaseComposition constraint solver

            const int knownPhaseIdx = setFirstPhaseMoleFractions ? phase0Idx : phase1Idx;
			// printf("The value of moleFracforchem1 is: %.5e\n", moleFrac[16]);
            for (int compIdx=numMajorComponents; compIdx<numComponents; ++compIdx)
            {
                fluidState.setMoleFraction(knownPhaseIdx, compIdx, priVars[compIdx]);
			
			// moleFrac[compIdx]=priVars[compIdx]; // added by du
		     // printf("The value of moleFrac1 is: %.e\n", fluidState.moleFraction(knownPhaseIdx,compIdx));
			// printf("The value of moleFracforchem is: %.e\n", moleFrac[compIdx]);
//              std::cout <<"both phases, priVars["<< FluidSystem::componentName(compIdx)<<"]"<<priVars[compIdx]<< std::endl;
            }
			
			// fluidState.setMoleFraction(phase0Idx, CO2aqIdx, moleFrac[CO2aqonlyIdx]); // added by Du, co2aqonly > co2aq
		    // printf("The value of moleFrac[CO2aqonlyIdx2] is: %.8e\n", moleFrac[CO2aqonlyIdx]); // 2.38639366e-07
			// printf("The value of moleFrac[CO2aqIdx1] is: %.8e\n", fluidState.moleFraction(wPhaseIdx, CO2aqIdx)); // 
			// printf("The value of  priVars[CO2aqIdx1] is: %.8e\n",  priVars[CO2aqIdx]); 
            //TODO: Only CO2 and water are present in the non-wetting phase, which is why we do not use the constraintsolver
			// printf("The value of moleFracforchem1 is: %.5e\n", moleFrac[16]);
            // ii += 1;
		    // printf("The value of ii is: %.8e\n", ii);
			// added by du
			// double fluidState.moleFraction(0, CO2aqonlyIdx) = 2.38640000e-07;
           for (int compIdx=numComponents; compIdx<numComponents+numSecComponents; ++compIdx)
            {
                fluidState.setMoleFractionSecComp(wPhaseIdx, compIdx, 2.38640000e-7);
                fluidState.setMoleFractionSecComp(nPhaseIdx, compIdx, 0);
            }
		   
            MiscibleMultiPhaseComposition::solve(fluidState,
                                                 paramCache,
                                                 knownPhaseIdx);
												 // ii);//,
												 // co2aqn);//,
												 //ratio
			// printf("The value of moleFracforchem2 is: %.5e\n", moleFrac[16]);
			// printf("The value of moleFrac[CO2aqIdx2] is: %.8e\n", fluidState.moleFraction(wPhaseIdx, CO2aqIdx)); // 
			// printf("The value of  priVars[CO2aqIdx2] is: %.8e\n",  priVars[CO2aqIdx]); 												 
			// fluidState.setMoleFractionSecComp(phase1Idx, CO2aqonlyIdx, fluidState.moleFraction(phase1Idx, CO2aqIdx));
			// fluidState.setMoleFraction(phase0Idx, CO2aqIdx, priVars[CO2aqIdx]);	// added by Du, co2aq > co2aq
			// printf("The value of moleFrac[CO2aqIdx3] is: %.8e\n", fluidState.moleFraction(wPhaseIdx, CO2aqIdx)); // 
			// printf("The value of  priVars[CO2aqIdx3] is: %.8e\n",  priVars[CO2aqIdx]); 	
			// printf("The value of moleFrac[CO2aqIdx]3 is: %.8e\n", moleFrac[CO2aqIdx]); //2.65687305e-06
            // if (useConstraintSolver) {
            //     MiscibleMultiPhaseComposition::solve(fluidState,
            //                                          paramCache,
            //                                          knownPhaseIdx);
            // }
            // // ... or calculated explicitly this way ...
            // // please note that we experienced some problems with un-regularized
            // // partial pressures due to their calculation from fugacity coefficients -
            // // that's why they are regularized below "within physically meaningful bounds"
            // else {
            //   Scalar partPressH2O = FluidSystem::fugacityCoefficient(fluidState_,
            //                                                         wPhaseIdx,
            //                                                         wCompIdx) * pw_;
            //   if (partPressH2O > pg_) partPressH2O = pg_;
            //   Scalar partPressNAPL = FluidSystem::fugacityCoefficient(fluidState_,
            //                                                          nPhaseIdx,
            //                                                          nCompIdx) * pn_;
            //   if (partPressNAPL > pg_) partPressNAPL = pg_;
            //   Scalar partPressAir = pg_ - partPressH2O - partPressNAPL;
			//
            //   Scalar xgn = partPressNAPL/pg_;
            //   Scalar xgw = partPressH2O/pg_;
            //   Scalar xgg = partPressAir/pg_;
			//
            //   // actually, it's nothing else than Henry coefficient
            //   Scalar xwn = partPressNAPL
            //                / (FluidSystem::fugacityCoefficient(fluidState_,
            //                                                    wPhaseIdx,nCompIdx)
            //                   * pw_);
            //   Scalar xwg = partPressAir
            //                / (FluidSystem::fugacityCoefficient(fluidState_,
            //                                                    wPhaseIdx,gCompIdx)
            //                   * pw_);
            //   Scalar xww = 1.-xwg-xwn;
			//
            //   Scalar xnn = 1.-2.e-10;
            //   Scalar xna = 1.e-10;
            //   Scalar xnw = 1.e-10;
			//
            //   fluidState_.setMoleFraction(wPhaseIdx, wCompIdx, xww);
            //   fluidState_.setMoleFraction(wPhaseIdx, gCompIdx, xwg);
            //   fluidState_.setMoleFraction(wPhaseIdx, nCompIdx, xwn);
            //   fluidState_.setMoleFraction(gPhaseIdx, wCompIdx, xgw);
            //   fluidState_.setMoleFraction(gPhaseIdx, gCompIdx, xgg);
            //   fluidState_.setMoleFraction(gPhaseIdx, nCompIdx, xgn);
            //   fluidState_.setMoleFraction(nPhaseIdx, wCompIdx, xnw);
            //   fluidState_.setMoleFraction(nPhaseIdx, gCompIdx, xna);
            //   fluidState_.setMoleFraction(nPhaseIdx, nCompIdx, xnn);
			//
            //   Scalar rhoW = FluidSystem::density(fluidState_, wPhaseIdx);
            //   Scalar rhoG = FluidSystem::density(fluidState_, gPhaseIdx);
            //   Scalar rhoN = FluidSystem::density(fluidState_, nPhaseIdx);
            //   Scalar rhoWMolar = FluidSystem::molarDensity(fluidState_, wPhaseIdx);
            //   Scalar rhoGMolar = FluidSystem::molarDensity(fluidState_, gPhaseIdx);
            //   Scalar rhoNMolar = FluidSystem::molarDensity(fluidState_, nPhaseIdx);
			//
            //   fluidState_.setDensity(wPhaseIdx, rhoW);
            //   fluidState_.setDensity(gPhaseIdx, rhoG);
            //   fluidState_.setDensity(nPhaseIdx, rhoN);
            //   fluidState_.setMolarDensity(wPhaseIdx, rhoWMolar);
            //   fluidState_.setMolarDensity(gPhaseIdx, rhoGMolar);
            //   fluidState_.setMolarDensity(nPhaseIdx, rhoNMolar);
            // }
												 
            for (int compIdx= 0; compIdx<numComponents; ++compIdx)
            {
                moleFrac[compIdx] =fluidState.moleFraction(wPhaseIdx, compIdx);
                // moleFrac[comp1Idx] =fluidState.moleFraction(nPhaseIdx, compIdx);
				//printf("The value of moleFrac1 is: %.6e\n", moleFrac[6]);
                // moleFrac[6] =1e-20*fluidState.moleFraction(wPhaseIdx, 6);
                // moleFrac[9] =1e-20*fluidState.moleFraction(wPhaseIdx, 9);
                // moleFrac[10] =1e-20*fluidState.moleFraction(wPhaseIdx, 10);
                // moleFrac[11] =1e-20*fluidState.moleFraction(wPhaseIdx, 11);	
                // moleFrac[13] =1e-20*fluidState.moleFraction(wPhaseIdx, 13);				
				// printf("The value of moleFrac2 is: %.6e\n", moleFrac[6]);
                // printf("The value of waterfugacity is: %.6e\n", fluidState.fugacityCoefficient(wPhaseIdx, compIdx));
		        // printf("The value of gasfugacity is: %.6e\n", fluidState.fugacityCoefficient(nPhaseIdx, compIdx));  
			}
			// printf("The value of moleFracforchem2 is: %.5e\n", moleFrac[16]);
		    // printf("The value of moleFrac is: %.e\n", moleFrac[comp0Idx]);
            // set the fluid state for secondary components to zero for now. will be calculated later in the chemistry


           // for (int compIdx=numComponents; compIdx<numComponents+numSecComponents; ++compIdx)
           //  {
           //      fluidState.setMoleFractionSecComp(wPhaseIdx, compIdx, 0);
           //      fluidState.setMoleFractionSecComp(nPhaseIdx, compIdx, 0);
           //  }
			

        }
        //we don't really care about the unrealistic nPhaseOnly-case, most comonents are only in the water phase
        else if (phasePresence == nPhaseOnly)
        {
            moleFrac[comp0Idx] = priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = moleFrac[comp0Idx];

            for (int compIdx = numMajorComponents; compIdx < numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
                sumMoleFracOtherComponents += moleFrac[compIdx];
            }

            moleFrac[comp1Idx] = 1 - sumMoleFracOtherComponents;

            // Set fluid state mole fractions
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState.setMoleFraction(phase1Idx, compIdx, moleFrac[compIdx]);

            // calculate the composition of the remaining phases (as
            // well as the densities of all phases). this is the job
            // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             phase1Idx);

            // set the fluid state for secondary components (only in water phase) to zero in the nPhaseOnly-case
            for (int compIdx=numComponents; compIdx<numComponents+numSecComponents; ++compIdx)
            {
                fluidState.setMoleFractionSecComp(wPhaseIdx, compIdx, 0);
                fluidState.setMoleFractionSecComp(nPhaseIdx, compIdx, 0);
            }

		}
		
        else if (phasePresence == wPhaseOnly)
        {
            // only the wetting phase is present, i.e. wetting phase
            // composition is stored explicitly.
            // extract _mass_ fractions in the non-wetting phase

            moleFrac[comp1Idx] = priVars[switchIdx];
            Scalar sumMoleFracOtherComponents = moleFrac[comp1Idx];
            for (int compIdx = numMajorComponents; compIdx < numComponents; ++compIdx)
            {
                moleFrac[compIdx] = priVars[compIdx];
//                 std::cout <<"w phase only, priVars["<< FluidSystem::componentName(compIdx)<<"]"<<priVars[compIdx]<< std::endl;

                sumMoleFracOtherComponents += moleFrac[compIdx];
            }

            moleFrac[comp0Idx] = 1 - sumMoleFracOtherComponents;

            // convert mass to mole fractions and set the fluid state
            for (int compIdx = 0; compIdx < numComponents; ++compIdx)
                fluidState.setMoleFraction(phase0Idx, compIdx, moleFrac[compIdx]);
            // set the fluid state for secondary components to zero for now. will be calculated later in the chemistry
            for (int compIdx=numComponents; compIdx<numComponents + numSecComponents; ++compIdx)
            {
                moleFrac[compIdx] = 0;
            }

//             // calculate the composition of the remaining phases (as
//             // well as the densities of all phases). this is the job
//             // of the "ComputeFromReferencePhase" constraint solver
            ComputeFromReferencePhase::solve(fluidState,
                                             paramCache,
                                             phase0Idx);

            //TODO:later use Brine_CO2 again
            Scalar XwSalinity= 0.0;
            for (int compIdx = NaIdx; compIdx<= CaIdx ; compIdx++)  //salinity = XlNa + XlCl + XlCa
            {
                if(fluidState.massFraction(wPhaseIdx, compIdx)>0)
                {
                    XwSalinity+= fluidState.massFraction(wPhaseIdx, compIdx);
                }
            }
                        Scalar xnH2O;
                        Scalar xwCO2;

             Brine_CO2::calculateMoleFractions(fluidState.temperature(),
                                        fluidState.pressure(nPhaseIdx),
                                        XwSalinity,
                                        /*knownPhaseIdx=*/-1,
                                        xwCO2,
                                        xnH2O);
            // normalize the phase compositions
            xwCO2 = std::max(0.0, std::min(1.0, xwCO2));
            xnH2O = std::max(0.0, std::min(1.0, xnH2O));

            // set the fluid state
            for (int compIdx=0; compIdx<numComponents+numSecComponents; ++compIdx)
            {
                fluidState.setMoleFraction(wPhaseIdx, compIdx, moleFrac[compIdx]);
                fluidState.setMoleFraction(nPhaseIdx, compIdx, 0);
            }
            //only CO2 and water are present in the non-wetting Phase
            fluidState.setMoleFraction(nPhaseIdx, wCompIdx, xnH2O);
            fluidState.setMoleFraction(nPhaseIdx, nCompIdx, 1-xnH2O);
        }
        // printf("The value of moleFracforchem is: %.e\n", moleFrac[comp0Idx]);
		// Scalar rhoMolar = 0.0; // added by du Declare rhoMolar outside the loop and initialize it
		Scalar rhoMolar = 0.0; // added by du Declare rhoMolar outside the loop and initialize it
		// printf("The value of moleFracforchem1 is: %.5e\n", moleFrac[16]);
		paramCache.updateAll(fluidState);

        for (int phaseIdx = 0; phaseIdx < ModelTraits::numPhases(); ++phaseIdx)
        {
            Scalar rho = FluidSystem::density(fluidState, paramCache, phaseIdx);
            Scalar mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
            // printf("The value of mu is: %.e\n", mu);
			Scalar h = ParentType::enthalpy(fluidState, paramCache, phaseIdx);
            // added by du
		    rhoMolar = FluidSystem::molarDensity(fluidState, paramCache, phase0Idx);
		    // printf("The value of rhoMolar1 is: %.8e\n", rhoMolar);
            // added by du           
			fluidState.setMolarDensity(phaseIdx, rhoMolar);
            fluidState.setDensity(phaseIdx, rho);
            fluidState.setViscosity(phaseIdx, mu);
            fluidState.setEnthalpy(phaseIdx, h);
        }
        // printf("The value of moleFracforctot1 is: %.e\n", moleFrac[CtotIdx]);
        //calculate the actual equilibrium aqueous phase composition including secondary components
//         if (phasePresence == bothPhases)
//         {
	    //   // printf("The value of moleFrac is: %.e\n", moleFrac[comp0Idx]);
           Chemistry chemistry;
		   // printf("The value of moleFracforchem is: %.e\n", moleFrac);
		   // printf("The value of moleFracforchem is: %.e\n", moleFrac[2]);
		   // printf("The value of moleFracforchem2 is: %.5e\n", moleFrac[16]);
		   // printf("reacehed here 2");
           chemistry.calculateEquilibriumChemistry(fluidState, phasePresence, moleFrac, rhoMolar);
		   
//         }
           // printf("The value of moleFracforctot2 is: %.e\n", moleFrac[CtotIdx]);
           for (int compIdx=numComponents; compIdx<numComponents + numSecComponents; ++compIdx)
           {
               fluidState.setMoleFractionSecComp(phase0Idx, compIdx, moleFrac[compIdx]);
               fluidState.setMoleFractionSecComp(phase1Idx, compIdx, 0);
           }
		   
		   Scalar knownPhaseIdx = phase0Idx;
		   MiscibleMultiPhaseComposition::solve(fluidState,
                                     paramCache,
                                     knownPhaseIdx);
        //  //all inorganic carbon in the gas phase is CO2:
           fluidState.setMoleFractionSecComp(phase1Idx, CO2aqonlyIdx, fluidState.moleFraction(phase1Idx, CO2aqIdx));
			

			// fluidState.setMoleFraction(wPhaseIdx, CtotIdx, moleFrac[CtotIdx]);
			// printf("The value of moleFracforchem is: %.e\n", moleFrac[CtotIdx]);		
			// fluidState.setMoleFraction(phase1Idx, CO2aqIdx, moleFrac[CO2aqonlyIdx]*fluidState.fugacityCoefficient(wPhaseIdx, CO2aqIdx));
			// fluidState.setMoleFraction(wPhaseIdx, HIdx, moleFrac[HIdx]);
            // all inorganic carbon in the gas phase is CO2 which is not right here
			
            // fluidState.setMoleFractionSecComp(phase1Idx, CO2aqIdx, fluidState.moleFraction(phase1Idx, CO2aqIdx));
			// fluidState.setMoleFraction(wPhaseIdx, CO2aqIdx, moleFrac[CO2aqIdx]);
			// fluidState.setMoleFraction(wPhaseIdx, HIdx, moleFrac[HIdx]);
		    // printf("The CO2aqIdx2 is: %.8e\n", moleFrac[CO2aqIdx]);
            // solidState.setVolumeFraction(gPhaseIdx, solidState.volumeFraction(gPhaseIdx)); // added by du
		    // printf("The vfofglass is: %.10f\n", solidState.volumeFraction(gPhaseIdx));
		    // paramCache.updateAll(fluidState); // added by du
		    // printf("The porosity2 is: %.10f\n", porosity_);
			
			// // added by du
            // const int knownPhaseIdx = setFirstPhaseMoleFractions ? phase0Idx : phase1Idx;
            // MiscibleMultiPhaseComposition::solve(fluidState,
            //                                      paramCache,
            //                                      knownPhaseIdx);
            // for (int compIdx= 0; compIdx<numComponents; ++compIdx)
            // {
            //     moleFrac[compIdx] =fluidState.moleFraction(wPhaseIdx, compIdx);
            //     // moleFrac[comp1Idx] =fluidState.moleFraction(nPhaseIdx, compIdx);
			// 	// printf("The value of moleFrac2 is: %.e\n", fluidState.moleFraction(knownPhaseIdx,compIdx));
            //     // printf("The value of waterfugacity is: %.6e\n", fluidState.fugacityCoefficient(wPhaseIdx, compIdx));
		    //     // printf("The value of gasfugacity is: %.6e\n", fluidState.fugacityCoefficient(nPhaseIdx, compIdx));  
			// }
		    // printf("The value of moleFrac is: %.e\n", moleFrac[comp0Idx]);
            // set the fluid state for secondary components to zero for now. will be calculated later in the chemistry

    }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const FluidState &fluidState() const
    { return fluidState_; }

    /*!
     * \brief Returns the phase state for the control-volume.
     */
    const SolidState &solidState() const
    { return solidState_; }

    /*!
     * \brief Returns the saturation of a given phase within
     *        the control volume in \f$[-]\f$.
     *
     * \param phaseIdx The phase index
     */
    Scalar saturation(int phaseIdx) const
    { return fluidState_.saturation(phaseIdx); }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar density(int phaseIdx) const
    {
        if (phaseIdx < ModelTraits::numPhases())
            return fluidState_.density(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the kinematic viscosity of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar viscosity(int phaseIdx) const
    {
        if (phaseIdx < ModelTraits::numPhases())
            return fluidState_.viscosity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the mass density of a given phase within the
     *        control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar molarDensity(int phaseIdx) const
    {
        if (phaseIdx < ModelTraits::numPhases())
            return fluidState_.molarDensity(phaseIdx);
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    }

    /*!
     * \brief Returns the effective pressure of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar pressure(int phaseIdx) const
    { return fluidState_.pressure(phaseIdx); }

    /*!
     * \brief Returns temperature inside the sub-control volume.
     *
     * Note that we assume thermodynamic equilibrium, i.e. the
     * temperature of the rock matrix and of all fluid phases are
     * identical.
     */
    Scalar temperature() const
    { return fluidState_.temperature(/*phaseIdx=*/0); }

    /*!
     * \brief Returns the effective mobility of a given phase within
     *        the control volume.
     *
     * \param phaseIdx The phase index
     */
    Scalar mobility(int phaseIdx) const
    { return mobility_[phaseIdx]; }

    /*!
     * \brief Returns the effective capillary pressure within the control volume
     *        in \f$[kg/(m*s^2)=N/m^2=Pa]\f$.
     */
    Scalar capillaryPressure() const
    { return fluidState_.pressure(FluidSystem::nPhaseIdx) - fluidState_.pressure(FluidSystem::wPhaseIdx); }

    /*!
     * \brief Returns the average porosity within the control volume.
     */
    Scalar porosity() const
    { return porosity_; }

    /*!
     * \brief Returns the permeability within the control volume.
     */
    const PermeabilityType& permeability() const
    { return permeability_; } 
	// added by du
    //{
    //    if (phaseIdx < ModelTraits::numPhases())
    //        return fluidState_.density(phaseIdx);
    //    else
    //        DUNE_THROW(Dune::InvalidStateException, "Invalid phase index " << phaseIdx);
    //}

    /*!
     * \brief Returns the diffusion coefficient
     */
    Scalar diffusionCoefficient(int phaseIdx, int compIdx) const
    {
        Scalar value;
		
		if (compIdx < phaseIdx)
            value = diffCoefficient_[phaseIdx][compIdx];
        else if (compIdx > phaseIdx)
            value = diffCoefficient_[phaseIdx][compIdx-1];
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient called for phaseIdx = compIdx");
        
		// printf("The value of diffusionCoefficient is: %.e\n", value);
        return value;
    }
   /*!
     * \brief Returns the molarity of a component in the phase
     *
     * \param phaseIdx the index of the fluid phase
     * \param compIdx the index of the component
     */
     Scalar molarity(int phaseIdx, int compIdx) const // [moles/m^3]
    { return fluidState_.molarity(phaseIdx, compIdx);}
	// molarDensity(phaseIdx)*moleFraction(phaseIdx, compIdx)

     /*!
      * \brief Returns the mass fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar massFraction(int phaseIdx, int compIdx) const
     { return fluidState_.massFraction(phaseIdx, compIdx); }

     /*!
      * \brief Returns the mole fraction of a component in the phase
      *
      * \param phaseIdx the index of the fluid phase
      * \param compIdx the index of the component
      */
     Scalar moleFraction(int phaseIdx, int compIdx) const
     { return fluidState_.moleFraction(phaseIdx, compIdx); }

//      int numSecFluidComponents()
//      { return numSecFluidComps;}

protected:
    FluidState fluidState_;
    SolidState solidState_;

private:
    void setDiffusionCoefficient_(int phaseIdx, int compIdx, Scalar d)
    {
        if (compIdx < phaseIdx)
            diffCoefficient_[phaseIdx][compIdx] = std::move(d);
        else if (compIdx > phaseIdx)
            diffCoefficient_[phaseIdx][compIdx-1] = std::move(d);
        else
            DUNE_THROW(Dune::InvalidStateException, "Diffusion coeffiecient for phaseIdx = compIdx doesn't exist");
        // printf("The value is: %.e\n", d);
	}

    Scalar refPC_;                     //!< The reference capillary pressure
    Scalar refPorosity_;               //!< The reference effective porosity within the control volume
    PermeabilityType refPermeability_; //!> The reference effective permeability within the control
    Scalar pc_;                        //!< The capillary pressure
    Scalar porosity_;                  //!< Effective porosity within the control volume
    PermeabilityType permeability_;    //!> Effective permeability within the control volume
    Scalar mobility_[ModelTraits::numPhases()]; //!< Effective mobility within the control volume
    std::array<std::array<Scalar, numComponents-1>, ModelTraits::numPhases()> diffCoefficient_;
};

} // end namespace Dumux

#endif
