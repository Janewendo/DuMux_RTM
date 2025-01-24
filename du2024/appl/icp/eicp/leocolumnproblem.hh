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

#ifndef DUMUX_LEO_COLUMN_PROBLEM_HH
#define DUMUX_LEO_COLUMN_PROBLEM_HH

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/leomin.hh>
#include <dumux/material/solidsystems/leominsolids.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/2picp/model.hh>

#include <dumux/material/binarycoefficients/brine_co2.hh>
#include <dumux/material/chemistry/biogeochemistry/leocarbonicacid.hh>

#include <appl/icp/icpspatialparams.hh>
#include <appl/icp/co2tableslaboratoryhightemp.hh>

#include "dumux/linear/seqsolverbackend.hh"

#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#define NONISOTHERMAL 0

namespace Dumux
{
template <class TypeTag>
class LEOColumnProblem;

namespace Properties
{
// Create new type tags
namespace TTag {
#if NONISOTHERMAL
struct LEOColumnTypeTag { using InheritsFrom = std::tuple<TwoPICPNI>; };
struct LEOColumnCCTpfaTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, CCTpfaModel>; };
struct LEOColumnBoxTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, BoxModel>; };
#else  LEO
struct LEOColumnTypeTag { using InheritsFrom = std::tuple<TwoPICP>; };
struct LEOColumnCCTpfaTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, CCTpfaModel>; };	
struct LEOColumnBoxTypeTag { using InheritsFrom = std::tuple<LEOColumnTypeTag, BoxModel>; };
#endif
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::LEOColumnTypeTag> { using type = Dune::YaspGrid<3>; }; 

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::LEOColumnTypeTag> { using type = LEOColumnProblem<TypeTag>; };

//Set the CO2 tables used.
SET_TYPE_PROP(LEOColumnTypeTag, CO2Tables, Dumux::ICP::CO2Tables);

// set the fluidSystem
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LEOColumnTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = GetPropType<TypeTag, Properties::CO2Tables>;
    using H2OTabulated = Components::TabulatedComponent<Components::H2O<Scalar>>;
    using type = Dumux::FluidSystems::LeoMinFluid<Scalar, CO2Tables, H2OTabulated>;
};

// set the solidSystem
template<class TypeTag>
struct SolidSystem<TypeTag, TTag::LEOColumnTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = SolidSystems::LeoMinSolidPhase<Scalar>;
};


//Set the problem chemistry
template<class TypeTag>
struct Chemistry<TypeTag, TTag::LEOColumnTypeTag>
{
    using CO2Tables = GetPropType<TypeTag, Properties::CO2Tables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using type = Dumux::LeoMinCarbonicAcid<TypeTag, CO2Tables, ModelTraits>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LEOColumnTypeTag> { using type = ICPSpatialParams<TypeTag>; };

template<class TypeTag>
struct Formulation<TypeTag, TTag::LEOColumnTypeTag>
// { static constexpr auto value = TwoPFormulation::p0s1; };
{ static constexpr auto value = TwoPFormulation::p1s0; };

}

/*!
 * \ingroup TwoPNCSecCompMinModel
 * \ingroup ImplicitTestProblems
 * \brief Problem for enzyme-induced calcium carbonate precipitation
 *  */
template <class TypeTag>
class LEOColumnProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        numComponents = FluidSystem::numComponents,

        pressureIdx = Indices::pressureIdx,
        switchIdx = Indices::switchIdx, //Saturation
		xwN2Idx = FluidSystem::N2Idx, 
        xwO2Idx = FluidSystem::O2Idx,   
        xwNaIdx = FluidSystem::NaIdx,
        xwHIdx = FluidSystem::HIdx,
        xwClIdx = FluidSystem::ClIdx,
        xwCaIdx = FluidSystem::CaIdx,
        xwFe2Idx = FluidSystem::Fe2Idx,
		xwCO2aqIdx = FluidSystem::CO2aqIdx,
        xwMgIdx = FluidSystem::MgIdx,
        xwAlIdx = FluidSystem::AlIdx,
        xwKIdx = FluidSystem::KIdx,
        xwMnIdx = FluidSystem::MnIdx,
        xwSiO2Idx = FluidSystem::SiO2Idx,
        xwHPO4Idx = FluidSystem::HPO4Idx,
        xwTiOH4Idx = FluidSystem::TiOH4Idx,		
        xwFe2totalIdx = FluidSystem::Fe2totalIdx,
		xwCO2aqtotalIdx = FluidSystem::CO2aqtotalIdx,
        xwHtotalIdx = FluidSystem::HtotalIdx,
        xwFe3Idx = FluidSystem::Fe3Idx,	
        xwCO3Idx = FluidSystem::CO3Idx,	
		xwHCO3Idx = FluidSystem::HCO3Idx,	
		xwOHIdx = FluidSystem::OHIdx,	
		
        phiFerrohydriteIdx = numComponents +1,
		phiGlassIdx = numComponents,
		phiProtoImogoliteIdx = numComponents +2,
		phiBirnessiteIdx = numComponents +3,
		phiHydroxyapatiteIdx = numComponents +4,
		phiSepioliteIdx = numComponents +5,
		
#if NONISOTHERMAL
        temperatureIdx = Indices::temperatureIdx,
        energyEqIdx = Indices::energyEqIdx,
#endif

        //Indices of the components
        wCompIdx = FluidSystem::wCompIdx,
        nCompIdx = FluidSystem::nCompIdx,
        O2Idx = FluidSystem::O2Idx,
        NaIdx = FluidSystem::NaIdx,
        HIdx = FluidSystem::HIdx,
        CaIdx = FluidSystem::CaIdx,
        Fe2Idx = FluidSystem::Fe2Idx,
		CO2aqIdx = FluidSystem::CO2aqIdx,
        MgIdx = FluidSystem::MgIdx,
        AlIdx = FluidSystem::AlIdx,
        KIdx = FluidSystem::KIdx,
        MnIdx = FluidSystem::MnIdx,
        SiO2Idx = FluidSystem::SiO2Idx,
        HPO4Idx = FluidSystem::HPO4Idx,
        TiOH4Idx = FluidSystem::TiOH4Idx,
        ClIdx = FluidSystem::ClIdx,	

        CO3Idx = FluidSystem::CO3Idx,
        HCO3Idx = FluidSystem::HCO3Idx,
        OHIdx = FluidSystem::OHIdx,

        //Index of the primary component of G and L phase
        conti0EqIdx = Indices::conti0EqIdx,
        conti1EqIdx = Indices::conti1EqIdx,
        // Phase State
        nPhaseOnly = Indices::secondPhaseOnly,
		// added by du
		wPhaseOnly = Indices::firstPhaseOnly,
        bothPhases = Indices::bothPhases,

        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
    };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using BoundaryTypes = GetPropType<TypeTag, Properties::BoundaryTypes>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Chemistry = GetPropType<TypeTag, Properties::Chemistry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView; // added
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace; // added

public:
    LEOColumnProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);

        name_  = getParam<std::string>("Problem.Name");

        //initial values
        densityW_ = getParam<Scalar>("Initial.initDensityW");
        initPressure_ = getParam<Scalar>("Initial.initPressure");
		
        initxwCO2aq_ = getParam<Scalar>("Initial.initxwCO2");
        initxwO2_ = getParam<Scalar>("Initial.initxwO2");
        initxwNa_ = getParam<Scalar>("Initial.initxwNa");
        initxwH_ = getParam<Scalar>("Initial.initxwH");
        initxwCl_ = getParam<Scalar>("Initial.initxwCl");
        initxwCa_ = getParam<Scalar>("Initial.initxwCa");
        initxwFe2_ = getParam<Scalar>("Initial.initxwFe2");
        initxwMg_ = getParam<Scalar>("Initial.initxwMg");
        initxwAl_ = getParam<Scalar>("Initial.initxwAl");
        initxwK_ = getParam<Scalar>("Initial.initxwK");
        initxwMn_ = getParam<Scalar>("Initial.initxwMn");
        initxwSiO2_ = getParam<Scalar>("Initial.initxwSiO2");
        initxwHPO4_ = getParam<Scalar>("Initial.initxwHPO4");
        initxwTiOH4_ = getParam<Scalar>("Initial.initxwTiOH4");
	

        initGlass_ = getParam<Scalar>("Initial.initGlass");
        // initCalcite_ = getParam<Scalar>("Initial.initCalcite");
        initFerrohydrite_ = getParam<Scalar>("Initial.initFerrohydrite");
		initProtoImogolite_ = getParam<Scalar>("Initial.initProtoImogolite");
		initBirnessite_ = getParam<Scalar>("Initial.initBirnessite");
		initHydroxyapatite_ = getParam<Scalar>("Initial.initHydroxyapatite");
		initSepiolite_ = getParam<Scalar>("Initial.initSepiolite");

        initTemperature_ = getParam<Scalar>("Initial.initTemperature");

        xwNaCorr_ = getParam<Scalar>("Initial.xwNaCorr");

        //injection values
        injQ_ = getParam<Scalar>("Injection.injVolumeflux");

        injO2_ = getParam<Scalar>("Injection.injO2");
        injCO2aq_ = getParam<Scalar>("Injection.injCO2");
        injH_ = getParam<Scalar>("Injection.injH");
        injNa_ = getParam<Scalar>("Injection.injNa");
        injCl_ = getParam<Scalar>("Injection.injCl");
        injCa_ = getParam<Scalar>("Injection.injCa");
        injFe2_ = getParam<Scalar>("Injection.injFe2");
		injMg_ = getParam<Scalar>("Injection.injMg");
		injAl_ = getParam<Scalar>("Injection.injAl");
		injK_ = getParam<Scalar>("Injection.injK");
		injMn_ = getParam<Scalar>("Injection.injMn");
		injSiO2_ = getParam<Scalar>("Injection.injSiO2");
		injHPO4_ = getParam<Scalar>("Injection.injHPO4");
		injTiOH4_ = getParam<Scalar>("Injection.injTiOH4");


        injNaCorr_ = getParam<Scalar>("Injection.injNaCorr");
        injTemperature_ = getParam<Scalar>("Injection.injTemperature");
        injPressure_ = getParam<Scalar>("Injection.injPressure");

        numInjections_ = getParam<int>("Injection.numInjections");
        injectionParameters_ = getParam<std::string>("Injection.InjectionParamFile");

        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;
        permeability_.resize(gridGeometry->gridView().size(codim));
        calcium_.resize(gridGeometry->gridView().size(codim));
        fe2_.resize(gridGeometry->gridView().size(codim));
		
		newdt_.resize(gridGeometry->gridView().size(codim));
		
        std::ifstream injectionData;
        std::string row;
        injectionData.open( injectionParameters_); // open the Injection data file
        if (not injectionData.is_open())
        {
            std::cerr << "\n\t -> Could not open file '"
                    << injectionParameters_
                    << "'. <- \n\n\n\n";
            exit(1) ;
        }
        // int tempType = 0;
		double tempType = 0.0;

        // print file to make sure it is the right file
        std::cout << "Read file: " << injectionParameters_ << " ..." << std::endl << std::endl;
        while(!injectionData.eof())
        {
            getline(injectionData, row);
            std::cout << row << std::endl;
        }
        injectionData.close();

      //read data from file
        injectionData.open(injectionParameters_);

        while(!injectionData.eof())
        {
            getline(injectionData, row);

	  
            if(row == "InjectionTypes")
            {
                // printf("reach here 3");
				getline(injectionData, row);
                while(row != "#")
                {
                    
					if (row != "#")
                        {
                        std::istringstream ist(row);
                        ist >> tempType;
                        injType_.push_back(tempType);
//                      std::cout << "size of injType: "<<injType_.size() << std::endl;
                        }
                    getline(injectionData, row);
                }
            }
        }

        injectionData.close();

//      check the injection data against the number of injections specified in the parameter file
        if (injType_.size() != numInjections_)
        {
            std::cerr <<  "numInjections from the parameterfile and the number of injection types specified in the injection data file do not match!"
                    <<"\n numInjections from parameter file = "<<numInjections_
                    <<"\n numInjTypes from injection data file = "<<injType_.size()
                    <<"\n Abort!\n";
            exit(1) ;
        }

#if NONISOTHERMAL
        FluidSystem::init(/*startTemp=*/295.15, /*endTemp=*/445.15, /*tempSteps=*/151,
             /*startPressure=*/1e4, /*endPressure=*/1e6, /*pressureSteps=*/500);
#else
        FluidSystem::init(/*startTemp=*/initTemperature_ -5.0, /*endTemp=*/initTemperature_ +5.0, /*tempSteps=*/5,
             /*startPressure=*/1e4, /*endPressure=*/1e6, /*pressureSteps=*/500);
#endif
    }

    void setTime( Scalar time )
    {
        time_ = time;
    }

    void setTimeStepSize( Scalar timeStepSize )
    {
        timeStepSize_ = timeStepSize;
    }

    void setEpisodeIdx( int epiIdx )
    {
        episodeIdx_ = epiIdx;
    }


    double injectionType(int episodeIdx)
    {
        return injType_[episodeIdx];
    }

   /*!
    * \name Problem parameters
    */


   /*!
    * \brief The problem name.
    *
    * This is used as a prefix for files generated by the simulation.
    */
//    const char *name() const
    const std::string name() const
    { return name_; }

 #if !NONISOTHERMAL
   /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 25 degrees Celsius.
    */
    Scalar temperature() const
    {
        return initTemperature_; //
    };
#endif

    // \}

   /*!
    * \name Boundary conditions
    */
    // \{

    /*!
    * \brief Specifies which kind of boundary condition should be
    *        used for which equation on a given boundary segment.
    */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const

    {
        BoundaryTypes bcTypes;
	    // // set all other as Neumann boundaries

	    if(globalPos[dim - 1]> this->gridGeometry().bBoxMax()[dim - 1] - eps_) 
	    {

			bcTypes.setAllNeumann();
	    	bcTypes.setDirichlet(xwCO2aqIdx,conti0EqIdx+xwCO2aqIdx); // du changed
	    	bcTypes.setDirichlet(xwO2Idx,conti0EqIdx+xwO2Idx); // du changed
	    	bcTypes.setDirichlet(xwHIdx,conti0EqIdx+xwHIdx); // du changed
		}

		
        // if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
		else if (globalPos[dim - 1]<= eps_)
	    { 
            bcTypes.setAllNeumann();
	    	bcTypes.setDirichlet(xwCO2aqIdx,conti0EqIdx+xwCO2aqIdx); // du changed
	    	bcTypes.setDirichlet(xwO2Idx,conti0EqIdx+xwO2Idx); // du changed
	    	bcTypes.setDirichlet(xwHIdx,conti0EqIdx+xwHIdx); // du changed
	    }

	    return bcTypes;
    }
   /*!
    * \brief Evaluate the boundary conditions for a dirichlet
    *        boundary segment.
    */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {

		if(globalPos[dim - 1]> this->gridGeometry().bBoxMax()[dim - 1] - eps_)
        {
			return top_(globalPos);
		}
		else if(globalPos[dim - 1] <= eps_)
        {
			return bottom_(globalPos);
		}
    }

   /*!
    * \brief Evaluate the initial value for a control volume.
    *
    * \param globalPos The global position
    *
    * For this method, the \a values parameter stores primary
    * variables.
    */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return initial_(globalPos);
    }

   /*!
    * \brief Evaluate the boundary conditions for a Neumann
    *        boundary segment.
    *
    * For this method, the \a values parameter stores the mass flux
    * in normal direction of each component. Negative values mean
    * influx.
    *
    * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
    */

    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf
						) const
    {  
        NumEqVector values(0.0);
        const auto xMax = this->gridGeometry().bBoxMax()[0];
        const auto& ipGlobal = scvf.ipGlobal(); 
		const auto& globalPos = scvf.center();
		const auto& volVars = elemVolVars[scvf.insideScvIdx()];

		 if (globalPos[dim - 1]> this->gridGeometry().bBoxMax()[dim - 1] - eps_)	
         {
             
			 // Scalar waterFlux =  -1*injType_[episodeIdx_]/3600; // injQ_/1; //[m/s]  0; //
             Scalar waterFlux =  -injQ_/1;
             values[conti0EqIdx + wCompIdx] = waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW_/FluidSystem::molarMass(wCompIdx);
             values[conti0EqIdx + nCompIdx] = waterFlux * volVars.moleFraction(conti0EqIdx, nCompIdx) * densityW_/FluidSystem::molarMass(wCompIdx);
			 values[conti0EqIdx + xwHIdx] = waterFlux * injH_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwCO2aqIdx] = waterFlux * injCO2aq_ * densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwO2Idx] = waterFlux * injO2_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwFe2Idx] = waterFlux * injFe2_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwCaIdx] = waterFlux * injCa_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwKIdx] = waterFlux * injK_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwMgIdx] = waterFlux * injMg_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwAlIdx] = waterFlux * injAl_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwMnIdx] = waterFlux * injMn_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwHPO4Idx] = waterFlux * injHPO4_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwTiOH4Idx] = waterFlux * injTiOH4_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwSiO2Idx] = waterFlux * injSiO2_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwNaIdx] = waterFlux * (injNa_+injNaCorr_)* densityW_ /FluidSystem::molarMass(wCompIdx); // FluidSystem::molarMass(NaIdx);
			 values[conti0EqIdx + xwClIdx] = waterFlux * injCl_ * densityW_ /FluidSystem::molarMass(wCompIdx);// FluidSystem::molarMass(NaIdx);               //NaCl ---> mol Cl = mol Na
	
             static const Scalar gasdirichletPressure = initPressure_;
             const auto& volVars = elemVolVars[scvf.insideScvIdx()];
             const auto& fluxVarsCache = elemFluxVarsCache[scvf];		   
			 
		 	 auto d = ipGlobal - element.geometry().center();
		     d /= d.two_norm2();
		 	 auto upwindgasTerm = volVars.mobility(nCompIdx);
			 
             const auto tij = vtmv(scvf.unitOuterNormal(), volVars.permeability(), d); // permeability m2 / m
			 
		 	 const Scalar densityN = FluidSystem::CO2::gasDensity(297.15, 1.01325e5);		 
	         
			 const auto gasFlux = 1.0*upwindgasTerm*tij*(-gasdirichletPressure+volVars.pressure(nCompIdx)+densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*(volVars.permeability()/tij)); // ipGlobal[dimWorld-1]); // [1/(Pa s)] [m2 / m] Pa = [m / s]

			  values[conti0EqIdx + nCompIdx] += gasFlux * volVars.moleFraction(conti1EqIdx, nCompIdx) *densityN/FluidSystem::molarMass(nCompIdx);
		 	  values[conti0EqIdx + wCompIdx] += gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx) *densityN/FluidSystem::molarMass(nCompIdx);
			  values[conti0EqIdx + xwCO2aqIdx] +=  gasFlux * volVars.moleFraction(conti1EqIdx, xwCO2aqIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	  values[conti0EqIdx + xwO2Idx] += gasFlux * volVars.moleFraction(conti1EqIdx, xwO2Idx) *densityN/FluidSystem::molarMass(nCompIdx);
	 }	
		 else if (ipGlobal[0] <= eps_)

         {
	 	  // set a fixed pressure on the right side of the domain
	 	  // q = -K (delta h/delta x) h=pw/(phow g) +z
             static const Scalar gasdirichletPressure = initPressure_;
             const auto& volVars = elemVolVars[scvf.insideScvIdx()];
			 
             const auto& fluxVarsCache = elemFluxVarsCache[scvf];		   
		 
		 	 auto d = ipGlobal - element.geometry().center();

		     d /= d.two_norm2();

		 	 auto upwindTerm = volVars.mobility(wCompIdx);

		 	 auto upwindgasTerm = volVars.mobility(nCompIdx);
			 
             const auto tij = vtmv(scvf.unitOuterNormal(), volVars.permeability(), d); // permeability m2 / m
		 	 const Scalar densityW = FluidSystem::H2O::liquidDensity(297.15, volVars.pressure(wCompIdx));//1.01325e5);
		 	 const Scalar densityN = FluidSystem::CO2::gasDensity(297.15, volVars.pressure(nCompIdx));// 1.01325e5);
			 const auto waterFlux = -1.0*upwindTerm*tij*(gasdirichletPressure- volVars.pressure(nCompIdx)+densityW_ *this->spatialParams().gravity(ipGlobal)[dimWorld-1]*(volVars.permeability()/tij));//ipGlobal[dimWorld-1]); 
			 const auto gasFlux = -1.0*upwindgasTerm*tij*(gasdirichletPressure-volVars.pressure(nCompIdx)+densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*(volVars.permeability()/tij)); //ipGlobal[dimWorld-1]); // [1/(Pa s)] [m2 / m] Pa = [m / s]

		 	 values[conti0EqIdx + nCompIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, nCompIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	 values[conti0EqIdx + wCompIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	 values[conti0EqIdx + xwCO2aqIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, xwCO2aqIdx)*densityN /FluidSystem::molarMass(nCompIdx);
 		 	 values[conti0EqIdx + xwO2Idx] = gasFlux * volVars.moleFraction(conti1EqIdx, xwO2Idx)*densityN /FluidSystem::molarMass(nCompIdx);
			 if (waterFlux >= 0)
			 {			 
             values[conti0EqIdx + wCompIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW_ /FluidSystem::molarMass(wCompIdx);		 				
             values[conti0EqIdx + nCompIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, nCompIdx) * densityW_/FluidSystem::molarMass(wCompIdx);
			 values[conti0EqIdx + xwCO2aqIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwCO2aqIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwHIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwHIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwFe2Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwFe2Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
	         
			 // values[conti0EqIdx + xwCO2aqIdx] += waterFlux * (volVars.moleFraction(conti0EqIdx, xwCO2aqtotalIdx)-volVars.moleFraction(conti0EqIdx, xwCO3Idx)-volVars.moleFraction(conti0EqIdx, xwHCO3Idx))* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 // values[conti0EqIdx + xwHIdx] += waterFlux * (volVars.moleFraction(conti0EqIdx, xwHtotalIdx)+volVars.moleFraction(conti0EqIdx, xwOHIdx)+2*volVars.moleFraction(conti0EqIdx, xwCO3Idx)+volVars.moleFraction(conti0EqIdx, xwHCO3Idx)-2*volVars.moleFraction(conti0EqIdx, xwFe3Idx))* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 // values[conti0EqIdx + xwFe2Idx] += waterFlux * (volVars.moleFraction(conti0EqIdx, xwFe2totalIdx)-volVars.moleFraction(conti0EqIdx, xwFe3Idx))* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
	         // printf("The value of volVars.moleFraction(conti0EqIdx, xwFe2totalIdx-xwFe3) is: %.6e\n", volVars.moleFraction(conti0EqIdx, xwFe2totalIdx)-volVars.moleFraction(conti0EqIdx, xwFe3Idx));
			 // printf("The value of volVars.moleFraction(conti0EqIdx, xwFe2Idx) is: %.6e\n", volVars.moleFraction(conti0EqIdx, xwFe2Idx));
			 // printf("The value of values[conti0EqIdx + xwCO2aqIdx] is: %.6e\n", values[conti0EqIdx + xwCO2aqIdx]);
	         // printf("The value of values[conti0EqIdx + xwHIdx] is: %.6e\n", values[conti0EqIdx + xwHIdx]);
	         // printf("The value of values[conti0EqIdx + xwFe2Idx] is: %.6e\n", values[conti0EqIdx + xwFe2Idx]);		
			 
			 values[conti0EqIdx + xwO2Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwO2Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwCaIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwCaIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwKIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwKIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwMgIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwMgIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwAlIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwAlIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwMnIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwMnIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwHPO4Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwHPO4Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwTiOH4Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwTiOH4Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwSiO2Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwSiO2Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwCaIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwCaIdx)* densityW_ /FluidSystem::molarMass(CaIdx);
			 values[conti0EqIdx + xwNaIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwNaIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); // FluidSystem::molarMass(NaIdx);
			 values[conti0EqIdx + xwClIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwClIdx) * densityW_ /FluidSystem::molarMass(wCompIdx);// FluidSystem::molarMass(NaIdx);               //NaCl ---> mol Cl = mol Na
			 }          
         }
		
         else
         {
           values = 0.0; //mol/m²/s
         }
         return values;
    }
   /*!
    * \name Volume terms
    */
    // \{

   /*!
    * \brief Evaluate the source term for all phases within a given
    *        sub-control-volume.
    *
    * This is the method for the case where the source term is
    * potentially solution dependent and requires some quantities that
    * are specific to the fully-implicit method.
    *
    * \param values The source and sink values for the conservation equations in units of
    *                 \f$ [ \textnormal{unit of conserved quantity} / (m^3 \cdot s )] \f$
    * \param element The finite element
    * \param fvGeometry The finite-volume geometry
    * \param elemVolVars All volume variables for the element
    * \param scv The subcontrolvolume
    *
    * For this method, the \a values parameter stores the conserved quantity rate
    * generated or annihilate per volume unit. Positive values mean
    * that the conserved quantity is created, negative ones mean that it vanishes.
    * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
    */
    NumEqVector source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        Chemistry chemistry;
		source = std::get<0>(chemistry.reactionSource(elemVolVars[scv],
                         timeStepSize_));		
				
        source = 0.0; //added by du
        return source;
    }

	double getnewdt(const SolutionVector& curSol)
    {
        Chemistry chemistry;
   
        // Iterate over elements in the grid
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
            auto fvGeometry = localView(this->gridGeometry());  // Geometry for finite volume
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
				VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
				const auto dofIdxGlobal = scv.dofIndex();  // Global degree of freedom index
				newdt_[dofIdxGlobal] = std::get<1>(chemistry.reactionSource(volVars, timeStepSize_));
			}		  
        }
		double minnewdt_ = *std::min_element(newdt_.begin(), newdt_.end());		
        return minnewdt_;  // Return the updated vector
		
    }


	
   /*!
    * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
    */

    const std::vector<Scalar>& getPermeability()
    {
        return permeability_;
    }
    const std::vector<Scalar>& getCalcium()
    {
        return calcium_;
    }

    const std::vector<Scalar>& getFe2()
    {
        return fe2_;
    }



    void updateVtkOutput(const SolutionVector& curSol)
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
                const auto dofIdxGlobal = scv.dofIndex();
                permeability_[dofIdxGlobal] = volVars.permeability();
                calcium_[dofIdxGlobal] = volVars.moleFraction(0,CaIdx)* volVars.molarDensity(0) * FluidSystem::molarMass(CaIdx);
                fe2_[dofIdxGlobal] = volVars.moleFraction(0,Fe2Idx)* volVars.molarDensity(0) * FluidSystem::molarMass(Fe2Idx);
 
            }
        }
    }

    void setGridVariables(std::shared_ptr<GridVariables> gv)
    { gridVariables_ = gv; }

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables initial_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);

		priVars.setState(bothPhases);
		
        priVars[pressureIdx] = initPressure_;
		
		priVars[switchIdx] = 0.3;
        priVars[xwHIdx] = initxwH_;
		priVars[xwCO2aqIdx] = initxwCO2aq_;
        priVars[xwO2Idx] = initxwO2_;
	    priVars[xwClIdx] = initxwCl_ + xwClCorr_;
        priVars[xwFe2Idx] = initxwFe2_;
        priVars[xwMgIdx] = initxwMg_;
        priVars[xwCaIdx] = initxwCa_;
        priVars[xwKIdx] = initxwK_;
        priVars[xwHPO4Idx] = initxwHPO4_;
		priVars[xwNaIdx] = initxwNa_;
        priVars[xwAlIdx] = initxwAl_; 
        priVars[xwSiO2Idx] = initxwSiO2_;
        priVars[xwTiOH4Idx] = initxwTiOH4_; 
        priVars[xwMnIdx] = initxwMn_;

        priVars[phiGlassIdx] = initGlass_; // [m^3/m^3]
        priVars[phiFerrohydriteIdx] = initFerrohydrite_; // [m^3/m^3]
		priVars[phiProtoImogoliteIdx] = initProtoImogolite_; // [m^3/m^3]
		priVars[phiBirnessiteIdx] = initBirnessite_; // [m^3/m^3]
		priVars[phiHydroxyapatiteIdx] = initHydroxyapatite_; // [m^3/m^3]
		priVars[phiSepioliteIdx] = initSepiolite_; // [m^3/m^3]
#if NONISOTHERMAL
        priVars[temperatureIdx] = initTemperature_;
#endif
        return priVars;
    }
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    PrimaryVariables top_(const GlobalPosition &globalPos) const
    {
	    // injProcess == 1;	
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = initPressure_;
        priVars.setState(bothPhases);
        priVars[xwHIdx] = initxwH_;
		priVars[xwCO2aqIdx] = initxwCO2aq_;
        priVars[xwO2Idx] = initxwO2_;
#if NONISOTHERMAL
        priVars[temperatureIdx] = initTemperature_;
#endif
        return priVars;
    }
    PrimaryVariables bottom_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        priVars[pressureIdx] = initPressure_;
        priVars.setState(bothPhases);
        priVars[xwO2Idx] = initxwO2_;
        priVars[xwHIdx] = initxwH_; 
		priVars[xwCO2aqIdx] = initxwCO2aq_; 
#if NONISOTHERMAL
        priVars[temperatureIdx] = initTemperature_;
#endif
        return priVars;
    }
    /*!
        * \brief Returns the molality of NaCl (mol NaCl / kg water) for a given mole fraction
        *
        * \param XwNaCl the XwNaCl [kg NaCl / kg solution]
        */
    static Scalar massTomoleFrac_(Scalar XwNaCl)
    {
        const Scalar Mw = FluidSystem::molarMass(wCompIdx);  // 18.015e-3; /* molecular weight of water [kg/mol] */
        const Scalar Ms = FluidSystem::molarMass(NaIdx) + FluidSystem::molarMass(ClIdx); // 58.44e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = XwNaCl;
        /* XwNaCl: conversion from mass fraction to mol fraction */
        const Scalar xwNaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return xwNaCl;
    }

    static constexpr Scalar eps_ = 1e-6;

    Scalar initPressure_;
    Scalar densityW_;

    Scalar densityW;
    Scalar densityN;
	
    Scalar initxwCO2aq_;
	Scalar initxwO2_;
    Scalar initxwNa_;
    Scalar initxwH_;
    Scalar initxwCa_;
    Scalar initxwCl_;	
    Scalar initxwFe2_;
    Scalar initxwK_;
	Scalar initxwMg_;
    Scalar initxwAl_;
    Scalar initxwMn_;
    Scalar initxwHPO4_;
    Scalar initxwTiOH4_;
    Scalar initxwSiO2_;	
    Scalar xwNaCorr_;
    Scalar xwClCorr_;

    Scalar initGlass_;
    Scalar initFerrohydrite_;
	Scalar initProtoImogolite_;
	Scalar initBirnessite_;
	Scalar initHydroxyapatite_;
	Scalar initSepiolite_;
	
    Scalar initPorosity_;
    Scalar initTemperature_;

    Scalar injQ_;

    Scalar injCO2aq_;
	Scalar injO2_;
    Scalar injNa_;
    Scalar injCa_;
    Scalar injCl_;
    Scalar injH_;
    Scalar injFe2_;
    Scalar injK_;
	Scalar injMg_;
    Scalar injAl_;
    Scalar injMn_;
    Scalar injHPO4_;
    Scalar injTiOH4_;
    Scalar injSiO2_;	
	
    Scalar injNaCorr_;
    Scalar injTemperature_;
    Scalar injPressure_;

    int numInjections_;
    std::string injectionParameters_;

    std::vector<double> injType_;	
    std::string name_;

    std::vector<Scalar> permeability_;
    std::vector<Scalar> calcium_;
    std::vector<Scalar> fe2_;
    std::vector<Scalar> newdt_;

    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    int episodeIdx_ = 0;
    std::shared_ptr<GridVariables> gridVariables_;
};
} //end namespace

#endif




