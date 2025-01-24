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
#include <dune/grid/uggrid.hh>
// #include <dune/grid/albertagrid.hh>
// #include <dune/alugrid/grid.hh>

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
// struct Grid<TypeTag, TTag::LEOColumnTypeTag> { using type = Dune::UGGrid<3>; }; //3D
// struct Grid<TypeTag, TTag::LEOColumnTypeTag> { using type = Dune::UGGrid<2>; }; //2D
struct Grid<TypeTag, TTag::LEOColumnTypeTag> { using type = Dune::YaspGrid<3>; }; //for 2D
// struct Grid<TypeTag, TTag::LEOColumnTypeTag> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };

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
		// xwwaterIdx = FluidSystem::BrineIdx,
		xwN2Idx = FluidSystem::N2Idx, //xwTCIdx = FluidSystem::TCIdx,//FluidSystem::TCaqIdx,// 
        xwO2Idx = FluidSystem::O2Idx,   
        xwNaIdx = FluidSystem::NaIdx,
        xwHIdx = FluidSystem::HIdx,
        xwClIdx = FluidSystem::ClIdx,
        xwCaIdx = FluidSystem::CaIdx,
        xwFe2Idx = FluidSystem::Fe2Idx,
        // xwFe3Idx = FluidSystem::Fe3Idx,
		xwCO2aqIdx = FluidSystem::CO2aqIdx,
		// xwCO2aqIdx = FluidSystem::CO2aqIdx,
        xwMgIdx = FluidSystem::MgIdx,
        xwAlIdx = FluidSystem::AlIdx,
        xwKIdx = FluidSystem::KIdx,
        xwMnIdx = FluidSystem::MnIdx,
        xwSiO2Idx = FluidSystem::SiO2Idx,
        xwHPO4Idx = FluidSystem::HPO4Idx,
        xwTiOH4Idx = FluidSystem::TiOH4Idx,		
		// xwTCaqIdx = FluidSystem::TCaqIdx,
        // xwUreaIdx = FluidSystem::UreaIdx,
        // xwTNHIdx = FluidSystem::TNHIdx,
        // xwUreaseIdx = FluidSystem::UreaseIdx,
        // phiImmUreaseIdx = numComponents,
        // phiCalciteIdx = numComponents+1 ,
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
		// waterIdx = FluidSystem::BrineIdx,
		// N2Idx = FluidSystem::N2Idx,// TCIdx = FluidSystem::TCIdx,
        O2Idx = FluidSystem::O2Idx,
        NaIdx = FluidSystem::NaIdx,
        HIdx = FluidSystem::HIdx,
        CaIdx = FluidSystem::CaIdx,
        Fe2Idx = FluidSystem::Fe2Idx,
		CO2aqIdx = FluidSystem::CO2aqIdx,
        // TNHIdx = FluidSystem::TNHIdx,
        MgIdx = FluidSystem::MgIdx,
        AlIdx = FluidSystem::AlIdx,
        KIdx = FluidSystem::KIdx,
        MnIdx = FluidSystem::MnIdx,
        SiO2Idx = FluidSystem::SiO2Idx,
        HPO4Idx = FluidSystem::HPO4Idx,
        TiOH4Idx = FluidSystem::TiOH4Idx,
        ClIdx = FluidSystem::ClIdx,	
		
        // Fe3Idx = FluidSystem::Fe3Idx,
        // // NH4Idx = FluidSystem::NH4Idx,
        CO3Idx = FluidSystem::CO3Idx,
        HCO3Idx = FluidSystem::HCO3Idx,
        // CO2aqIdx = FluidSystem::CO2aqIdx,
		// // CO2gIdx = FluidSystem::CO2gIdx,

        OHIdx = FluidSystem::OHIdx,

        //Index of the primary component of G and L phase
        conti0EqIdx = Indices::conti0EqIdx,
        conti1EqIdx = Indices::conti1EqIdx,
        // conti0EqIdx = Indices::conti0EqIdx,
        // conti1EqIdx = Indices::conti1EqIdx,
        // conti0EqIdx = Indices::conti0EqIdx,
        // conti1EqIdx = Indices::conti1EqIdx,
		// conti0EqIdx = Indices::conti0EqIdx,
        // conti1EqIdx = Indices::conti1EqIdx,
		// conti0EqIdx = Indices::conti0EqIdx,
        // conti1EqIdx = Indices::conti1EqIdx,
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
        // initSaturation_ = getParam<Scalar>("Initial.initgasSaturation");
        // initgasPressure_ = getParam<Scalar>("Initial.initgasPressure");
		
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
	
        // initxwTNH_ = getParam<Scalar>("Initial.initxwTNH");;
        initGlass_ = getParam<Scalar>("Initial.initGlass");
        // initCalcite_ = getParam<Scalar>("Initial.initCalcite");
        initFerrohydrite_ = getParam<Scalar>("Initial.initFerrohydrite");
		initProtoImogolite_ = getParam<Scalar>("Initial.initProtoImogolite");
		initBirnessite_ = getParam<Scalar>("Initial.initBirnessite");
		initHydroxyapatite_ = getParam<Scalar>("Initial.initHydroxyapatite");
		initSepiolite_ = getParam<Scalar>("Initial.initSepiolite");
        // initImmUrease_ = getParam<Scalar>("Initial.initImmUrease");
        initTemperature_ = getParam<Scalar>("Initial.initTemperature");

        xwNaCorr_ = getParam<Scalar>("Initial.xwNaCorr");
        // xwClCorr_ = getParam<Scalar>("Initial.xwClCorr");

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
		// injUrea_ = getParam<Scalar>("Injection.injUrea");
        // injTNH_ = getParam<Scalar>("Injection.injTNH");
        // injEnzymeSource_= getParam<Scalar>("Injection.injEnzymeSource");

        injNaCorr_ = getParam<Scalar>("Injection.injNaCorr");
        injTemperature_ = getParam<Scalar>("Injection.injTemperature");
        injPressure_ = getParam<Scalar>("Injection.injPressure");

        numInjections_ = getParam<int>("Injection.numInjections");
        injectionParameters_ = getParam<std::string>("Injection.InjectionParamFile");

        unsigned int codim = GetPropType<TypeTag, Properties::GridGeometry>::discMethod == DiscretizationMethod::box ? dim : 0;
        permeability_.resize(gridGeometry->gridView().size(codim));
        calcium_.resize(gridGeometry->gridView().size(codim));
        // urea_.resize(gridGeometry->gridView().size(codim));
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
	    // printf("reach here 1");
        while(!injectionData.eof())
        {
            getline(injectionData, row);
	//	    printf("reach here 2");
	  
            if(row == "InjectionTypes")
            {
                // printf("reach here 3");
				getline(injectionData, row);
                while(row != "#")
                {
                    // printf("reach here 4");
					if (row != "#")
                        {
                        std::istringstream ist(row);
                        ist >> tempType;
                        injType_.push_back(tempType);
	                    // printf("reach here 5");
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

    // int injectionType(int episodeIdx)
    // {
    //     return injType_[episodeIdx];
    // }

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
    // {
    //     BoundaryTypes bcTypes;
	// 
    //     Scalar zmax = this->gridGeometry().bBoxMax()[dim - 1];
    //     bcTypes.setAllNeumann();
    //     if(globalPos[dim - 1] > zmax - eps_)
    //         bcTypes.setAllDirichlet();
	// 
    //     return bcTypes;
    // }
    {
        BoundaryTypes bcTypes;
	    // // set all other as Neumann boundaries
        // bcTypes.setAllNeumann();   // this will make it use the neumann condition 	
		// set the left of the domain (with the global position in "0 = x" direction as a Dirichlet boundary
        // printf("The value of dim is: %.e\n", dim); // added by du
	    if(globalPos[dim - 1]> this->gridGeometry().bBoxMax()[dim - 1] - eps_) // if(globalPos[1]<= eps_) bottom // if (globalPos[0] < eps_) left
	    {
			// printf("top boundary.\n");
			// printf("The value of globalPos[0] is: %.e\n", globalPos[0]);
			// printf("dim.\n");
			// printf("The value of dim is: %.e\n", dim);
			// printf("dimWorld.\n");
			// printf("The value of dimWorld is: %.e\n", dimWorld);
    		// bcTypes.setDirichlet(pressureIdx,conti1EqIdx); // du changed
            // bcTypes.setDirichlet(switchIdx,conti0EqIdx); // du changed
			// bcTypes.setNeumann(conti0EqIdx);
			// bcTypes.setNeumann(conti1EqIdx);
		    // bcTypes.setAllDirichlet();
			bcTypes.setAllNeumann();
	    	bcTypes.setDirichlet(xwCO2aqIdx,conti0EqIdx+xwCO2aqIdx); // du changed
	    	bcTypes.setDirichlet(xwO2Idx,conti0EqIdx+xwO2Idx); // du changed
	    	bcTypes.setDirichlet(xwHIdx,conti0EqIdx+xwHIdx); // du changed
	    	// bcTypes.setDirichlet(pressureIdx,conti1EqIdx); // du changed
	    	// bcTypes.setDirichlet(xwCO2aqIdx,conti1EqIdx); // du changed
	    	// bcTypes.setDirichlet(xwO2Idx,conti1EqIdx); // du changed
	    	// bcTypes.setDirichlet(xwO2Idx,conti1EqIdx); 
			// bcTypes.setDirichlet(switchIdx,conti0EqIdx); 
			// bcTypes.setNeumann(conti1EqIdx);
            // bcTypes.setNeumann(conti0EqIdx);
		}
	    // set the top of the domain as Neumann boundaries
        // else if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)  0.055586<TAN(10/180)
        
	    // 3D 
	    // else if (globalPos[2] > globalPos[0] * 0.05558 + 1.0 - eps_)
	    // { 
	    // 	// printf("top boundary.\n");
		// 	// printf("The value of globalPos[1] is: %.2lf\n", globalPos[1]);
		// 	bcTypes.setDirichlet(pressureIdx,conti1EqIdx); // du changed
	    // 	// printf("top boundary.\n");
        //     bcTypes.setNeumann(conti0EqIdx);
		// 	// bcTypes.setAllDirichlet();
	    // }
		
        // if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_)
		else if (globalPos[dim - 1]<= eps_)
	    { 
            bcTypes.setAllNeumann();
            // bcTypes.setNeumann(conti0EqIdx);
	    	// bcTypes.setDirichlet(pressureIdx,conti1EqIdx); // du changed
	    	bcTypes.setDirichlet(xwCO2aqIdx,conti0EqIdx+xwCO2aqIdx); // du changed
	    	bcTypes.setDirichlet(xwO2Idx,conti0EqIdx+xwO2Idx); // du changed
	    	bcTypes.setDirichlet(xwHIdx,conti0EqIdx+xwHIdx); // du changed
	    	// bcTypes.setDirichlet(xwCO2aqIdx,conti1EqIdx); // du changed
	    	// bcTypes.setDirichlet(xwO2Idx,conti1EqIdx); // du changed
			// bcTypes.setDirichlet(switchIdx,conti0EqIdx); 
			// bcTypes.setNeumann(conti1EqIdx);
		    // bcTypes.setAllDirichlet();
			// bcTypes.setAllNeumann();
	    }

	    return bcTypes;
    }
   /*!
    * \brief Evaluate the boundary conditions for a dirichlet
    *        boundary segment.
    */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // return top_(globalPos);
		// return initial_(globalPos);
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
//    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
//    {
//        NumEqVector values(0.0);
//
//        Scalar waterFlux = injQ_/4.640675e-4; //[m/s]
//        //4.640675e-4 m^2 = 0.957" diameter column cross-sectional area.
//
//        int injProcess = injType_[episodeIdx_];
//
////         negative values for injection
//        // if(globalPos[1]<= eps_)
//		if(globalPos[1] > gridGeometry().bBoxMax()[1] - eps_)
//        {
//           //basic rinse injection (injProcess == -1 )
//           values[conti0EqIdx + wCompIdx] = -waterFlux * 996/FluidSystem::molarMass(wCompIdx);
//           values[conti0EqIdx + nCompIdx] = -waterFlux * injTC_*996 /FluidSystem::molarMass(nCompIdx);
//           values[conti0EqIdx + xwCaIdx] = 0;
//           values[conti0EqIdx + xwUreaIdx] = 0;
//           values[conti0EqIdx + xwUreaseIdx] = 0;
//           values[conti0EqIdx + xwTNHIdx] = 0;
//           values[conti0EqIdx + phiCalciteIdx] = 0;
//           values[conti0EqIdx + xwFe2Idx] = 0;
//           values[conti0EqIdx + phiFerrohydriteIdx] = 0;
//           values[conti0EqIdx + phiImmUreaseIdx] = 0;
//           values[conti0EqIdx + xwNaIdx] = -waterFlux * (injNa_+injNaCorr_) /FluidSystem::molarMass(NaIdx);
//           values[conti0EqIdx + xwClIdx] = -waterFlux *injNa_ /FluidSystem::molarMass(NaIdx);               //NaCl ---> mol Cl = mol Na
//#if NONISOTHERMAL
//            values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
//                                    injTemperature_, injPressure_,
//                                   (injNa_
//                                   +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
//#endif
//
//            if (injProcess == -1)       // rinse,
//            {
//                //only NH4Cl
//                values[conti0EqIdx + xwTNHIdx] += -waterFlux * injTNH_ /FluidSystem::molarMass(TNHIdx);
//                values[conti0EqIdx + xwClIdx] += -waterFlux * injTNH_ /FluidSystem::molarMass(TNHIdx);
//#if NONISOTHERMAL
//                values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
//                                    injTemperature_, injPressure_,
//                                    (injTNH_/FluidSystem::molarMass(TNHIdx)*FluidSystem::molarMass(ClIdx)
//                                    +injNa_
//                                    +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
//#endif
//            }
//
//           else if (injProcess <-89)       // no injection
//           {
//            values = 0.0; //mol/m²/s
//           }
//
//           else if (injProcess == 1)              //ca-rich injection: ca and urea injected Na (pH) and Cl are different(CaCl2)
//           {
//               values[conti0EqIdx + wCompIdx] = - waterFlux * 0.8716 * densityW_ /FluidSystem::molarMass(wCompIdx);       //TODO 0.8716 check factor!!!
//               values[conti0EqIdx + nCompIdx] = - waterFlux * injTC_ * densityW_ /FluidSystem::molarMass(nCompIdx);
//               values[conti0EqIdx + xwCaIdx] = - waterFlux * injCa_/FluidSystem::molarMass(CaIdx);
//               values[conti0EqIdx + xwFe2Idx] = - waterFlux * injFe2_/FluidSystem::molarMass(Fe2Idx);			   
//               values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
//               values[conti0EqIdx + xwClIdx] += - waterFlux * 2 * injCa_/FluidSystem::molarMass(CaIdx);
//
//#if NONISOTHERMAL
//            values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
//                                    injTemperature_, injPressure_,
//                                   (injCa_
//                                   +2*injCa_ /FluidSystem::molarMass(CaIdx)*FluidSystem::molarMass(ClIdx)
//                                   +injNa_
//                                   +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
//#endif
//       }
//              else if(injProcess == 2)              //injection of Urease and Jack Bean Meal
//       {
//                values[conti0EqIdx + xwUreaseIdx] = -waterFlux * injEnzymeSource_ /FluidSystem::molarMass(UreaseIdx);
//       }
//               else
//               {
//                      DUNE_THROW(Dune::InvalidStateException, "Invalid injection process " << injProcess);
//               }
//        }
//        else
//        {
//               values = 0.0; //mol/m²/s
//        }
//        return values;
//    }

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
		// no-flow at the right boundary and bottom and left boundary
        // if ((ipGlobal[0] > this->gridGeometry().bBoxMax()[0] - eps_) or ipGlobal[1] < (ipGlobal[0] * 0.05558 + eps_) or ipGlobal[1] > ipGlobal[0] * 0.05558 + 1.0 - eps_)
        // {
		// 	return values;
		// }
		// 
		// 
		// // else if (ipGlobal[1] > ipGlobal[0] * 0.05558 + 1.0 - eps_)
		// else
		 if (globalPos[dim - 1]> this->gridGeometry().bBoxMax()[dim - 1] - eps_)	
 		 // if (ipGlobal[0]> this->gridGeometry().bBoxMax()[0] - eps_) // (ipGlobal[0] <= eps_) //  // bottom boundary, left boundary if (ipGlobal[0] <= eps_)
         {
			 // printf("The value of episodeIdx_ is: %e\n", episodeIdx_);
			 // printf("neumann boundary.\n");
		     // printf("The value of ipGlobal[0] is: %.e\n", ipGlobal[0]);
             
			 Scalar waterFlux =  -1*injType_[episodeIdx_]/3600; // injQ_/1; //[m/s]  0; //
             // Scalar waterFlux =  -1*injQ_/1; //[m/s]  0; //
			 
             // printf("The value of waterFlux1 is: %e\n", waterFlux);
			 // printf("The value of episodeIdx_ is: %e\n", episodeIdx_);
			 // printf("The value of injType_[episodeIdx_] is: %e\n", injType_[episodeIdx_]);
			 // Scalar waterFlux = 0; //[m/s] tested fine, the Sw increase with depth
             //4.640675e-4 m^2 = 0.957" diameter column cross-sectional area.
 		     // values[conti1EqIdx] = - waterFlux * 996/FluidSystem::molarMass(wCompIdx);
             // int injProcess = injType_[episodeIdx_];
 		 
             //         negative values for injection upward
             // if(globalPos[1]<= eps_)      
   		       //basic rinse injection (injProcess == -1 )
             values[conti0EqIdx + wCompIdx] = waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW_/FluidSystem::molarMass(wCompIdx);
             values[conti0EqIdx + nCompIdx] = waterFlux * volVars.moleFraction(conti0EqIdx, nCompIdx) * densityW_/FluidSystem::molarMass(wCompIdx);
             // printf("The value of conti0EqIdx + nCompIdx is: %.8e\n", values[conti0EqIdx + nCompIdx]);
			 // printf("The value of conti0EqIdx + wCompIdx is: %.10e\n",  values[conti0EqIdx + wCompIdx]);
			 // printf("The value of waterFlux1 is: %.8e\n", waterFlux);
			 // printf("The value of volVars.moleFraction(conti0EqIdx, wCompIdx)1 is: %.8e\n", volVars.moleFraction(conti0EqIdx, wCompIdx));
			 // printf("The value of densityW1 is: %.8e\n", densityW_);
			 // printf("The value of FluidSystem::molarMass(wCompIdx)1 is: %.8e\n", FluidSystem::molarMass(wCompIdx));
			 // printf("The value of values[conti0EqIdx + wCompIdx]1 is: %.8e\n", values[conti0EqIdx + wCompIdx]);
			 // printf("The value of 1waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW/FluidSystem::molarMass(wCompIdx) is: %.8e\n", waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW_/FluidSystem::molarMass(wCompIdx));

	         // printf("The value of values[conti0EqIdx + wCompIdx] is: %.2lf\n", values[conti0EqIdx + wCompIdx]);
	         // printf("The value of values[conti1EqIdx + nCompIdx] is: %.e\n", values[conti1EqIdx + nCompIdx]);
	         // // printf("The value of values[conti0EqIdx + nCompIdx] is: %.2lf\n", values[conti0EqIdx + nCompIdx]);
	         // printf("The value of values[conti1EqIdx + wCompIdx] is: %.2lf\n", values[conti1EqIdx + wCompIdx]);

	         // values[conti0EqIdx + nCompIdx] = waterFlux * injTC_* densityW_ /FluidSystem::molarMass(nCompIdx); //du deleted
             // values[conti0EqIdx + xwCaIdx] = 0;
             // // values[conti0EqIdx + xwUreaIdx] = 0;
             // // values[conti0EqIdx + xwUreaseIdx] = 0;
             // // values[conti0EqIdx + xwTNHIdx] = 0;
             // // values[conti0EqIdx + phiCalciteIdx] = 0;
             // values[conti0EqIdx + xwFe2Idx] = 0;
             // // values[conti0EqIdx + phiFerrohydriteIdx] = 0;
             // // values[conti0EqIdx + phiImmUreaseIdx] = 0;
             // injH_ = 1.80686E-05*(injH_/1.80686E-05-1E-14/(injH_/1.80686E-05)-2*4.52168E-7*4.69029E-11*injCO2aq_/1.80686E-05/(injH_/1.80686E-05)/(injH_/1.80686E-05)-4.52168E-7*injCO2aq_/1.80686E-05/(injH_/1.80686E-05));
		     //injCO2aq_ = 1.80686E-05*(injCO2aq_/1.80686E-05+4.52168E-7*injCO2aq_/1.80686E-05/(injH_/1.80686E-05)+4.52168E-7*4.69029E-11*injCO2aq_/1.80686E-05/(injH_/1.80686E-05)/(injH_/1.80686E-05)); // give to S_w

			 values[conti0EqIdx + xwHIdx] = waterFlux * (1.80686E-05*(injH_/1.80686E-05-1E-14/(injH_/1.80686E-05)-2*4.52168E-7*4.69029E-11*injCO2aq_/1.80686E-05/(injH_/1.80686E-05)/(injH_/1.80686E-05)-4.52168E-7*injCO2aq_/1.80686E-05/(injH_/1.80686E-05)))* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwCO2aqIdx] = waterFlux * (1.80686E-05*(injCO2aq_/1.80686E-05+4.52168E-7*injCO2aq_/1.80686E-05/(injH_/1.80686E-05)+4.52168E-7*4.69029E-11*injCO2aq_/1.80686E-05/(injH_/1.80686E-05)/(injH_/1.80686E-05)))* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
		     // values[conti0EqIdx + xwCO2aqIdx] = waterFlux * injCO2aq_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 // values[conti0EqIdx + xwHIdx] = waterFlux * injH_* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
	         // values[conti0EqIdx + xwCO2aqIdx] = waterFlux * volVars.moleFraction(conti0EqIdx, CO2aqIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 // values[conti0EqIdx + xwHIdx] = waterFlux * volVars.moleFraction(conti0EqIdx, HIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);

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

			 // values[conti0EqIdx + xwCaIdx] = waterFlux * 0.0428* densityW_ /FluidSystem::molarMass(CaIdx);

			 values[conti0EqIdx + xwNaIdx] = waterFlux * (injNa_+injNaCorr_)* densityW_ /FluidSystem::molarMass(wCompIdx); // FluidSystem::molarMass(NaIdx);
             // printf("The value of values[conti0EqIdx + xwNaIdx] is: %.2lf\n", values[conti0EqIdx + xwNaIdx]);
	        
			 values[conti0EqIdx + xwClIdx] = waterFlux * injCl_ * densityW_ /FluidSystem::molarMass(wCompIdx);// FluidSystem::molarMass(NaIdx);               //NaCl ---> mol Cl = mol Na
 // #if  NONISOTHERMAL
 //             values[energyEqIdx] = waterFlux*Brine::liquidEnthalpy(
 //                                      injTemperature_, injPressure_,
 //                                     (injNa_
 //                                     +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
 // #end if 	
             static const Scalar gasdirichletPressure = initPressure_;
             const auto& volVars = elemVolVars[scvf.insideScvIdx()];
             const auto& fluxVarsCache = elemFluxVarsCache[scvf];		   
			 
		 	 auto d = ipGlobal - element.geometry().center();
		     d /= d.two_norm2();
		 	 auto upwindgasTerm = volVars.mobility(nCompIdx);
			 
             const auto tij = vtmv(scvf.unitOuterNormal(), volVars.permeability(), d); // permeability m2 / m
			 
		 	 const Scalar densityN = FluidSystem::CO2::gasDensity(297.15, 1.01325e5);		 
	         
			 const auto gasFlux = 1.0*upwindgasTerm*tij*(-gasdirichletPressure+volVars.pressure(nCompIdx)+densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*(volVars.permeability()/tij)); // ipGlobal[dimWorld-1]); // [1/(Pa s)] [m2 / m] Pa = [m / s]
			 // const auto gasFlux = 0.0;
			 
			 // printf("The value of upwindgasTerm is: %.10e\n",  upwindgasTerm);
			 // printf("The value of tij is: %.10e\n",  tij);
			 // printf("The value of -gasdirichletPressure+volVars.pressure(nCompIdx) is: %.10e\n",  -gasdirichletPressure+volVars.pressure(nCompIdx));
			 // printf("The value of densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*ipGlobal[dimWorld-1] is: %.10e\n",  densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*ipGlobal[dimWorld-1]);
		 	 // printf("The value of volVars.pressure(nCompIdx) is: %.8e\n", volVars.pressure(nCompIdx));
			 // printf("The value of gasdirichletPressure is: %.8e\n", gasdirichletPressure);
			 // printf("The value of gasFlux1 is: %.8e\n", gasFlux);
			 // if (gasFlux >= 0)
			 // {
			  values[conti0EqIdx + nCompIdx] += gasFlux * volVars.moleFraction(conti1EqIdx, nCompIdx) *densityN/FluidSystem::molarMass(nCompIdx);
		 	  values[conti0EqIdx + wCompIdx] += gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx) *densityN/FluidSystem::molarMass(nCompIdx);
              // values[conti0EqIdx + wCompIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx) *densityN/FluidSystem::molarMass(nCompIdx) + waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW/FluidSystem::molarMass(wCompIdx);
		 	  values[conti0EqIdx + xwCO2aqIdx] +=  gasFlux * volVars.moleFraction(conti1EqIdx, xwCO2aqIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	  values[conti0EqIdx + xwO2Idx] += gasFlux * volVars.moleFraction(conti1EqIdx, xwO2Idx) *densityN/FluidSystem::molarMass(nCompIdx);
			 // }
		 	  // printf("The value of xwCO2aqIdx is: %.8e\n", volVars.moleFraction(conti1EqIdx, xwCO2aqIdx));
			  // printf("The value of xwO2Idx is: %.8e\n", volVars.moleFraction(conti1EqIdx, xwO2Idx));
			 // if (gasFlux < 0)
			 // {
			  // values[conti0EqIdx + nCompIdx] += gasFlux * volVars.moleFraction(conti1EqIdx, nCompIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	  // values[conti0EqIdx + wCompIdx] += gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx) *densityN/FluidSystem::molarMass(nCompIdx);
		 	  // values[conti0EqIdx + xwCO2aqIdx] +=  gasFlux * initxwCO2aq_ *densityN_/FluidSystem::molarMass(nCompIdx);
 		 	  // values[conti0EqIdx + xwO2Idx] += gasFlux * initxwO2_ *densityN_/FluidSystem::molarMass(nCompIdx);
			 // }
			 // printf("The value of conti0EqIdx + wCompIdx is: %.10e\n",  values[conti0EqIdx + wCompIdx]);
			 // printf("The value of volVars.moleFraction(conti1EqIdx, wCompIdx is: %.10e\n",  volVars.moleFraction(conti1EqIdx, wCompIdx));
			 // printf("The value of densityN is: %.10e\n",  densityN);
			 // printf("The value of FluidSystem::molarMass(nCompIdx) is: %.10e\n",  FluidSystem::molarMass(nCompIdx));
			 // printf("The value of conti0EqIdx + wCompIdx is: %.10e\n",  gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx) *densityN/FluidSystem::molarMass(nCompIdx));
             // // // Check the current CO2 concentration
             // double currentCO2molefraction = volVars.moleFraction(conti1EqIdx, xwCO2aqIdx);
             // // double currentnmolefraction = volVars.moleFraction(conti1EqIdx, nCompIdx);
             // // // double currentwConcentration = volVars.moleFraction(conti1EqIdx, wCompIdx);
             // double currentO2molefraction = volVars.moleFraction(conti1EqIdx, xwO2Idx);
             // // 
             // if (currentCO2molefraction < initxwCO2aq_) 
			 // {
             //         // If the current concentration is less than the target, increase flux to restore concentration
             //         double molefractionCO2aqDifference = initxwCO2aq_ - currentCO2molefraction;
             //         
             //         // Adjust the flux proportionally to the difference
             //         values[conti1EqIdx + xwCO2aqIdx] = gasFlux * molefractionCO2aqDifference / FluidSystem::molarMass(nCompIdx);
             // 
             //         // printf("CO2 concentration below target, adjusting flux to compensate: %.8e\n", values[conti1EqIdx + xwCO2aqIdx]);
             // } 
		     // else 
		     // {
             //         // If concentration is at or above the target, no CO2 flux adjustment is needed
             //         values[conti1EqIdx + xwCO2aqIdx] = 0;  // No additional CO2 flux
             // }
			 // 
             // 
             // if (currentO2molefraction < initxwO2_) 
			 // {
             //         // If the current concentration is less than the target, increase flux to restore concentration
             //         double molefractionO2Difference = initxwO2_ - currentO2molefraction;
             //         
             //         // Adjust the flux proportionally to the difference
             //         values[conti1EqIdx + xwO2Idx] = gasFlux * molefractionO2Difference / FluidSystem::molarMass(nCompIdx);
             // 
             //         // printf("CO2 concentration below target, adjusting flux to compensate: %.8e\n", values[conti1EqIdx + xwCO2aqIdx]);
             // } 
		     // else 
		     // {
             //         // If concentration is at or above the target, no CO2 flux adjustment is needed
             //         values[conti1EqIdx + xwO2Idx] = 0;  // No additional CO2 flux
             // }
			 // 
             // if (currentnmolefraction < initxwO2_) 
			 // {
             //         // If the current concentration is less than the target, increase flux to restore concentration
             //         double molefractionO2Difference = initxwO2 - currentO2molefraction;
             //         
             //         // Adjust the flux proportionally to the difference
             //         values[conti1EqIdx + xwO2Idx] = gasFlux * molefractionO2Difference / FluidSystem::molarMass(nCompIdx);
             // 
             //         // printf("CO2 concentration below target, adjusting flux to compensate: %.8e\n", values[conti1EqIdx + xwCO2aqIdx]);
             //     } 
			 //   else 
			 //   {
             //         // If concentration is at or above the target, no CO2 flux adjustment is needed
             //         values[conti1EqIdx + xwO2Idx] = 0;  // No additional CO2 flux
             //     }
             // } 			 
 		     // printf("The value of values[conti1EqIdx + nCompIdx] is: %.8e\n", volVars.moleFraction(conti1EqIdx, nCompIdx));
 		     // printf("The value of values[conti1EqIdx + wCompIdx] is: %.8e\n", volVars.moleFraction(conti1EqIdx, wCompIdx));	
 		     // printf("The value of values[conti1EqIdx + xwCO2aqIdx] is: %.8e\n", volVars.moleFraction(conti1EqIdx, xwCO2aqIdx));		
 		     // printf("The value of values[conti1EqIdx + xwO2Idx] is: %.8e\n", volVars.moleFraction(conti1EqIdx, xwO2Idx));
	 }	
			 // else if (ipGlobal[1] > ipGlobal[0] * 0.05558 + 1.0 - eps_)
        // {
        //    // printf("The value of ipGlobal[1] is: %.2lf\n", ipGlobal[1]);
		//    Scalar waterFlux = injQ_/4.640675e-4; //[m/s]
        //    //4.640675e-4 m^2 = 0.957" diameter column cross-sectional area.
		// 
        //    // int injProcess = injType_[episodeIdx_];
		// 
        //    //         negative values for injection upward
        //    // if(globalPos[1]<= eps_)      
  		//    //basic rinse injection (injProcess == -1 )
        //    values[conti1EqIdx] = waterFlux * 996/FluidSystem::molarMass(wCompIdx);
    	//    
        //    // values[conti0EqIdx + wCompIdx] = waterFlux * 996/FluidSystem::molarMass(wCompIdx);
        //    // values[conti0EqIdx + nCompIdx] = waterFlux * injTC_*996 /FluidSystem::molarMass(nCompIdx); //du deleted
        //    // values[conti0EqIdx + xwCaIdx] = 0;
        //    // values[conti0EqIdx + xwUreaIdx] = 0;
        //    // values[conti0EqIdx + xwUreaseIdx] = 0;
        //    // values[conti0EqIdx + xwTNHIdx] = 0;
        //    // values[conti0EqIdx + phiCalciteIdx] = 0;
        //    // values[conti0EqIdx + xwFe2Idx] = 0;
        //    // values[conti0EqIdx + phiFerrohydriteIdx] = 0;
        //    // values[conti0EqIdx + phiImmUreaseIdx] = 0;
        //    // values[conti0EqIdx + xwNaIdx] = waterFlux * (injNa_+injNaCorr_) /FluidSystem::molarMass(NaIdx);
        //    // values[conti0EqIdx + xwClIdx] = waterFlux *injNa_ /FluidSystem::molarMass(NaIdx);               //NaCl ---> mol Cl = mol Na
#if NONI// SOTHERMAL
        //    values[energyEqIdx] = waterFlux*Brine::liquidEnthalpy(
        //                             injTemperature_, injPressure_,
        //                            (injNa_
        //                            +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)						
#endif  // 
        // }  
  // 
		// //
        // //     if (injProcess == -1)       // rinse,
        // //     {
        // //         //only NH4Cl
        // //         values[conti0EqIdx + xwTNHIdx] += -waterFlux * injTNH_ /FluidSystem::molarMass(TNHIdx);
        // //         values[conti0EqIdx + xwClIdx] += -waterFlux * injTNH_ /FluidSystem::molarMass(TNHIdx);
#if NONI// // SOTHERMAL
        // //         values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
        // //                             injTemperature_, injPressure_,
        // //                             (injTNH_/FluidSystem::molarMass(TNHIdx)*FluidSystem::molarMass(ClIdx)
        // //                             +injNa_
        // //                             +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
#endif  // // 
        // //     }
		// // 
        // //    else if (injProcess <-89)       // no injection
        // //    {
        // //     values = 0.0; //mol/m²/s
        // //    }
		// // 
        // //    else if (injProcess == 1)              //ca-rich injection: ca and urea injected Na (pH) and Cl are different(CaCl2)
        // //    {
        // //        values[conti0EqIdx + wCompIdx] = - waterFlux * 0.8716 * densityW_ /FluidSystem::molarMass(wCompIdx);       //TODO 0.8716 check factor!!!
        // //        values[conti0EqIdx + nCompIdx] = - waterFlux * injTC_ * densityW_ /FluidSystem::molarMass(nCompIdx);
        // //        values[conti0EqIdx + xwCaIdx] = - waterFlux * injCa_/FluidSystem::molarMass(CaIdx);
        // //        values[conti0EqIdx + xwFe2Idx] = - waterFlux * injFe2_/FluidSystem::molarMass(Fe2Idx);			   
        // //        values[conti0EqIdx + xwUreaIdx] = - waterFlux * injUrea_ /FluidSystem::molarMass(UreaIdx);
        // //        values[conti0EqIdx + xwClIdx] += - waterFlux * 2 * injCa_/FluidSystem::molarMass(CaIdx);
		// // 
#if NONI// // SOTHERMAL
        // //     values[energyEqIdx] = -waterFlux*Brine::liquidEnthalpy(
        // //                             injTemperature_, injPressure_,
        // //                            (injCa_
        // //                            +2*injCa_ /FluidSystem::molarMass(CaIdx)*FluidSystem::molarMass(ClIdx)
        // //                            +injNa_
        // //                            +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
#endif  // // 
        // //    }
        // //       else if(injProcess == 2)              //injection of Urease and Jack Bean Meal
        // //       {
        // //         values[conti0EqIdx + xwUreaseIdx] = -waterFlux * injEnzymeSource_ /FluidSystem::molarMass(UreaseIdx);
        // //       }
        // //       else
        // //       {
        // //               DUNE_THROW(Dune::InvalidStateException, "Invalid injection process " << injProcess);
        // //       }
	    // }
		// 
        // else
		 // else if (ipGlobal[0] <= eps_)
		 else if (ipGlobal[0] <= eps_)
		 // if (ipGlobal[0]> this->gridGeometry().bBoxMax()[0] - eps_) //(ipGlobal[0] <= eps_)
         {
	 	  // set a fixed pressure on the right side of the domain
	 	  // q = -K (delta h/delta x) h=pw/(phow g) +z
             static const Scalar gasdirichletPressure = initPressure_;
             const auto& volVars = elemVolVars[scvf.insideScvIdx()];
			 
             const auto& fluxVarsCache = elemFluxVarsCache[scvf];		   
	         // if (volVars.pressure(nCompIdx) > 1.013251e5)
		 	 // {
			 // printf("The value of volVars.pressure(nCompIdx) is: %10e\n", volVars.pressure(nCompIdx));
             // printf("The value of waterFlux is: %10e\n", waterFlux);// 3.088461e-06
			 // }			 
		 	 auto d = ipGlobal - element.geometry().center();
 		 	 // printf("The value of d is: %e\n", d); // 4.678068e-310
		 	 // printf("The value of d is: %e\n", d);
		     d /= d.two_norm2();
		 	 // printf("d.\n");
		 	 // printf("The value of d is: %e\n", d);		
			 // printf("The value of volVars.pressure(nCompIdx) is: %8e\n", volVars.pressure(nCompIdx));			 
             // evaluate the pressure gradient
             // GlobalPosition d(0.0);
             // for (const auto& scv : scvs(fvGeometry))
             // {
             //     const auto xIp = scv.dofPosition()[0];
             //     auto tmp = fluxVarsCache.gradN(scv.localDofIndex());
             //     // tmp *= xIp > xMax - eps_ ? gasdirichletPressure
             //     tmp *= xIp < eps_ ? gasdirichletPressure
             //                       : elemVolVars[scv].pressure(nCompIdx);
             //     d += tmp;
             // }	
		 	 // printf("d.\n");
		 	 // printf("The value of d is: %.2lf\n", d);				
             // auto upwindTerm = useMoles ? volVars.molarDensity(FluidSystem::N2Idx) : volVars.density(FluidSystem::N2Idx);
             // auto upwindTerm = volVars.density(wCompIdx); // 1000 kg/m³ 
	         // upwindTerm *= volVars.mobility(wCompIdx); // s*m/kg mass flux per unit pressure gradient
		 	 auto upwindTerm = volVars.mobility(wCompIdx);
		 	 // Kr/viscosity = 1/(Pa s)
		 	 // This one has value, around 0.0112428678
		 	 // printf("upwindTerm.\n");
		 	 // printf("The value of upwindTerm is: %.2lf\n", upwindTerm);
		 	 auto upwindgasTerm = volVars.mobility(nCompIdx);
			 
             const auto tij = vtmv(scvf.unitOuterNormal(), volVars.permeability(), d); // permeability m2 / m
		 	 // printf("tij.\n");
		 	 // printf("The value of tij is: %10e\n", tij);
			 // printf("The value of volVars.permeability() is: %10e\n", volVars.permeability());
	 	     // scvf.unitOuterNormal() unit normal vector to a surface or boundary (such as a control volume face) in the domain
	 	     // d represents a directional component, such as the distance vector between points or elements in the domain
             // // vector-tensor and vector-vector multiplication 
	 	     // permeability 2e-10, tij 3.2e-8, 1/0.0125 = 160
			 
	 	     // const auto tij = vtmv(scvf.unitOuterNormal(), volVars.permeability(FluidSystem::N2Idx), d);
	         // 
	 	     // const auto phaseFlux = -1.0*upwindTerm*tij*(dirichletPressure - volVars.pressure());
             // const auto phaseFlux = -1.0*upwindTerm*tij*(volVars.pressure(FluidSystem::N2Idx)+volVars.pressure(FluidSystem::H2OIdx)-waterdirichletPressure - volVars.pressure(N2Idx));
	         // pw-pw = pa-pc - (pa-pc) >> pa - pa
		 	 const Scalar densityW = FluidSystem::H2O::liquidDensity(297.15, volVars.pressure(wCompIdx));//1.01325e5);
		 	 const Scalar densityN = FluidSystem::CO2::gasDensity(297.15, volVars.pressure(nCompIdx));// 1.01325e5);
			 // printf("The value of densityW is: %e\n", densityW);
			 // printf("The value of densityN is: %e\n", densityN);
		  // printf("The value of densityW is: %e\n", densityW);
		  // printf("The value of gravity is: %e\n", this->spatialParams().gravity(ipGlobal)[dimWorld-1]);			 
	         // const auto waterFlux = -1.0*upwindTerm*tij*(gasdirichletPressure-volVars.pressure(nCompIdx)+densityW*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*0.005); // [1/(Pa s)] [m2 / m] Pa = [m / s]
	         
			 const auto waterFlux = -1.0*upwindTerm*tij*(gasdirichletPressure- volVars.pressure(nCompIdx)+densityW_ *this->spatialParams().gravity(ipGlobal)[dimWorld-1]*(volVars.permeability()/tij));//ipGlobal[dimWorld-1]); 
	         // const auto waterFlux = 0.0;
			 
			 // const auto waterFlux = 1.0*upwindTerm*tij*113; // [1/(Pa s)] [m2 / m] Pa = [m / s] 113 is pho g h 
	         // if (volVars.pressure(nCompIdx) > 1.013251e5)
		 	 // {
			 // // printf("The value of volVars.pressure(nCompIdx) is: %10e\n", volVars.pressure(nCompIdx));
		 	 // printf("The value of upwindTerm is: %.8e\n", upwindTerm);
			 // printf("The value of tij is: %.8e\n", tij);
			 // printf("The value of gasdirichletPressure-volVars.pressure(nCompIdx) is: %.8e\n", gasdirichletPressure-volVars.pressure(nCompIdx));
			 // printf("The value of densityW*this->spatialParams().gravity(ipGlobal)[dimWorld-1] is: %.8e\n", densityW*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*(volVars.permeability()/tij)); // *ipGlobal[dimWorld-1] );
		     // printf("The value of ipGlobal[dimWorld-1] is: %.8e\n", ipGlobal[dimWorld-1] );
             // printf("The value of waterFlux is: %10e\n", waterFlux);// 3.088461e-06
			 // }
		 	 // printf("The value of waterFlux2 is: %10e\n", waterFlux);
		 
	         // const auto gasFlux = 1.0*upwindgasTerm*tij*(-gasdirichletPressure+volVars.pressure(nCompIdx)); //+densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]); // [1/(Pa s)] [m2 / m] Pa = [m / s]
		 	 // const Scalar densityN = FluidSystem::CO2::gasDensity(297.15, 1.01325e5);		 
	         // const auto gasFlux = -1.0*upwindgasTerm*tij*(gasdirichletPressure-volVars.pressure(nCompIdx)+densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*0.005); // [1/(Pa s)] [m2 / m] Pa = [m / s]
	         
			 const auto gasFlux = -1.0*upwindgasTerm*tij*(gasdirichletPressure-volVars.pressure(nCompIdx)+densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1]*(volVars.permeability()/tij)); //ipGlobal[dimWorld-1]); // [1/(Pa s)] [m2 / m] Pa = [m / s]
		 	 // const auto gasFlux = 0.0;
			 
			 // printf("The value of upwindgasTerm is: %.8e\n", upwindgasTerm);
			 // printf("The value of tij is: %.8e\n", tij);
			 // printf("The value of gasdirichletPressure-volVars.pressure(nCompIdx) is: %.8e\n", gasdirichletPressure-volVars.pressure(nCompIdx));
			 // printf("The value of densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1] is: %.8e\n", densityN*this->spatialParams().gravity(ipGlobal)[dimWorld-1] );
		 	 // // printf("The value of volVars.pressure2(nCompIdx) is: %.8e\n", volVars.pressure(nCompIdx));
			 // // printf("The value of gasdirichletPressure2 is: %.8e\n", gasdirichletPressure);
			 // printf("The value of gasFlux2 is: %.8e\n", gasFlux);
			 // printf("The value of waterFlux2 is: %.8e\n", waterFlux);
		 	 // // printf("gasFlux.\n"); // 1.948533e+00
		 	 // printf("The value of gasFlux is: %e\n", gasFlux);
		     // printf("The value of gasdirichletPressure-volVars.pressure(nCompIdx) is: %e\n", gasdirichletPressure-volVars.pressure(nCompIdx));
		 	 // printf("volVars.pressure(nCompIdx).\n");
		 	 // printf("The value of volVars.pressure(nCompIdx) is: %e\n", volVars.pressure(nCompIdx));			
             // set Neumann bc values
             // Scalar waterFlux = injQ_; //[m/s]
			 // if (gasFlux >= 0)
			 // {
		 	 values[conti0EqIdx + nCompIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, nCompIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	 values[conti0EqIdx + wCompIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	 values[conti0EqIdx + xwCO2aqIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, xwCO2aqIdx)*densityN /FluidSystem::molarMass(nCompIdx);
 		 	 values[conti0EqIdx + xwO2Idx] = gasFlux * volVars.moleFraction(conti1EqIdx, xwO2Idx)*densityN /FluidSystem::molarMass(nCompIdx);
			 // }
		 	 // printf("The value of xwCO2aqIdx is: %.8e\n", volVars.moleFraction(conti0EqIdx, xwCO2aqIdx));
			 // printf("The value of xwO2Idx is: %.8e\n", volVars.moleFraction(conti0EqIdx, xwO2Idx));
			 // printf("The value of xwwIdx is: %.8e\n", gasFlux);
			 // printf("The value of xwwIdx is: %.8e\n", volVars.moleFraction(conti1EqIdx, wCompIdx));
			 // printf("The value of xwwIdx is: %.8e\n", densityN);
			 // printf("The value of xwwIdx is: %.8e\n", FluidSystem::molarMass(nCompIdx));
			 // printf("The value of xwwIdx is: %.8e\n", values[conti0EqIdx + wCompIdx]);
			 // printf("The value of xwwIdx is: %.8e\n", gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx)*densityN /FluidSystem::molarMass(nCompIdx));
			 // if (gasFlux < 0)
			 // {
		 	 // values[conti0EqIdx + nCompIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, nCompIdx)*densityN /FluidSystem::molarMass(nCompIdx);
		 	 // values[conti0EqIdx + wCompIdx] = gasFlux * volVars.moleFraction(conti1EqIdx, wCompIdx) *densityN /FluidSystem::molarMass(nCompIdx);
		 	 // values[conti0EqIdx + xwCO2aqIdx] = gasFlux * injCO2aq_ *densityN/FluidSystem::molarMass(nCompIdx);
 		 	 // values[conti0EqIdx + xwO2Idx] = gasFlux * injO2_*densityN /FluidSystem::molarMass(nCompIdx);
			 // }	
			 // printf("The value of conti0EqIdx + wCompIdx2 is: %.10e\n",  values[conti0EqIdx + wCompIdx]);
			 if (waterFlux >= 0)
			 {			 
             values[conti0EqIdx + wCompIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW_ /FluidSystem::molarMass(wCompIdx);		 				
             values[conti0EqIdx + nCompIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, nCompIdx) * densityW_/FluidSystem::molarMass(wCompIdx);
			 values[conti0EqIdx + xwCO2aqIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwCO2aqIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwHIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwHIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwO2Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwO2Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwFe2Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwFe2Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwCaIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwCaIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwKIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwKIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwMgIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwMgIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwAlIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwAlIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwMnIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwMnIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwHPO4Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwHPO4Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwTiOH4Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwTiOH4Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 values[conti0EqIdx + xwSiO2Idx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwSiO2Idx)* densityW_ /FluidSystem::molarMass(wCompIdx); //FluidSystem::molarMass(nCompIdx);
			 // 
			 values[conti0EqIdx + xwCaIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwCaIdx)* densityW_ /FluidSystem::molarMass(CaIdx);
             // // printf("The value of conti0EqIdx + xwCaIdx is: %.2lf\n", xwCaIdx);
			 values[conti0EqIdx + xwNaIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwNaIdx)* densityW_ /FluidSystem::molarMass(wCompIdx); // FluidSystem::molarMass(NaIdx);
             // // printf("The value of values[conti0EqIdx + xwNaIdx] is: %.2lf\n", values[conti0EqIdx + xwNaIdx]);
			 // 
			 values[conti0EqIdx + xwClIdx] += waterFlux * volVars.moleFraction(conti0EqIdx, xwClIdx) * densityW_ /FluidSystem::molarMass(wCompIdx);// FluidSystem::molarMass(NaIdx);               //NaCl ---> mol Cl = mol Na
			 }
			 // printf("The value of conti0EqIdx + wCompIdx3 is: %.10e\n",  values[conti0EqIdx + wCompIdx]);
			 // 
			 // printf("The value of waterFlux is: %.8e\n", waterFlux);
			 // printf("The value of volVars.moleFraction(conti0EqIdx, wCompIdx) is: %.8e\n", volVars.moleFraction(conti0EqIdx, wCompIdx));
			 // printf("The value of densityW is: %.8e\n", densityW);
			 // printf("The value of FluidSystem::molarMass(wCompIdx) is: %.8e\n", FluidSystem::molarMass(wCompIdx));
			 // printf("The value of values[conti0EqIdx + wCompIdx] is: %.8e\n", values[conti0EqIdx + wCompIdx]);
			 // printf("The value of waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW/FluidSystem::molarMass(wCompIdx) is: %.8e\n", waterFlux * volVars.moleFraction(conti0EqIdx, wCompIdx) * densityW/FluidSystem::molarMass(wCompIdx));
// #if NON ISOTHERMAL
 //             values[energyEqIdx] = waterFlux*Brine::liquidEnthalpy(
 //                                      injTemperature_, injPressure_,
 //                                     (injNa_
 //                                     +injNa_ /FluidSystem::molarMass(NaIdx)*FluidSystem::molarMass(ClIdx))/densityW_); // W/(m^2)
 // #endif  
	 	 // //4.640675e-4 m^2 = 0.957" diameter column cross-sectional area. mol/m2 s2
         //     // emulate an outflow condition for the component transport on the right side
         //     // values[contiN2EqIdx] = phaseFlux * (useMoles ? volVars.moleFraction(0, N2Idx) : volVars.massFraction(0, N2Idx));
	     //           
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
		// NumEqVector newdt(0.0);

        Chemistry chemistry;
		// if (episodeIdx_ > 50)
        // {
			
		// source = chemistry.reactionSource(elemVolVars[scv],
                        // timeStepSize_);
		source = std::get<0>(chemistry.reactionSource(elemVolVars[scv],
                        timeStepSize_));
		// newdt = std::get<1>(chemistry.reactionSource(elemVolVars[scv],
                         // timeStepSize_));
		
		// }
		// else				
		// {
		source = 0.0; //added by du
		// }
		
 
        return source;
    }

    // std::vector<Scalar> getnewdt(const GridVariables& gridVariables)
	double getnewdt(const SolutionVector& curSol)
    {
        Chemistry chemistry;
   
        // Iterate over elements in the grid
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
            auto fvGeometry = localView(this->gridGeometry());  // Geometry for finite volume
            fvGeometry.bindElement(element);
			// auto elemVolVars = localView(gridVariables.curGridVolVars()); // Volume variables for the current grid
			
			// auto scvsRange = scvs(fvGeometry);
            // printf("Number of SCVs: %zu\n", scvsRange);
            // Iterate over the sub-control volumes (scvs) within the element
            for (auto&& scv : scvs(fvGeometry))
            {
				VolumeVariables volVars;
                volVars.update(elemSol, *this, element, scv);
				const auto dofIdxGlobal = scv.dofIndex();  // Global degree of freedom index
				newdt_[dofIdxGlobal] = std::get<1>(chemistry.reactionSource(volVars, timeStepSize_));
			}		  
        }
		// printf("The value of newdt_ is: %.12e\n", newdt_);
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
    // const std::vector<Scalar>& getUrea()
    // {
    //     return urea_;
    // }
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
                // urea_[dofIdxGlobal] = volVars.moleFraction(0,UreaIdx)* volVars.molarDensity(0) * FluidSystem::molarMass(UreaIdx);
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


		// priVars.setState(wPhaseOnly); doesn't work, negative pressure
		// priVars.setState(bothPhases); doesn't work, negative pressure	
		priVars.setState(bothPhases);
		
        priVars[pressureIdx] = initPressure_;
		
		// priVars.setState(nPhaseOnly), the following is in water phase 

        // values[0] = 1.0e5*(1.1 - globalPos[dimWorld-1]*0.1);
		// Du added
        // get the water density at atmospheric conditions
		
        // const Scalar densityW = FluidSystem::H2O::liquidDensity(298.15, 1.01325e5);
		
        // assume an initially hydrostatic liquid pressure profile
        // note: we subtract rho_w*g*h because g is defined negative
        // const Scalar pw = 1.0e5 - densityW*this->spatialParams().gravity(globalPos)[dimWorld-1]*(aquiferDepth_ - globalPos[dimWorld-1]);
        // const Scalar pw = 90000 + densityW*this->spatialParams().gravity(globalPos)[dimWorld-1]*(globalPos[dimWorld-1]);
		// const Scalar pw = 89274.3;
		// if (globalPos[1] <= 0.1+eps_)
		// {	
		//    // typename this->spatialParams::MaterialLaw MaterialLaw;
        //    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
        //    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
        //    const auto& materialParams = this->spatialParams().materialLawParamsAtPos(globalPos);
		//    
        //    priVars[switchIdx] = MaterialLaw::sw(materialParams, initPressure_ - pw);
		//    
		//    if (priVars[switchIdx] <= 0.2)
	    //    {
		//    priVars[switchIdx] = 0.2; //Residual water saturation;
		//    }       
		//    if (priVars[switchIdx] >= 1.0)
	    //    {
		//    priVars[switchIdx] = 1.0; //Maximum water saturation;
		//    }		
		// }
		// else
		// {
			priVars[switchIdx] = 0.3;
		// }
		// priVars[switchIdx] = 0.4; //initxwTC_;
		
		// if (globalPos[0] <= eps_)
	    // {
		//	priVars[switchIdx] = 1.0; //initxwTC_;
		// }
        // priVars.setState(wPhaseOnly);
		// initH_ = initxwH_/1.80686E-05; 
		// initCO2aq_ = initxwCO2aq_/1.80686E-05;
		// k1_ = 4.52168E-7;
		// k2_ = 4.69029E-11;
		// kw_ = 1e-14;
        // priVars[xwNaIdx] = initxwNa_ + xwNaCorr_;
        // priVars[xwHIdx] = 1.80686E-05*(initH_-kw_/(initH)-2*k1_*k2_*initCO2aq_ /(initH_)/(initH_)-k1_*initCO2aq_/(initH));
		// priVars[xwCO2aqIdx] = 1.80686E-05*(initCO2aq_+k1_*initCO2aq_/(initH_)+k1_*k2_*initCO2aq_/(initH_)/(initH_)); // give to S_w
        priVars[xwHIdx] = 1.80686E-05*(initxwH_/1.80686E-05-1E-14/(initxwH_/1.80686E-05)-2*4.52168E-7*4.69029E-11*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)/(initxwH_/1.80686E-05)-4.52168E-7*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05));
		priVars[xwCO2aqIdx] = 1.80686E-05*(initxwCO2aq_/1.80686E-05+4.52168E-7*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)+4.52168E-7*4.69029E-11*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)/(initxwH_/1.80686E-05)); // give to S_w

		// priVars[switchIdx+1] = initxwTC_;
        // priVars[xwClIdx] = initxwCl_ + initxwTNH_ + 2*initxwCa_ + xwClCorr_;
        priVars[xwO2Idx] = initxwO2_;
	    priVars[xwClIdx] = initxwCl_ + xwClCorr_;
        // priVars[xwCaIdx] = initxwCa_;
        // priVars[xwTNHIdx] = initxwTNH_;
        priVars[xwFe2Idx] = 1.80686E-05*(initxwFe2_/1.80686E-05+initxwFe2_/1.80686E-05*initxwH_/1.80686E-05*pow((initxwO2_/1.80686E-05),0.25)/1.71633E-08);
        // printf("The value of priVars[xwFe2Idx] is: %.7e\n", priVars[xwFe2Idx]);
        priVars[xwMgIdx] = initxwMg_;
        priVars[xwCaIdx] = initxwCa_;
        priVars[xwKIdx] = initxwK_;
        priVars[xwHPO4Idx] = initxwHPO4_;
		priVars[xwNaIdx] = initxwNa_*1E-20;
        priVars[xwAlIdx] = initxwAl_*1E-20;
        priVars[xwSiO2Idx] = initxwSiO2_*1E-20;
        priVars[xwTiOH4Idx] = initxwTiOH4_*1E-20;
        priVars[xwMnIdx] = initxwMn_*1E-20;
        // priVars[xwHIdx] = initxwH_;
        priVars[phiGlassIdx] = initGlass_; // [m^3/m^3]
        // printf("The value of initxwTiOH4_ is: %.e\n", initxwTiOH4_);
        // priVars[phiCalciteIdx] = initCalcite_; // [m^3/m^3]
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
        // priVars.setState(nPhaseOnly);
        priVars[pressureIdx] = initPressure_;
		// // need to be hydrostatic
        // // priVars[switchIdx] = initxwTC_;
		// priVars[switchIdx] = 0.3;
        priVars.setState(bothPhases);
		// // priVars[xwNaIdx] = initxwNa_ + xwNaCorr_;
		// priVars[xwCO2aqIdx] = initxwCO2aq_; //*2; //initxwCO2_;
        priVars[xwHIdx] = 1.80686E-05*(initxwH_/1.80686E-05-1E-14/(initxwH_/1.80686E-05)-2*4.52168E-7*4.69029E-11*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)/(initxwH_/1.80686E-05)-4.52168E-7*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05));
		priVars[xwCO2aqIdx] = 1.80686E-05*(initxwCO2aq_/1.80686E-05+4.52168E-7*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)+4.52168E-7*4.69029E-11*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)/(initxwH_/1.80686E-05)); // give to S_w
        priVars[xwO2Idx] = initxwO2_;

        // priVars[xwNaIdx] = initxwNa_ + xwNaCorr_;
        // // priVars[xwClIdx] = initxwCl_ + initxwTNH_ + 2*initxwCa_ + xwClCorr_;
	    // priVars[xwClIdx] = initxwCl_ +   xwClCorr_;
        // // priVars[xwHIdx] = initxwH_;
        // // priVars[xwCaIdx] = initxwCa_;
        // // priVars[xwTNHIdx] = initxwTNH_;
        // priVars[xwFe2Idx] = initxwFe2_;
        // priVars[xwMgIdx] = initxwMg_;
        // priVars[xwCaIdx] = initxwCa_;
        // priVars[xwKIdx] = initxwK_;
        // priVars[xwHPO4Idx] = initxwHPO4_;
        // priVars[xwAlIdx] = initxwAl_;
        // priVars[xwSiO2Idx] = initxwSiO2_;
        // priVars[xwTiOH4Idx] = initxwTiOH4_;
        // priVars[xwMnIdx] = initxwMn_;
        // priVars[phiGlassIdx] = initGlass_;
        // priVars[xwHIdx] = initxwH_;
        // priVars[phiGlassIdx] = initGlass_; // [m^3/m^3]
        // priVars[phiCalciteIdx] = initCalcite_; // [m^3/m^3]
        // priVars[phiFerrohydriteIdx] = initFerrohydrite_; // [m^3/m^3]
#if NONISOTHERMAL
        priVars[temperatureIdx] = initTemperature_;
#endif
        return priVars;
    }
    PrimaryVariables bottom_(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        // priVars.setState(nPhaseOnly);
        priVars[pressureIdx] = initPressure_;
		// // need to be hydrostatic
        // // priVars[switchIdx] = initxwTC_;
		// priVars[switchIdx] = 0.3;
        priVars.setState(bothPhases);
		// priVars[xwCO2aqIdx] = initxwCO2aq_;
        priVars[xwO2Idx] = initxwO2_;
        priVars[xwHIdx] = 1.80686E-05*(initxwH_/1.80686E-05-1E-14/(initxwH_/1.80686E-05)-2*4.52168E-7*4.69029E-11*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)/(initxwH_/1.80686E-05)-4.52168E-7*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05));
		priVars[xwCO2aqIdx] = 1.80686E-05*(initxwCO2aq_/1.80686E-05+4.52168E-7*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)+4.52168E-7*4.69029E-11*initxwCO2aq_/1.80686E-05/(initxwH_/1.80686E-05)/(initxwH_/1.80686E-05)); // give to S_w

        // // priVars[xwHIdx] = initxwH_;
        // priVars[xwNaIdx] = initxwNa_ + xwNaCorr_;
        // priVars[xwClIdx] = initxwCl_ + xwClCorr_;;
        // // priVars[xwClIdx] = initxwCl_ + initxwTNH_ + 2*initxwCa_ + xwClCorr_;
        // // priVars[xwCaIdx] = 0.0428;
        // // priVars[xwTNHIdx] = initxwTNH_;
        // priVars[xwFe2Idx] = initxwFe2_;
        // priVars[xwMgIdx] = initxwMg_;
        // priVars[xwCaIdx] = initxwCa_;
        // priVars[xwKIdx] = initxwK_;
        // priVars[xwHPO4Idx] = initxwHPO4_;
        // priVars[phiGlassIdx] = initGlass_;
        // priVars[xwAlIdx] = initxwAl_;
        // priVars[xwSiO2Idx] = initxwSiO2_;
        // priVars[xwTiOH4Idx] = initxwTiOH4_;
        // priVars[xwMnIdx] = initxwMn_;
        // priVars[xwHIdx] = initxwH_;
        // priVars[phiGlassIdx] = initGlass_; // [m^3/m^3]
        // priVars[phiCalciteIdx] = initCalcite_; // [m^3/m^3]
        // priVars[phiFerrohydriteIdx] = initFerrohydrite_; // [m^3/m^3]
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
    // Scalar initSaturation_;
    // Scalar initgasPressure_;
    Scalar densityW_;
    // Scalar densityN_;

    Scalar densityW;
    Scalar densityN;
	
    Scalar initxwCO2aq_;
	Scalar initxwO2_;
    Scalar initxwNa_;
    Scalar initxwH_;
    Scalar initxwCa_;
    Scalar initxwCl_;	
    // Scalar initxwUrea_;
    // Scalar initxwTNH_;
    Scalar initxwFe2_;
    Scalar initxwK_;
	Scalar initxwMg_;
    Scalar initxwAl_;
    Scalar initxwMn_;
    Scalar initxwHPO4_;
    Scalar initxwTiOH4_;
    Scalar initxwSiO2_;	
    // Scalar initxwEnzymeSource_;
    Scalar xwNaCorr_;
    Scalar xwClCorr_;

    Scalar initGlass_;
    // Scalar initCalcite_;
    Scalar initFerrohydrite_;
	Scalar initProtoImogolite_;
	Scalar initBirnessite_;
	Scalar initHydroxyapatite_;
	Scalar initSepiolite_;
	
    // Scalar initImmUrease_;
    // Scalar initImmJBM_;
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
	
    // Scalar injUrea_;
    // Scalar injTNH_;
    Scalar injNaCorr_;
    // Scalar injEnzymeSource_;
    Scalar injTemperature_;
    Scalar injPressure_;

    // Scalar ureaseInJBM_;

    int numInjections_;
    std::string injectionParameters_;

    // std::vector<int> injType_;
    std::vector<double> injType_;	
    std::string name_;

    std::vector<Scalar> permeability_;
    std::vector<Scalar> calcium_;
    // std::vector<Scalar> urea_;
    std::vector<Scalar> fe2_;
    std::vector<Scalar> newdt_;

    Scalar time_ = 0.0;
    Scalar timeStepSize_ = 0.0;
    int episodeIdx_ = 0;
    std::shared_ptr<GridVariables> gridVariables_;
};
} //end namespace

#endif




