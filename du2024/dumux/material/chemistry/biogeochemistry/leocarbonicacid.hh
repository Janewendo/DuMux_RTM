
#ifndef LEO_CARBONIC_ACID_HH_
#define LEO_CARBONIC_ACID_HH_

#include <dumux/common/exceptions.hh>
#include <dumux/material/fluidsystems/leomin.hh>

#include <cmath>
#include <iostream>
#include <dumux/common/math.hh>
#include <stdio.h>

namespace Dumux
{
/*!
 * \brief The equilibrium chemistry is calculated in this class. The function calculateEquilbriumChemistry is used to
 * control the Newton Solver "newton1D". The chemical functions and derivations are implemented in the private part of
 * class.
 */
template <class TypeTag, class CO2Tables, class ModelTraits>
class LeoMinCarbonicAcid
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using SolidSystem = GetPropType<TypeTag, Properties::SolidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Sources = GetPropType<TypeTag, Properties::NumEqVector>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using ThisType = LeoMinCarbonicAcid<TypeTag, CO2Tables, ModelTraits>;
    using H2O = Components::H2O<Scalar>;

    enum
    {
        // phase presence enums
        secondPhaseOnly = ModelTraits::Indices::secondPhaseOnly,
        firstPhaseOnly = ModelTraits::Indices::firstPhaseOnly,
        bothPhases = ModelTraits::Indices::bothPhases,
        nPhaseOnly = secondPhaseOnly,
        wPhaseOnly = firstPhaseOnly,
    };

public:

    LeoMinCarbonicAcid()
{

    pKaFactor_ = getParam<Scalar>("Geochem.pKaFactor");
    glAsw_      = getParam<Scalar>("GlassCoefficients.glAsw");
    glrc_      = getParam<Scalar>("GlassCoefficients.glrc");	
    glbeta_      = getParam<Scalar>("GlassCoefficients.glbeta");
    glsigma_      = getParam<Scalar>("GlassCoefficients.glsigma");
    glp_      = getParam<Scalar>("GlassCoefficients.glp");	

    ferAsw_      = getParam<Scalar>("GlassCoefficients.ferAsw");
    ferrc_      = getParam<Scalar>("GlassCoefficients.ferrc");	
    ferbeta_      = getParam<Scalar>("GlassCoefficients.ferbeta");
    fersigma_      = getParam<Scalar>("GlassCoefficients.fersigma");
    ferp_      = getParam<Scalar>("GlassCoefficients.ferp");

    proAsw_      = getParam<Scalar>("GlassCoefficients.proAsw");
    prorc_      = getParam<Scalar>("GlassCoefficients.prorc");	
    probeta_      = getParam<Scalar>("GlassCoefficients.probeta");
    prosigma_      = getParam<Scalar>("GlassCoefficients.prosigma");
    prop_      = getParam<Scalar>("GlassCoefficients.prop");

    birAsw_      = getParam<Scalar>("GlassCoefficients.birAsw");
    birrc_      = getParam<Scalar>("GlassCoefficients.birrc");	
    birbeta_      = getParam<Scalar>("GlassCoefficients.birbeta");
    birsigma_      = getParam<Scalar>("GlassCoefficients.birsigma");
    birp_      = getParam<Scalar>("GlassCoefficients.birp");

    hydAsw_      = getParam<Scalar>("GlassCoefficients.hydAsw");
    hydrc_      = getParam<Scalar>("GlassCoefficients.hydrc");	
    hydbeta_      = getParam<Scalar>("GlassCoefficients.hydbeta");
    hydsigma_      = getParam<Scalar>("GlassCoefficients.hydsigma");
    hydp_      = getParam<Scalar>("GlassCoefficients.hydp");	
	
    sepAsw_      = getParam<Scalar>("GlassCoefficients.sepAsw");
    seprc_      = getParam<Scalar>("GlassCoefficients.seprc");	
    sepbeta_      = getParam<Scalar>("GlassCoefficients.sepbeta");
    sepsigma_      = getParam<Scalar>("GlassCoefficients.sepsigma");
    sepp_      = getParam<Scalar>("GlassCoefficients.sepp");	

}

    static const int wPhaseIdx    = FluidSystem::wPhaseIdx;
    static const int nPhaseIdx    = FluidSystem::nPhaseIdx;

    static const int wCompIdx     = FluidSystem::wCompIdx;
    static const int nCompIdx     = FluidSystem::nCompIdx;

    static const int H2OIdx       = FluidSystem::H2OIdx;
    static const int N2Idx        = FluidSystem::N2Idx;
    static const int O2Idx        = FluidSystem::O2Idx;

    static const int CaIdx        = FluidSystem::CaIdx;
    static const int NaIdx        = FluidSystem::NaIdx;
    static const int HIdx        = FluidSystem::HIdx;
    static const int HtotalIdx        = FluidSystem::HtotalIdx;

    static const int MgIdx        = FluidSystem::MgIdx;
    static const int AlIdx        = FluidSystem::AlIdx;
    static const int SiO2Idx       = FluidSystem::SiO2Idx;
    static const int HPO4Idx        = FluidSystem::HPO4Idx;
    static const int KIdx        = FluidSystem::KIdx;
    static const int MnIdx       = FluidSystem::MnIdx;
    static const int TiOH4Idx       = FluidSystem::TiOH4Idx;
	
    static const int Fe2Idx       = FluidSystem::Fe2Idx;

    static const int ClIdx         = FluidSystem::ClIdx;
    static const int OHIdx        = FluidSystem::OHIdx;

    static const int CO2aqIdx       = FluidSystem::CO2aqIdx;
    static const int CO2aqtotalIdx       = FluidSystem::CO2aqtotalIdx;
    static const int HCO3Idx      = FluidSystem::HCO3Idx;
    static const int CO3Idx       = FluidSystem::CO3Idx;
    static const int Fe2totalIdx       = FluidSystem::Fe2totalIdx;
    static const int Fe3Idx       = FluidSystem::Fe3Idx;
	

    static const int numComponents      = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numMajorComponents;
    static const int numSecComponents   = FluidSystem::numSecComponents;
    static const int numTotComponents   = numComponents + numSecComponents;
    static const int numPhases          = FluidSystem::numPhases;

    static const int fPhaseIdx          = SolidSystem::FerrohydriteIdx;	
	static const int gPhaseIdx          = SolidSystem::GlassIdx;
	static const int pPhaseIdx          = SolidSystem::ProtoImogoliteIdx;
	static const int hPhaseIdx          = SolidSystem::HydroxyapatiteIdx;
	static const int bPhaseIdx          = SolidSystem::BirnessiteIdx;
	static const int sPhaseIdx          = SolidSystem::SepioliteIdx;
    static const int numSolidComponents = SolidSystem::numComponents;
    static const int numInertComponents = SolidSystem::numInertComponents;

    static const int phiGlassIdx      = numComponents + gPhaseIdx;
    static const int phiFerrohydriteIdx    = numComponents + fPhaseIdx;
	static const int phiProtoImogoliteIdx    = numComponents + pPhaseIdx;
	static const int phiBirnessiteIdx    = numComponents + bPhaseIdx;
	static const int phiHydroxyapatiteIdx    = numComponents + hPhaseIdx;
	static const int phiSepioliteIdx    = numComponents + sPhaseIdx;

    typedef Dune::FieldVector<Scalar, 4> Vector;   // Ionic Strength with NH4/totalnh
    typedef Dune::FieldVector<Scalar, 2> SolVector;
    typedef Dune::FieldVector<Scalar, numTotComponents> CompVector;

    typedef CompositionalSecCompFluidState<Scalar, FluidSystem> FluidState;

    template <class FluidState>
    void calculateEquilibriumChemistry(const FluidState &fluidState, int phaseState, CompVector &variable, Scalar rhoMolar)
 {
        const VolumeVariables volVars{};
        gammaCO2_ = 1.0;
        // h2o_ = 55.508; //molH2O/kgH2O
        Scalar h2o_ = moleFracToMolarity(variable[H2OIdx], rhoMolar);
        pressure_ = fluidState.pressure(wPhaseIdx);
        temperature_ = fluidState.temperature();
		Scalar moleFracWater = variable[H2OIdx]; //
		Scalar n2_ = moleFracToMolarity(variable[N2Idx], rhoMolar);		
        for (int i = 0; i < numComponents + numSecComponents; ++i)
            {				
                if(std::isnan(variable[i]))
                {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid component mole fraction " << variable); break;
                }
            }
        if (phaseState == bothPhases) //both Phases: solve an open system with co2 concentration constant
        {
			co2aq_ = moleFracToMolarity(variable[CO2aqIdx],rhoMolar);//-moleFracToMolarity(variable[CO3Idx],rhoMolar)-moleFracToMolarity(variable[HCO3Idx],rhoMolar);		
            o2_ = moleFracToMolarity(variable[O2Idx],rhoMolar);		
			ca_ = moleFracToMolarity(variable[CaIdx],rhoMolar);
            cl_ = moleFracToMolarity(variable[ClIdx], rhoMolar);//moleFracWater);
            fe2_ = moleFracToMolarity(variable[Fe2Idx], rhoMolar);
			na_ = moleFracToMolarity(variable[NaIdx],rhoMolar);
            h_ = moleFracToMolarity(variable[HIdx],rhoMolar);//+moleFracToMolarity(variable[OHIdx],rhoMolar)+moleFracToMolarity(variable[CO3Idx],rhoMolar)+2*moleFracToMolarity(variable[HCO3Idx],rhoMolar);
            mg_ = moleFracToMolarity(variable[MgIdx],rhoMolar);
            al_ = moleFracToMolarity(variable[AlIdx],rhoMolar);
            sio2_ = moleFracToMolarity(variable[SiO2Idx],rhoMolar);
            hpo4_ = moleFracToMolarity(variable[HPO4Idx],rhoMolar);
            k_ = moleFracToMolarity(variable[KIdx],rhoMolar);
            mn_ = moleFracToMolarity(variable[MnIdx],rhoMolar);
            tioh4_ = moleFracToMolarity(variable[TiOH4Idx],rhoMolar);

		   initCO2_ = co2aq_;//h_; //Initial guess
           Scalar activityCO2 = initCO2_;		
	       k1_ = const1(pressure_, temperature_);
           k2_ = const2(pressure_, temperature_);
           k3_ = const3(pressure_, temperature_);
           k4_ = const4(pressure_, temperature_);
           k5_ = const5(pressure_, temperature_);  
	       k6_ = const6(pressure_, temperature_);
           k7_ = const7(pressure_, temperature_);
           k8_ = const8(pressure_, temperature_);
           k9_ = const9(pressure_, temperature_);
           k10_ = const10(pressure_, temperature_); 
		   k11_ = const11(pressure_, temperature_);
           kw_ = constW(pressure_, temperature_);
		   
           Scalar tolAbs = 1e-20; //1e-20;
           Scalar tolRel = 1e-20;// 1e-15;
           int maxIter = 40;
            // printf("The value of CO2_ is: %.10e\n", co2aq_);		   
			// if(newton1D(activityCO2, &ThisType::H_Conly, tolAbs, tolRel, maxIter) == false) //Alex' Newton
            // {
            //     initCO2_ = co2aq_; //8.07096e-13 ;//h_;
			// 	activityCO2 = initCO2_;		
            //     Scalar a0 = 0.0;
            //     Scalar b0 = 1e-1;
            //     Scalar tol = 1e-15;//1e-15;
			// 	if(bisection1D(activityCO2, &ThisType::H_Conly, a0, b0, tol) == false) //Alex' bisection
			// 
            //     {
            //         DUNE_THROW(Dune::InvalidStateException, "in Chemistry: Bisection did not converge!" );
            //     }
            // }

		    // H_Conly(activityCO2);
			// printf("The value of CO22_ is: %.10e\n", co2aq_);
            // printf("The value of variable[HIdx] is: %.10e\n", variable[HIdx]);
			// printf("The value of H_ is: %.10e\n", h_);
            oh_ = kw_/h_;
		    hco3_ = k1_*co2aq_/(h_);
            co3_ = k1_*k2_*co2aq_/(h_*h_);
		    fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(kw_/h_, 4)), 0.25);
		
            //update mole fractions in the variable vector for the open system
		    htotal_ = h_ - oh_ - 2*co3_ - hco3_ + 2*fe3_ ;
		    co2aqtotal_ = co2aq_ + co3_ + hco3_;
			fe2total_ = fe2_ + fe3_ ; 		
            Scalar totalMolarity = h2o_ + n2_ + o2_ +  co2aqtotal_  + na_ + cl_ + ca_ + fe2total_ + al_ + mg_ +k_ + mn_ + hpo4_ + tioh4_ + sio2_ + htotal_; // du do not know
            // printf("The value of totalMolarity is: %.10e\n", totalMolarity);			
			variable[CO2aqIdx] = co2aq_/totalMolarity;
			variable[CO2aqtotalIdx] = co2aqtotal_/totalMolarity;
            variable[HCO3Idx] = hco3_/totalMolarity;
            variable[CO3Idx] = co3_/totalMolarity;
            variable[OHIdx] = oh_/totalMolarity;
            variable[HIdx] = h_/totalMolarity;
            variable[HtotalIdx] = htotal_/totalMolarity;
			variable[Fe2totalIdx] = fe2total_/totalMolarity;
			variable[Fe2Idx] = fe2_/totalMolarity;
			variable[Fe3Idx] = fe3_/totalMolarity;
			// printf("The value of H_2 is: %.10e\n", h_);
            // printf("The value of variable[HIdx]2 is: %.10e\n", variable[HIdx]);	
            // printf("The value of variable[Fe2totalIdx]2 is: %.10e\n", variable[Fe2totalIdx]);
            // printf("The value of variable[Fe2Idx]2 is: %.10e\n", variable[Fe2Idx]);	
            // printf("The value of variable[Fe3Idx]2 is: %.10e\n", variable[Fe3Idx]);					
        }

        else if (phaseState == wPhaseOnly) //wPhaseOnly: solve a closed system with cTot concentration constant
        {
			co2aq_ = moleFracToMolarity(variable[CO2aqIdx],rhoMolar);				
			ca_ = moleFracToMolarity(variable[CaIdx],rhoMolar);
			na_ = moleFracToMolarity(variable[NaIdx],rhoMolar);
            cl_ = moleFracToMolarity(variable[ClIdx],rhoMolar);
            fe2_ = moleFracToMolarity(variable[Fe2Idx],rhoMolar);			
            kw_ = constW(pressure_, temperature_);


            Scalar tolAbs = 1e-20; //1e-11;old
            Scalar tolRel = 1e-20; //1e-11;old
            int maxIter = 40;
            initH_ = 1e-7;
            Scalar activityH = initH_;

           for (int i = 0; i < numComponents + numSecComponents; ++i)
           {
              if(std::isnan(variable[i]))
              {
                  std::cout<<"wPhaseOnly, moleFrac of  "<< FluidSystem::componentName(i) << "  in chemistry is: "<< variable[i]<< " all moleFracs are: \n"  << variable<<"\n" <<std::endl;
              }
           }

        }

        else if (phaseState == nPhaseOnly) //no secondary components in the gas phase, except CTot = CO2!
        {
            // variable[CO2aqIdx] = variable[CO2aqIdx]; 
			// // variable[CO2Idx] = variable[CO2Idx];
			// // variable[CO2Idx] = variable[CO2Idx];//all cTot in the gas is CO2
            // // variable[HCO3Idx] = 0.0;
            // // variable[CO3Idx] = 0.0;
            // // variable[NH4Idx] = 0.0;
            // variable[OHIdx] = 0.0;
            // variable[HIdx] = 0.0;
        }

        else
        {
            DUNE_THROW(Dune::InvalidStateException, "Invalid phaseState" );
        }
		// return variable;    // added by du  
    }

    //Return equlibrium constant for chemical equation:
    //H2CO3 <--> H + HCO3
   // static Scalar const1(const Scalar pw, const Scalar T)
   // {
   //     return 4.52168e-7;//return(pow(10,-6.3447));
   // }
 
    // H2O +CO2  <--> H + HCO3
    static Scalar const1(const Scalar pw, const Scalar T)
    {
        return 4.52168e-7;//return(pow(10,-6.3447));
    }
   
    //Return equlibrium constant for chemical equation:
    //HCO3 <--> H + CO3
    static Scalar const2(const Scalar pw, const Scalar T)
    {
        return 4.69029e-11;//return(pow(10,-10.3288));
    }
	
	//Return equlibrium constant for chemical equation:
    //4Fe2+ + O2 + 2H2O = 4Fe3+ + 4OH-
    static Scalar const3(const Scalar pw, const Scalar T)
    {
        return 5.01187e-11;//return(pow(10,-7.7654));
    }
	//Return equlibrium constant for chemical equation:
    //'Glass_FB' 92.6645  12  1.80 'SiO2(aq)' 2.05 'H2O' 0.07 
	// 'Ti(OH)4(aq)' 0.62 'Al+++' 0.38 'Fe++' 0.01 'Mn++' 0.30 
	// 'Mg++' 0.44 'Ca++' 0.26  'Na+' 0.06 'K+' 0.03 'HPO4--' 
	// -4.38 'H+'  17.61  
    static Scalar const4(const Scalar pw, const Scalar T)
    {
        return 4.0738E+17;//return(pow(10,17.61));
    }	
	//return(pow(10,-7.7654));
    // //Return equlibrium constant for chemical equation:
    // //CO2(g) <--> CO2(aq)
    // static Scalar constco2(const Scalar pw, const Scalar T)
    // {
    //     return 3.40408e-2;//return(pow(10,-1.468));
    // }
    // Return equlibrium constant for dissolution reaction:
    //  CaCO3(s) +2H <--> Ca + CO2 + H2O
	// CO2(g)' 0.0000 3 -1.0000 'H2O' 1.0000 'H+' 1.0000 'HCO3-' -7.8136 
	// 'CO2(g)*' 0.0000 1 1.0000 'CO2(aq)' -1.468
	// 'Calcite' 36.9340 3 -1.0000 'H+' 1.0000 'Ca++' 1.0000 'HCO3-' 1.71
    static Scalar const5(const Scalar pw, const Scalar T)
    {
        return 2.3055E+08;//return(pow(10,8.3628)); //may not right
    }
    // ferrohydrite
    static Scalar const6(const Scalar pw, const Scalar T)
    {
        return 7.7625E+04;//return(pow(10,4.89)); //may not right
    }
// proto
    static Scalar const7(const Scalar pw, const Scalar T)
    {
        return 1.0471E+07;//return(pow(10,7.02)); //may not right
    }
// birnessite
    static Scalar const8(const Scalar pw, const Scalar T)
    {
        return 4.7863E+11;//return(pow(10,11.68)); //may not right
    }
// Hydroxyapatite
    static Scalar const9(const Scalar pw, const Scalar T)
    {
        return 8.5114E-04;//return(pow(10,-3.07)); //may not right
    }
// sepiolite
    static Scalar const10(const Scalar pw, const Scalar T)
    {
        return 2.7542E+30;//return(pow(10,30.44)); //may not right
    }	
// Fe
// 'Fe+++' 4 -0.5000 'H2O' 0.2500 'O2(g)' 1.0000 'Fe++' 1.0000 'H+' -7.7654
    static Scalar const11(const Scalar pw, const Scalar T)
    {
        return 1.71633E-08;//return(pow(10,-7.7654)); //may not right
    }	
//  //    CaCO3(s) +H <--> Ca + HCO3
    // static Scalar const5(const Scalar pw, const Scalar T)
    // {
    //     return 51.2861384;//return(pow(10,1.71));
    // }	
// //    Return equlibrium constant for dissolution reaction:
// //    CaCO3(s) <--> Ca + CO3
    // static Scalar solubilityProductCaCO(const Scalar pw, const Scalar T)
    // {
    //     return 4.8e-9;
    // }
// 
// //    Return equlibrium constant for dissolution reaction:
// //    Fe(OH)2(s) <--> Fe + 2OH
//     static Scalar solubilityProductFeOH2(const Scalar pw, const Scalar T)
//     {
//         return(pow(10,-4.89));
// 		 // added by DU
//     }

    //Return equlibrium constant for chemical equation:
    // H2O <--> H + OH
    static Scalar constW(const Scalar pw, const Scalar T)
    {
        return 1e-14;
    }


    static Scalar massFracToMolality(const Scalar massFracX, const Scalar molarMassX, const Scalar massFracSalinity,
            const Scalar massFracC)
    {
        Scalar molalityX = massFracX/molarMassX/(1- massFracSalinity - massFracC);
        return molalityX;
    }


    /*!
     * \brief Returns the mass fraction of a component x (kg x / kg solution) for a given
     * molality fraction (mol x / mol solution)
     * The salinity and the mole Fraction of CO2 are considered
     *
     */

    static Scalar molalityToMassFrac(Scalar molalityX, Scalar molarMassX, Scalar massFracSalinity, Scalar massFracCTot)
    {
        Scalar massFracX = molalityX * molarMassX * (1 - massFracSalinity - massFracCTot);
        return massFracX;
    }

    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        if(moleFracX<0)
            moleFracX=0;
        if(moleFracSalinity<0)
            moleFracSalinity=0;
        if(moleFracCTot<0)
            moleFracCTot=0;

        Scalar molalityX = moleFracX / (1 - moleFracSalinity - moleFracCTot) / FluidSystem::molarMass(H2OIdx);
        return molalityX;
    }

    static Scalar molalityToMoleFrac(Scalar molalityX, Scalar moleFracSalinity, Scalar moleFracCTot)
    {
        Scalar moleFracX = molalityX * (1 - moleFracSalinity - moleFracCTot) * FluidSystem::molarMass(H2OIdx);
        return moleFracX;
    }

    static Scalar moleFracToMolality(Scalar moleFracX, Scalar moleFracWater)
    {
        if(moleFracX<0)
            moleFracX=0;
        if(moleFracWater<0)
            moleFracWater=0;
        Scalar molalityX = moleFracX / moleFracWater / FluidSystem::molarMass(H2OIdx);
        return molalityX;
    }

	
    static Scalar molalityToMoleFrac(Scalar molalityX, Scalar moleFracWater)
    {
        Scalar moleFracX = molalityX * moleFracWater * FluidSystem::molarMass(H2OIdx);
        return moleFracX;
    }

    // added by du
    static Scalar moleFracToMolarity(Scalar moleFracX, Scalar molarDensity)
    {
        // if(moleFracX<0)
        //     moleFracX=0;
        Scalar molarityX = moleFracX * molarDensity; 
        return molarityX;
    }

    static Scalar molarityToMoleFrac(Scalar molarityX, Scalar molarDensity)
    {
        Scalar moleFracX = molarityX / molarDensity; 
        return moleFracX;
    }
    /*!
     * \brief The ionic strength of a substance is calculated only with the salinity until now!
     */

    static Scalar ionicStrength(Scalar molaritySalinity)
    {
        Scalar ionicStrength = 0.0;
        //Scalar Product
        for (int compIdx = 0; compIdx < 2; ++compIdx)
        {
            ionicStrength += molaritySalinity;
        }
        ionicStrength *= 0.5;

        return ionicStrength;
    }
    // static Scalar ionicStrength(Scalar mNa, Scalar mCl, Scalar mCa, Scalar mNH4, Scalar mFe2 )
    static Scalar ionicStrength(Scalar mNa, Scalar mCl, Scalar mCa, Scalar mFe2 )
     {
        Scalar ionicStrength = 0.5*( mNa    * FluidSystem::charge(NaIdx) * FluidSystem::charge(NaIdx)
        + mCl   * FluidSystem::charge(ClIdx) * FluidSystem::charge(ClIdx)
        + mCa   * FluidSystem::charge(CaIdx) * FluidSystem::charge(CaIdx)
		+ mFe2  * FluidSystem::charge(Fe2Idx) * FluidSystem::charge(Fe2Idx));
        return ionicStrength;
    }

    static Scalar ionicStrength(Scalar mNa, Scalar mCl, Scalar mCa )
    {
        Scalar ionicStrength = 0.5*( mNa    * FluidSystem::charge(NaIdx) * FluidSystem::charge(NaIdx)
        + mCl   * FluidSystem::charge(ClIdx) * FluidSystem::charge(ClIdx)
        + mCa   * FluidSystem::charge(CaIdx) * FluidSystem::charge(CaIdx));

        return ionicStrength;
    }

    void ionicStrength()
    {
        ionicStrength_ = 0.0;
        //Scalar Product
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            ionicStrength_ += molarity_[compIdx] * charge_[compIdx] * charge_[compIdx];
        }
        ionicStrength_ *= 0.5;
    }

    //Calculates the activity with a modified Debye-Hückel equation after Parkhurst (1990) for
    //ionic strengths up to 2.
    static Scalar activityCoefficient(Scalar ionicStrength, Scalar temperatureK, int compIdx)
    {
        if (ionicStrength<0)
        {
            ionicStrength = 0;
        }
        Scalar charge = FluidSystem::charge(compIdx);
        Scalar ai = FluidSystem::ai(compIdx);
        Scalar bi = FluidSystem::bi(compIdx);
        Scalar A = 0.5085;
        Scalar B = 0.3285e10;
        // The actual modified Debye Hückel equation
        Scalar logActivityCoefficient = -A*(charge*charge)*sqrt(ionicStrength)/(1 + B*ai*sqrt(ionicStrength))
                + bi*ionicStrength;

        return pow(10, logActivityCoefficient);
//          return 1.0;
    }

    static Scalar J(Scalar x)
    {

      Scalar c[5], res;

      /*Pitzer 1974, Thermodaynamics of Electrolytes V*/

      c[1]=4.581;  c[2]=0.7237;  c[3]=0.0120;  c[4]=0.528;

      res = x/(4. + c[1]*pow(x,-c[2])*exp(-c[3]*pow(x,c[4])));


      return(res);
    }
    static Scalar Jprime(Scalar x)
    {

      Scalar res, eps;

      eps = 1.E-3;

      res = (J(x+eps) - J(x))/eps;


      return(res);
    }

    // static Scalar Appa_Ksp(Scalar mNa, Scalar mCa, Scalar mNH4, Scalar mHCO3, Scalar mCO3, Scalar mCl, Scalar temp)
    static Scalar Appa_Ksp(Scalar mNa, Scalar mCa, Scalar mCl, Scalar temp)
    {

      Scalar f, B_cacl, C_cacl, B1_cacl, I, sqrt_I, gamma_Ca, gamma_CO3, Ksp;
      Scalar beta_cacl_0, beta_cacl_1, C_cacl_phi;
      Scalar beta_nacl_0, beta_nacl_1, C_nacl_phi;
      Scalar psi_canacl, theta_naca;
      Scalar B1_nacl, C_nacl;
      Scalar A_phi, a[6], T,x_clcl, x_cana,x_caca,x_nana;
      Scalar E_theta_cana,  E1_theta_cana;
      I = 0.5*( mNa + 4.*mCa + mCl) + 1.E-20;
      sqrt_I = sqrt(I);

      T = temp;
      a[0]=-8.1765300E-1; a[1]=-8.6852760E-1; a[2]=1.9251000E+4; a[3]=5.2514840E-3; a[4]=-7.1493970E-6; a[5]=9.3385590E-12;

      A_phi = a[0] + a[1]/(T-222.) + a[2]/(T*T) + a[3]*T + a[4]*T*T + a[5]*T*T*T*T;
      /*MODELING AND NUMERICAL SIMULATION OF SALT TRANSPORT AND PHASE TRANSITIONS IN UNSATURATED POROUS BUILDING MATERIALS By Andreas Nicolai*/

      beta_cacl_0 = 0.3159;  beta_cacl_1 = 1.614; C_cacl_phi = -0.00034;
      beta_nacl_0 = 0.0765; beta_nacl_1 = 0.2664; C_nacl_phi = 0.00127;
      psi_canacl = -0.014; // psi_co3nacl = 0.016;
      theta_naca = 0.07; // theta_clco3 = -0.053;



      x_clcl = 6.*(-1.)*(-1.)*A_phi*sqrt_I;
      x_cana = 6.*(+2.)*(+1.)*A_phi*sqrt_I;
      x_caca = 6.*(+2.)*(+2.)*A_phi*sqrt_I;
      x_nana = 6.*(+1.)*(+1.)*A_phi*sqrt_I;


      E_theta_cana = ((+2.)*(+1.)/(4.*I))*( J(x_cana) - 0.5*J(x_caca) - 0.5*J(x_nana) );

      E1_theta_cana = -(E_theta_cana/I) + ((+2)*(+1)/(8*I*I))*( x_cana*Jprime(x_cana) - 0.5*x_caca*Jprime(x_caca) - 0.5*x_nana*Jprime(x_nana) );

        f = -A_phi * ( sqrt_I/(1. + 1.2*sqrt_I) + (2./1.2)*log(1. + 1.2*sqrt_I) );
        B_cacl = beta_cacl_0 + (beta_cacl_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));
        B1_cacl = (beta_cacl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
        C_cacl = C_cacl_phi  / (2.*sqrt(2.*1.));

        B1_nacl = (beta_nacl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
        C_nacl = C_nacl_phi  / (2.*sqrt(1.*1.));

        gamma_Ca = exp (
                4.*f
                + mCl*(2.*B_cacl + mCl*C_cacl)
                + mNa*mCl*(4.*B1_nacl + 2.*C_nacl)
                + mCa*mCl*(4.*B1_cacl + 2.*C_cacl)



                + mNa*(2.*theta_naca + 2.*E_theta_cana + mCl*psi_canacl)
                + 4.*mNa*mCa*E1_theta_cana // + 4.*mCl*mCO3*E1_theta_clco3
                );



        Ksp = 3.31131e-9/(gamma_Ca);
		 // Activity coefficients are calculated using an exponential function with 
		 // various parameters that account for interactions between ions in the solution.
      return(Ksp);
    }

   Scalar pH(const VolumeVariables &volVars)
   {
      //Scalar mH = moleFracToMolality(volVars.moleFraction(wPhaseIdx,HIdx), volVars.moleFracSalinity(), volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_H/kg_H2O]
      Scalar mH = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,HIdx));  //[mol_H/kg_H2O]

      Scalar pH = -log10(mH);
         return pH;
	  // printf("The value of pH is: %.e\n", pH); 
   }


   std::pair<Sources, Scalar> reactionSource(const VolumeVariables &volVars, const Scalar dt)
   // Sources reactionSource(const VolumeVariables &volVars,
            // const Scalar dt)
    {
        Sources q(0.0);
		Scalar newdt;
    //  //define and compute some parameters for simplicity:
        Scalar temperature = volVars.temperature(); //temperature
        Scalar temperatureC = volVars.temperature() - 273.15; //temperature in °C
        Scalar porosity = volVars.porosity();

        
		const Scalar initialvolFracGlass = 6.1e-1; // volVars.solidVolumeFraction(gPhaseIdx);
        const Scalar initialvolFracFerrohydrite = 1e-30;
		const Scalar initialvolFracProtoImogolite = 1e-30;
		const Scalar initialvolFracBirnessite = 1e-30;
		const Scalar initialvolFracHydroxyapatite = 1e-30;
		const Scalar initialvolFracSepiolite = 1e-30;

        const Scalar initialPorosity = 1.0-initialvolFracGlass-initialvolFracFerrohydrite-initialvolFracProtoImogolite- initialvolFracBirnessite-initialvolFracHydroxyapatite-initialvolFracSepiolite; // tempPorosity;
		   
        Scalar Sw   =  volVars.saturation(wPhaseIdx);
        Scalar xWWater = volVars.moleFraction(wPhaseIdx,wCompIdx);
		
        Scalar volFracGlass = volVars.solidVolumeFraction(gPhaseIdx);
        if (volFracGlass < 1e-30) // Added by du
        {
		volFracGlass = 1e-30;}
        
        Scalar volFracFerrohydrite = volVars.solidVolumeFraction(fPhaseIdx);
        if (volFracFerrohydrite < 1e-30) // Added by du
        { 
		volFracFerrohydrite = 1e-30;}
		
        Scalar volFracProtoImogolite = volVars.solidVolumeFraction(pPhaseIdx);
        if (volFracProtoImogolite < 1e-30) // Added by du
        {
		volFracProtoImogolite = 1e-30;}		
		
        Scalar volFracBirnessite = volVars.solidVolumeFraction(bPhaseIdx);
		if (volFracBirnessite < 1e-30) // Added by du
		{
        volFracBirnessite = 1e-30;
	    }
		
		Scalar volFracHydroxyapatite = volVars.solidVolumeFraction(hPhaseIdx);
        if (volFracHydroxyapatite < 1e-30) // Added by du
        {
		volFracHydroxyapatite = 1e-30;}
        
		Scalar volFracSepiolite = volVars.solidVolumeFraction(sPhaseIdx);
        if (volFracSepiolite < 1e-30) // Added by du
        {
		volFracSepiolite = 1e-30;}		

		Scalar moleFracWater = volVars.moleFraction(wPhaseIdx,H2OIdx);

		// htotal_ = h_ - oh_ - 2*co3_ - hco3_ + fe3_ ;
		// co2aqtotal_ = co2aq_ + co3_ + hco3_;
		// fe2total_ = fe2_ + fe3_ ; 
			
        Scalar mH = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,HIdx),volVars.molarDensity(wPhaseIdx));  //[mol_H/kg_H2O]
        Scalar mHtotal = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,HtotalIdx),volVars.molarDensity(wPhaseIdx));  //[mol_H/kg_H2O]
        Scalar mCO2aq = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,CO2aqIdx), volVars.molarDensity(wPhaseIdx));  //[mol_CO3/kg_H2O]
        Scalar mCO2aqtotal = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,CO2aqtotalIdx), volVars.molarDensity(wPhaseIdx));  //[mol_CO3/kg_H2O]
        Scalar mCO3 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,CO3Idx), volVars.molarDensity(wPhaseIdx));  //[mol_CO3/kg_H2O]
        Scalar mHCO3= moleFracToMolarity(volVars.moleFraction(wPhaseIdx,HCO3Idx), volVars.molarDensity(wPhaseIdx));  //[mol_CO3/kg_H2O]
        Scalar mFe2total = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,Fe2totalIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        Scalar mFe2 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,Fe2Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        Scalar mFe3 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,Fe3Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        Scalar mOH = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,OHIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        // printf("The value of mH is: %.7e\n", mH);
		// printf("The value of mCO2aq is: %.7e\n", mCO2aq);
		// printf("The value of mFe2 is: %.7e\n", mFe2);
		mH = mHtotal + mOH + 2* mCO3 + mHCO3 - 2*mFe3;
		mCO2aq = mCO2aqtotal - mCO3 - mHCO3;	
		mFe2 = mFe2total - mFe3;
		// printf("The value of mH2 is: %.7e\n", mH);
		// printf("The value of mCO2aq2 is: %.7e\n", mCO2aq);
		// printf("The value of mFe22 is: %.7e\n", mFe2);
        if (mOH < 0)
		{
			mOH = 0;
		}       
		if (mFe2 < 0 )//|| std::isnan(mFe2)) // added by DU
		{
			mFe2 = 0;
		}
        if (mFe3 < 0 )//|| std::isnan(mFe2)) // added by DU
		{
			mFe3 = 0;
		}
        if (mCO2aq < 0)
        {
			mCO2aq = 0;
		} 
        if (mH < 0 )//|| std::isnan(mH))
		{
			mH = 0;
		}
		
		Scalar mNa = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,NaIdx),volVars.molarDensity(wPhaseIdx));  //[mol_sodium/kg_H2O]
        if (mNa < 0)// || std::isnan(mNa))
		{
			mNa = 0;
		}

        Scalar mCa = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,CaIdx),volVars.molarDensity(wPhaseIdx));  //[mol_calcium/kg_H2O]
        if (mCa < 0 )//|| std::isnan(mCa))
		{
			mCa = 0;
			
		}

        Scalar mO2 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,O2Idx), volVars.molarDensity(wPhaseIdx));  //[mol_CO3/kg_H2O]
        if (mO2 < 0)
        {
			mO2 = 0;
		}


        Scalar mSiO2 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,SiO2Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mSiO2 < 0 )//|| std::isnan(mSiO2)) // added by DU
		{
			mSiO2 = 0;
		}
        Scalar mTiOH4 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,TiOH4Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mTiOH4 < 0 )//|| std::isnan(mTiOH4)) // added by DU
		{
			mTiOH4 = 0;
		}
        Scalar mAl = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,AlIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mAl < 0 )//|| std::isnan(mAl)) // added by DU
		{
			mAl = 0;
		}
        Scalar mMn = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,MnIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mMn < 0 )//|| std::isnan(mMn)) // added by DU
		{
			mMn = 0;
		}
        Scalar mMg = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,MgIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mMg < 0 )// || std::isnan(mMg)) // added by DU
		{
			mMg = 0;
		}
        Scalar mK = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,KIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mK < 0 )//|| std::isnan(mK)) // added by DU
		{
			mK = 0;
		}
        Scalar mHPO4 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,HPO4Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mHPO4 < 0)// || std::isnan(mHPO4)) // added by DU
		{
			mHPO4 = 0; 
		}

        // compute dissolution and precipitation rate of basalt glass
        Scalar glsp = const4(pressure_, temperature_);
        Scalar glOmegaApprox_ = pow(mSiO2,1.8) * pow(mTiOH4,0.07) * pow(mAl,0.62) * pow(mFe2,0.38) * pow(mMn,0.01)* pow(mMg,0.3)
		* pow(mCa,0.44)* pow(mNa,0.26)* pow(mK,0.06)* pow(mHPO4,0.03)/ pow (mH,4.38) /glsp;
		Scalar glAw0 = glAsw_ * initialvolFracGlass * volVars.solidComponentDensity(gPhaseIdx);
 		Scalar glAwcd = glAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracGlass/initialvolFracGlass)*(volFracGlass/initialvolFracGlass));   // TODO Asw should be a function of Sw, too!
 		Scalar glAwcp = glAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!

        if (glAwcd < 0 || std::isnan(glAwcd))
        {
        std::cout<< "glAwcd = "<<glAwcd<<std::endl;
        glAwcd = 0;
        std::cout<< "glAwcd, corrected = "<<glAwcd<<std::endl;
        }

        if (glAwcp < 0 || std::isnan(glAwcp))
        {
        std::cout<< "glAwcp = "<<glAwcp<<std::endl;
        glAwcp = 0;
        std::cout<< "glAwcp, corrected = "<<glAwcp<<std::endl;
        }
		
		Scalar glr = 0;
			if ( glOmegaApprox_ < 1 && glr > 1e-6)
			{
			}		

        Scalar glrdiss = 0;
        Scalar glrprec = 0;
		
        if (glOmegaApprox_ > 1)
        {
		    glr = glAwcp * pow(10,glrc_) * glp_ * pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_);  
            glrdiss = 0;
            glrprec = glr;//[mol/m³s] // [mol/kg s]

	
        }
        else
        {
		    glr = glAwcd * pow(10,glrc_) * glp_ * pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_); 
			glrdiss = glr; //[mol/kg s]
            glrprec = 0;
			if (std::isnan(glr))
			{

			}
        }

		
        // ferrohydrite
        Scalar fersp = const6(pressure_, temperature_);
        Scalar ferOmegaApprox_ = pow(mFe3,1) / pow (mH,3) /fersp;
        Scalar ferAw0 = ferAsw_ * initialvolFracFerrohydrite * volVars.solidComponentDensity(fPhaseIdx);
 		Scalar ferAwcp = ferAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!
 		Scalar ferAwcd = ferAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracFerrohydrite*volFracFerrohydrite));   // TODO Asw should be a function of Sw, too!
		Scalar ferAwf = ferAsw_ * volFracFerrohydrite * volVars.solidComponentDensity(fPhaseIdx);
		
		Scalar ferAwp = ferAwcp;
		Scalar ferAwd = ferAwcd;
				

        if (ferAwd < 0 || std::isnan(ferAwd))
        {
        std::cout<< "ferAwd = "<<ferAwd<<std::endl;
        ferAwd = 0;
        std::cout<< "ferAwd, corrected = "<<ferAwd<<std::endl;
        }
		
        if (ferAwp < 0 || std::isnan(ferAwp))
        {
        std::cout<< "ferAwp = "<<ferAwp<<std::endl;
        ferAwp = 0;
        std::cout<< "ferAwp, corrected = "<<ferAwp<<std::endl;
        }		

        Scalar ferrdiss = 0;
        Scalar ferrprec = 0;
		
        if (ferOmegaApprox_ > 1)
        {
		    Scalar ferr = ferAwp * pow(10,ferrc_) * ferp_ * pow(abs(1-pow(ferOmegaApprox_,1/fersigma_)),ferbeta_);
            ferrdiss = 0;
            ferrprec = ferr;//[mol/kg s]
        }
        else
        {
		    Scalar ferr = ferAwd * pow(10,ferrc_) * ferp_ * pow(abs(1-pow(ferOmegaApprox_,1/fersigma_)),ferbeta_);
			ferrdiss = ferr; //[mol/kg s]
            ferrprec = 0;
			if (std::isnan(ferr))
			{
			}
        }
		
		
		// if (ferOmegaApprox_ > 1)
		// {
		// if(ferrprec > volVars.moleFraction(wPhaseIdx,Fe3Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     ferrprec =  volVars.moleFraction(wPhaseIdx,Fe3Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
		// else
		// {
		// if(ferrdiss >  volFracFerrohydrite * volVars.solidComponentDensity(fPhaseIdx) / dt)
		// {
		// 	ferrdiss =  volFracFerrohydrite * volVars.solidComponentDensity(fPhaseIdx) / dt;
		// }
		// }
		
        // if (ferOmegaApprox_ > 1)
        // {		 
		// ferrprecdt =  volVars.moleFraction(wPhaseIdx,Fe3Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / ferrprec;
		// ferrdt = ferrprecdt;
		// }
		// else
		// {
		// ferrdissdt =  volFracFerrohydrite * volVars.solidComponentDensity(fPhaseIdx) / ferrdiss;
		// ferrdt =  ferrdissdt;	
		// }		
		
	    // ferrdiss = 0;
		// ferrprec = 0;

		
        // ProtoImogolite
        Scalar prosp = const7(pressure_, temperature_);
        Scalar proOmegaApprox_ = pow(mAl,2) * pow(mSiO2,1) / pow (mH,6) /prosp;
        // printf("The value of mAl is: %.12e\n", mAl);
		// printf("The value of mSiO2 is: %.12e\n", mSiO2);
		// printf("The value of mH is: %.12e\n", mH);
		// printf("The value of prosp is: %.12e\n", prosp);
		// printf("The value of proOmegaApprox_ is: %.12e\n", proOmegaApprox_);
        Scalar proAw0 = proAsw_ * initialvolFracProtoImogolite * volVars.solidComponentDensity(pPhaseIdx);
		Scalar proAwf = proAsw_ * volFracProtoImogolite * volVars.solidComponentDensity(pPhaseIdx);
 		Scalar proAwcd = proAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracProtoImogolite*volFracProtoImogolite));   // TODO Asw should be a function of Sw, too!
 		Scalar proAwcp = proAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!

		Scalar proAwd = proAwcd;
		Scalar proAwp = proAwcp;		
	
	
        if (proAwp < 0 || std::isnan(proAwp))
        {
        std::cout<< "proAwp = "<<proAwp<<std::endl;
        proAwp = 0;
        std::cout<< "proAwp, corrected = "<<proAwp<<std::endl;
        }
		
        if (proAwd < 0 || std::isnan(proAwd))
        {
        std::cout<< "proAwd = "<<proAwd<<std::endl;
        proAwd = 0;
        std::cout<< "proAwd, corrected = "<<proAwd<<std::endl;
        }
			

        Scalar prordiss = 0;
        Scalar prorprec = 0;
		Scalar pror = 0;
		
        if (proOmegaApprox_ > 1)
        {
			pror = proAwp * pow(10,prorc_) * prop_ * pow(abs(1-pow(proOmegaApprox_,1/prosigma_)),probeta_);
            prordiss = 0;
            prorprec = pror;//[mol/kg s]
        }
        else
        {
			pror = proAwd * pow(10,prorc_) * prop_ * pow(abs(1-pow(proOmegaApprox_,1/prosigma_)),probeta_);
			prordiss = pror; //[mol/kg s]
            prorprec = 0;
			if (std::isnan(pror))
			{
			}
        }
		
		
        // if (proOmegaApprox_ > 1)
        // {
		// if (volVars.moleFraction(wPhaseIdx,SiO2Idx) > volVars.moleFraction(wPhaseIdx,AlIdx))
		// {	
		// if (prorprec > volVars.moleFraction(wPhaseIdx,AlIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     // printf("The value of prorprec1 is: %.12e\n", prorprec);
		// 	prorprec =  volVars.moleFraction(wPhaseIdx,AlIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        //     // printf("The value of prorprec2 is: %.12e\n", prorprec);
		// }
		// }
		// else
		// {	
		// if (prorprec > volVars.moleFraction(wPhaseIdx,SiO2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     prorprec =  volVars.moleFraction(wPhaseIdx,SiO2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
		// }
		// else
		// {
		// if (prordiss >  volFracProtoImogolite * volVars.solidComponentDensity(pPhaseIdx) / dt)
		// {
		// prordiss =  volFracProtoImogolite * volVars.solidComponentDensity(pPhaseIdx) / dt;
		// }
		// }	

		
        //   if (proOmegaApprox_ > 1)
        //   {
		//    if (volVars.moleFraction(wPhaseIdx,SiO2Idx) > volVars.moleFraction(wPhaseIdx,AlIdx))
		//    {
		//    prorprecdt =  volVars.moleFraction(wPhaseIdx,AlIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / prorprec;
		//    }
		//    else
		//    {
		//    prorprecdt =  volVars.moleFraction(wPhaseIdx,SiO2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / prorprec;
		//    }
		//    prordt = prorprecdt;	
		//  }
		// else
		// {
		//   prordissdt =  volFracProtoImogolite * volVars.solidComponentDensity(pPhaseIdx) / prordiss;
		//   prordt =  prordissdt;	
		// }
		
		
        // Birnessite
        Scalar birsp = const8(pressure_, temperature_);
        Scalar birOmegaApprox_ = pow(mMn,8) * pow(mO2,3) / pow (mH,16) /birsp;
		Scalar birAw0 = birAsw_ * initialvolFracBirnessite * volVars.solidComponentDensity(bPhaseIdx);
 		Scalar birAwf = birAsw_ * volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx);
 		Scalar birAwcp = birAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!
 		Scalar birAwcd = birAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracBirnessite)*(volFracBirnessite));   // TODO Asw should be a function of Sw, too!

		Scalar birAwd = birAwcd;
		Scalar birAwp = birAwcp;		

        if (birAwd < 0 || std::isnan(birAwd))
        {
        std::cout<< "birAwd = "<<birAwd<<std::endl;
        birAwd = 0;
        std::cout<< "birAwd, corrected = "<<birAwd<<std::endl;
        }

        if (birAwp < 0 || std::isnan(birAwp))
        {
        std::cout<< "birAwp = "<<birAwp<<std::endl;
        birAwp = 0;
        std::cout<< "birAwp, corrected = "<<birAwp<<std::endl;
        }		

        Scalar birrdiss = 0;
        Scalar birrprec = 0;
		
        if (birOmegaApprox_ > 1)
        {
			Scalar birr = birAwp * pow(10,birrc_) * birp_ * pow(abs(1-pow(birOmegaApprox_,1/birsigma_)),birbeta_);
            birrdiss = 0;
            birrprec = birr;//[mol/kg s]
        }
        else
        {
		    Scalar birr = birAwd * pow(10,birrc_) * birp_ * pow(abs(1-pow(birOmegaApprox_,1/birsigma_)),birbeta_);
			birrdiss = birr; //[mol/kg s]
            birrprec = 0;
			if (std::isnan(birr))
			{
			}
        }
		
		// if (birOmegaApprox_ > 1)
		// {
		// if (volVars.moleFraction(wPhaseIdx,O2Idx) > volVars.moleFraction(wPhaseIdx,MnIdx))
		// {	
		// if(birrprec > volVars.moleFraction(wPhaseIdx,MnIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     birrprec =  volVars.moleFraction(wPhaseIdx,MnIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
		// else
		// {	
		// if(birrprec > volVars.moleFraction(wPhaseIdx,O2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     birrprec =  volVars.moleFraction(wPhaseIdx,O2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
		// }
		// else
		// {	
	    //   if(birrdiss >  volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / dt)
		//   {
		//   birrdiss =  volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / dt;
		//   }
		// }	
		
		// if (birOmegaApprox_ > 1)
		// {
		//   if (volVars.moleFraction(wPhaseIdx,O2Idx) > volVars.moleFraction(wPhaseIdx,MnIdx))
		//   {	
        //       birrprecdt =  volVars.moleFraction(wPhaseIdx,MnIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / birrprec;
        //   }
		//   else
		//   {	
        //       birrprecdt =  volVars.moleFraction(wPhaseIdx,O2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / birrprec;
        //   }	
		//   birrdt =  birrprecdt;
		// }
        // else
		// {		
		//   birrdissdt =  volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / birrdiss;
		//   birrdt =  birrdissdt;
		// }
		
        // Hydroxyapatite
        Scalar hydsp = const9(pressure_, temperature_);
        Scalar hydOmegaApprox_ = pow(mCa,5) * pow(mHPO4,3) / pow (mH,4) /hydsp;
		Scalar hydAw0 = hydAsw_ * initialvolFracHydroxyapatite * volVars.solidComponentDensity(hPhaseIdx);
 		Scalar hydAwf = hydAsw_ * volFracHydroxyapatite * volVars.solidComponentDensity(hPhaseIdx);
 		Scalar hydAwcd = hydAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracHydroxyapatite)*(volFracHydroxyapatite));   
 		Scalar hydAwcp = hydAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   

		Scalar hydAwd = hydAwcd;
		Scalar hydAwp = hydAwcp;

        if (hydAwd < 0 || std::isnan(hydAwd))
        {
        std::cout<< "hydAwd = "<<hydAwd<<std::endl;
        hydAwd = 0;
        std::cout<< "hydAwd, corrected = "<<hydAwd<<std::endl;
        }

        if (hydAwp < 0 || std::isnan(hydAwp))
        {
        std::cout<< "hydAwp = "<<hydAwp<<std::endl;
        hydAwp = 0;
        std::cout<< "hydAwp, corrected = "<<hydAwp<<std::endl;
        }
			

        Scalar hydrdiss = 0;
        Scalar hydrprec = 0;
		
        if (hydOmegaApprox_ > 1)
        {
		    Scalar hydr = hydAwp * pow(10,hydrc_) * hydp_ * pow(abs(1-pow(hydOmegaApprox_,1/hydsigma_)),hydbeta_);
			// printf("Precipitation");   
            hydrdiss = 0;
            hydrprec = hydr;//[mol/kg s]
        }
        else
        {
		    Scalar hydr = hydAwd * pow(10,hydrc_) * hydp_ * pow(abs(1-pow(hydOmegaApprox_,1/hydsigma_)),hydbeta_);
			// printf("Dissolution");   
			hydrdiss = hydr; //[mol/kg s]
            hydrprec = 0;
			if (std::isnan(hydr))
			{
			}
        }

		// if (hydOmegaApprox_ > 1)
        // {
		// if (volVars.moleFraction(wPhaseIdx,HPO4Idx) > volVars.moleFraction(wPhaseIdx,CaIdx))
		// {		
		// if(hydrprec > volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     hydrprec =  volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
		// else
		// {		
		// if(hydrprec > volVars.moleFraction(wPhaseIdx,HPO4Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     hydrprec =  volVars.moleFraction(wPhaseIdx,HPO4Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
		// }
		// else
		// {
		// if(hydrdiss > volFracHydroxyapatite * volVars.solidComponentDensity(hPhaseIdx) / dt)
        // {
		//   hydrdiss =  volFracHydroxyapatite * volVars.solidComponentDensity(hPhaseIdx) / dt;
		// }
		// }	
		
        // if (hydOmegaApprox_ > 1)
        // {
		//   if (volVars.moleFraction(wPhaseIdx,HPO4Idx) > volVars.moleFraction(wPhaseIdx,CaIdx))
        //   {
        //       hydrprecdt =  volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / hydrprec;
        //   }
		//   else
        //   {
        //       hydrprecdt =  volVars.moleFraction(wPhaseIdx,HPO4Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / hydrprec;
        //   }
		//   hydrdt =  hydrprecdt;
		// }
		// else
		// {
		//   hydrdissdt =  volFracHydroxyapatite * volVars.solidComponentDensity(hPhaseIdx) / hydrdiss;
        //   hydrdt =  hydrdissdt;
		// }
		
        // Sepiolite
        Scalar sepsp = const10(pressure_, temperature_);
        Scalar sepOmegaApprox_ = pow(mMg,4) * pow(mSiO2,6) / pow (mH,8) /sepsp;
		Scalar sepAw0 = sepAsw_ * initialvolFracSepiolite * volVars.solidComponentDensity(sPhaseIdx);
 		Scalar sepAwp = sepAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   
        Scalar sepAwd = sepAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracSepiolite)*(volFracSepiolite));   

		
		if (sepAwd < 0 || std::isnan(sepAwd))
        {
        std::cout<< "sepAwd = "<<sepAwd<<std::endl;
        sepAwd = 0;
        std::cout<< "sepAwd, corrected = "<<sepAwd<<std::endl;
        }
	
		if (sepAwp < 0 || std::isnan(sepAwp))
        {
        std::cout<< "sepAwp = "<<sepAwp<<std::endl;
        sepAwp = 0;
        std::cout<< "sepAwp, corrected = "<<sepAwp<<std::endl;
        }
	

        Scalar seprdiss = 0;
        Scalar seprprec = 0;
		
        if (sepOmegaApprox_ > 1)
        {
		    Scalar sepr = sepAwp * pow(10,seprc_) * sepp_ * pow(abs(1-pow(sepOmegaApprox_,1/sepsigma_)),sepbeta_);
            seprdiss = 0;
            seprprec = sepr;//[mol/kg s]
        }
        else
        {
		    Scalar sepr = sepAwd * pow(10,seprc_) * sepp_ * pow(abs(1-pow(sepOmegaApprox_,1/sepsigma_)),sepbeta_);
			seprdiss = sepr; //[mol/kg s]
            seprprec = 0;
			if (std::isnan(sepr))
			{
			}
        }
		
        // if (sepOmegaApprox_ > 1)
        // {
	    // if (volVars.moleFraction(wPhaseIdx,MgIdx) > volVars.moleFraction(wPhaseIdx,SiO2Idx))
		// {
		// if(seprprec > volVars.moleFraction(wPhaseIdx,SiO2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     seprprec =  volVars.moleFraction(wPhaseIdx,SiO2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
        // else	
		// {
		// if(seprprec > volVars.moleFraction(wPhaseIdx,MgIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt)
        // {
        //     seprprec =  volVars.moleFraction(wPhaseIdx,MgIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / dt;
        // }
		// }
		// }
		// else
		// {
		// if(seprdiss > volFracSepiolite * volVars.solidComponentDensity(sPhaseIdx) / dt)
        // {
        //     seprdiss =  volFracSepiolite * volVars.solidComponentDensity(sPhaseIdx) / dt;
		// 	 
        // }
		// }		
		
        // if (sepOmegaApprox_ > 1)
        // {
		//   if (volVars.moleFraction(wPhaseIdx,MgIdx) > volVars.moleFraction(wPhaseIdx,SiO2Idx))
        //   {
        //       seprprecdt =  volVars.moleFraction(wPhaseIdx,SiO2Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / seprprec;
        //   }
        //   else	
        //   {
        //       seprprecdt =  volVars.moleFraction(wPhaseIdx,MgIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx) / seprprec;
        //   }
		//   seprdt =  seprprecdt;
		// }
        // else
        // {
        // seprdissdt =  volFracSepiolite * volVars.solidComponentDensity(sPhaseIdx) / seprdiss;
        // seprdt =  seprdissdt;
		// }	
	

		
		// ferrprec= 0;
		// ferrdiss= 0;
		// prorprec= 0; //issue lies here
		// prordiss= 0;
		// birrprec= 0;
		// birrdiss= 0;
		// hydrprec= 0;
		// hydrdiss= 0;
		// seprprec= 0;
		// seprdiss= 0;
	
		// printf("The value of glrprec is: %.12e\n", glrprec); 
		// printf("The value of glrdiss is: %.12e\n", glrdiss);
		// printf("The value of ferrprec is: %.12e\n", ferrprec); 
		// printf("The value of ferrdiss is: %.12e\n", ferrdiss); 
		// printf("The value of prorprec is: %.12e\n", prorprec); 
		// printf("The value of prordiss is: %.12e\n", prordiss); 
		// printf("The value of birrprec is: %.12e\n", birrprec); 
		// printf("The value of birrdiss is: %.12e\n", birrdiss); 
		// printf("The value of hydrprec is: %.12e\n", hydrprec); 
		// printf("The value of hydrdiss is: %.12e\n", hydrdiss); 
		// printf("The value of seprprec is: %.12e\n", seprprec); 
		// printf("The value of seprdiss is: %.12e\n", seprdiss); 
		
        q[wCompIdx] += 0;// 2.05 *(- glrprec + glrdiss)+13*(- birrprec + birrdiss) +1*(- hydrprec + hydrdiss)+11*(- seprprec + seprdiss)+ 2 * (- ferrprec + ferrdiss) +5*(- prorprec + prordiss) ;// +(- calrprec + calrdiss);
        q[nCompIdx] += 0;
        q[CO2aqIdx] += 0;// (- calrprec + calrdiss);		
        q[NaIdx] += 0.26*(- glrprec + glrdiss);
        q[ClIdx] += 0;
        q[O2Idx] += 3*(- birrprec + birrdiss);
        q[CaIdx] += 0.44 *(- glrprec + glrdiss)+5*(- hydrprec + hydrdiss);//+ 1 *(- calrprec + calrdiss);
		q[Fe2Idx] += 0.38*(- glrprec + glrdiss)+1* (- ferrprec + ferrdiss) ; // added by du
        q[SiO2Idx] += 1.8*(- glrprec + glrdiss)+6*(- seprprec + seprdiss)+1*(- prorprec + prordiss) ; // added by du
        q[TiOH4Idx] += 0.07*(- glrprec + glrdiss);// added by du
        q[AlIdx]+= 0.62*(- glrprec + glrdiss)+2*(- prorprec + prordiss) ; // added by du
        q[MnIdx]+= 0.01*(- glrprec + glrdiss)+8*(- birrprec + birrdiss); // added by du
        q[MgIdx]+=0.3*(- glrprec + glrdiss)+4*(- seprprec + seprdiss); // added by du
        q[KIdx]+= 0.06*(- glrprec + glrdiss); // added by du
        q[HPO4Idx]+= 0.03*(- glrprec + glrdiss)+3*(- hydrprec + hydrdiss); // added by du
		q[HIdx]+= 4.38*(+ glrprec - glrdiss)+16*(birrprec - birrdiss)+4*(hydrprec - hydrdiss)+8*(seprprec - seprdiss)+ 2* (ferrprec - ferrdiss) +6*(prorprec - prordiss) ; // added by du +2*(calrprec - calrdiss) 
        q[phiGlassIdx] += (glrprec - glrdiss);
		q[phiFerrohydriteIdx] += (ferrprec - ferrdiss) ;
		q[phiProtoImogoliteIdx] += (prorprec - prordiss) ;
		q[phiBirnessiteIdx] += (birrprec - birrdiss);
		q[phiHydroxyapatiteIdx] += (hydrprec - hydrdiss);
		q[phiSepioliteIdx] += (seprprec - seprdiss);
		
		// printf("The value of q[HIdx] is: %.7e\n", q[HIdx]);
		// printf("The value of mH is: %.7e\n", mH);		
		// printf("The value of q[AlIdx] is: %.7e\n", q[AlIdx]);
		// printf("The value of mAl is: %.7e\n", mAl); 
		// printf("The value of q[Fe2Idx] is: %.7e\n", q[Fe2Idx]);
		// printf("The value of mFe2is: %.7e\n", mFe2); 		
		// printf("The value of q[phiGlassIdx] is: %.7e\n", q[phiGlassIdx]);		
        // printf("The value of q[phiProtoImogoliteIdx] is: %.7e\n", q[phiProtoImogoliteIdx]);
		// printf("The value of q[phiFerrohydriteIdx] is: %.7e\n", q[phiFerrohydriteIdx]);		
        // printf("The value of q[phiBirnessiteIdx] is: %.7e\n", q[phiBirnessiteIdx]);
        // printf("The value of q[phiHydroxyapatiteIdx] is: %.7e\n", q[phiHydroxyapatiteIdx]);
        // printf("The value of q[phiSepioliteIdx] is: %.7e\n", q[phiSepioliteIdx]);
		// 
		Scalar Hdt = 0;
		
	    if (0 > q[AlIdx])
		{
		prorprecdt = mAl / abs(q[AlIdx]);
		newdt = prorprecdt;
		}
	    // else if (0 > q[HIdx])
		// {
		// Hdt = mH / abs(q[HIdx]);
		// newdt = Hdt;
		// }
		else 
		{newdt = 0;}


		return std::make_pair(q, newdt);
    }

private:
    //Newton Solver which returns true if convergence is reached and false if not.
    // x(i+1) = x(i) - f(i)/df(i) = x(i) + h(i)

    bool newton1D(Scalar &xVar, const Scalar tolAbs, const int maxIter)
    {
        // printf("reached newton"); unreached
		/*NEWTON*/
        bool converge = false;
          Scalar eps = 1.e-3;
          Scalar eps2 = 1.e-10;
          Scalar b=0;
          Scalar c=100;
          Scalar r;
          Scalar pHc = - log(xVar);
          Scalar pHb = pHc+eps;
          Scalar Hb,Hc;
          Scalar CO3l,CO3r,CO3b, CO3c;
          // Scalar NH3l,NH3r,NH3b, NH3c;
          Scalar error =100;
          iter_ = 0;
          int i = 0;
//        *HCO3 = *NH4 = CO3c = 0.;
          Scalar oh,hco3,co3,co2aq;
          oh=hco3=co3=co2aq=0;
          // Scalar oh,co2;
          // oh=co2=0;
          while (absolute(c) > tolAbs)
            {
              Hb = pow(10.,-pHb);
              CO3l = 0.; CO3r = ctot_;
              CO3b = (CO3l+CO3r)/2.;
              while (absolute(error)>1.E-11)
            {
              CO3r = CO3b + eps2;
              r = ctot_ - (Hb * CO3r / k2_)  - (Hb * (Hb*CO3r/k2_) / k1_) - CO3r ;
              error = ctot_ - (Hb * CO3b / k2_)  - (Hb * (Hb*CO3b/k2_) / k1_) - CO3b;
              CO3b = CO3b - (eps2*error)/(r-error);
              i++;  if (i>1.E2) break;
            }
			
              hco3 = Hb * CO3b / k2_;
              co2aq = Hb * hco3 / k1_;
			  // co2g = co2aq / kh_;

            //   NH3l = 0.;
            //   NH3r = totalnh_;
            //   error =100.;
            //   i = 0;
            //   NH3b = (NH3l+NH3r)/2.;
            //   while (absolute(error)>1.E-11)
            // {
            //   NH3r = NH3b + eps2;
            //   r = totalnh_ - Hb * NH3r / ka_ - NH3r;
            //   error = totalnh_ - Hb * NH3b / ka_ - NH3b;
            //   NH3b = NH3b - (eps2*error)/(r-error);
            //   i++;  if (i>1.E2) break;
            // }
            //   nh4 = Hb * NH3b / ka_;

              oh = kw_ / Hb;

              //b = - Hb + 2*CO3b + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_ - 2.*fe2_;
              b = - Hb + 2*CO3b + hco3 + oh - na_ + cl_ - 2.*ca_ - 2.*fe2_;
              pHc = pHc - (eps*c)/(b-c);

              Hc = pow(10.,-pHc);
              CO3l = 0.;
              CO3r = ctot_;
              error =100.; i = 0;
              CO3c = (CO3l+CO3r)/2.;
              while (absolute(error)>1.E-11)
            {
              CO3r = CO3c + eps2;
              r = ctot_ - (Hc * CO3r / k2_)  - (Hc * (Hc*CO3r/k2_) / k1_) - CO3r ;
              error = ctot_ - (Hc * CO3c / k2_)  - (Hc * (Hc*CO3c/k2_) / k1_) - CO3c ;
              CO3c = CO3c - (eps2*error)/(r-error);
              i++; if (i>1.E2) break;
            }
              hco3 = Hc * CO3c / k2_;
              co2aq = Hc * hco3 / k1_;
			  // co2g = co2aq / kh_;

            //   NH3l = 0.;
            //   NH3r = totalnh_;
            //   error =100.;
            //   i = 0;
            //   NH3c = (NH3l+NH3r)/2.;
            //   while (absolute(error)>1.E-11)
            // {
            //   NH3r = NH3c + eps2;
            //   r = totalnh_ - Hc * NH3r / ka_ - NH3r;
            //   error = totalnh_ - Hc * NH3c / ka_ - NH3c;
            //   NH3c = NH3c - (eps2*error)/(r-error);
            //   i++; if (i>1.E2) break;
            // }
            //   nh4 = Hc * NH3c / ka_;

              oh = kw_ / Hc;
              // c = - Hc + 2*CO3c + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;
              c = - Hc + 2*CO3c + hco3 + oh - na_ + cl_ - 2.*ca_- 2.*fe2_;
              pHb = pHc+eps;
              iter_+=1;
              if (iter_>maxIter || isnan(error) || isnan(c))
            {
              /*sprintf(buf, "Bisection pH: %4.2f \n", pHc);
              UserWrite(buf);*/
              break;
            }
            }
          h_ = Hc;
          oh_ = kw_ / Hc;
          // nh4_ = Hc * NH3c / ka_; added by du
		  // // nh4_ = 0;
          co3_ = CO3c;
          hco3_ = Hc * CO3c / k2_;
          co2aq_ = Hc * hco3 / k1_;
		  // co2g_ = Hc * hco3 / k1_ / kh_;
		  
//        (this->*funcPtr)(xVar);
//        Scalar h = -fdf_[0]/fdf_[1]; // h = x(i) - x(i-1)
//        Scalar hLast = h*0.5; //initial Step
//        iter_ = 0;
//        bool converge = false;
//        if (std::isnan(h))
//        {
//            return converge = false;
//        }
//
//        while(absolute(h) > tolAbs || absolute(h/hLast)  > 1 + tolRel)
//        {
//            if(iter_ > maxIter){break;}
//
//            if(iter_ > 0)
//            {
//                (this->*funcPtr)(xVar);
//                hLast = h;
//                h = -fdf_[0]/fdf_[1];
//            }
//            if (std::isnan(h))
//            {
//                return converge = false;
//            }
//
//            xVar = xVar + h;
//            iter_ = iter_ + 1;
//        }
        if(Hc < 0.0) {return converge = false;}
        if(iter_ <= maxIter) {converge = true;}
        return converge;

    }

    bool newton1D(Scalar &xVar, void (ThisType::*funcPtr)(Scalar), const Scalar tolAbs, const Scalar tolRel, const int maxIter)
    {
        // printf("reached newton2"); // reached
        if (!Valgrind::CheckDefined(xVar))
        {
            std::cout << "----!Valgrind::CheckDefined(xVar) in chemistry \n";
            DUNE_THROW(Dune::InvalidStateException, "xVar is not defined.");
        }
        (this->*funcPtr)(xVar);

        Scalar h = -fdf_[0]/fdf_[1]; // h = x(i) - x(i-1)
        Scalar hLast = h*0.5; //initial Step
        iter_ = 0;
        bool converge = false;
        if (std::isnan(h))
        {
            return converge = false;
        }

        while(absolute(h) > tolAbs || absolute(h/hLast)  > 1 + tolRel)
        {
            if(iter_ > maxIter){break;}

            if(iter_ > 0)
            {
                (this->*funcPtr)(xVar);
                hLast = h;
                h = -fdf_[0]/fdf_[1];
            }
            if (std::isnan(h))
            {
                return converge = false;
            }

            xVar = xVar + h;
            iter_ = iter_ + 1;
        }
        if(xVar < 0.0) {return converge = false;}
        if(iter_ <= maxIter) {converge = true; newtonOrBisection_ = true; }
        return converge;

    }

    //Bisection Method Solver returns true if convergence is reached and false if not.
    //xVar is the variable for which the system is solved
    //funcPtr is the pointer to the function which is to be solved
    //a0 is the lower starting value, b0 is the upper starting value. The root must be inside the interval [a0, b0]
    //tol is the stopping critium a-b
    bool bisection1D(Scalar &xVar, void (ThisType::*funcPtr)(Scalar), const Scalar a0, const Scalar b0, const Scalar tol)
    {
        // printf("reached bisection1D"); // reached
		Scalar iterNo = 0;
        int maxIter = 200;
        bool converge = false;
        int sfb, sfx;

        Scalar a = a0;
        Scalar b = b0;
        (this->*funcPtr)(b);
        sfb = sign(fdf_[0]);

        while(b-a > tol)
        {
            if(iterNo > maxIter)
            {
                return converge;
            }
            xVar = (b + a)/2;
            (this->*funcPtr)(xVar);
            sfx = sign(fdf_[0]);
            iterNo = iterNo + 1;
            if (sfx == 0)
                break;
            else
                {
                    if(sfx == sfb)
                    {
                        b = xVar;
                    }
                    else
                    {
                        a = xVar;
                    }
                }
        }
        newtonOrBisection_ = false;
        converge = true;
        return converge;
    }

    bool bisection1D(const Scalar tol)
    {
        // printf("reached bisection1D 2"); unreached
		bool converge = false;
        Scalar eps = 1.e-3;
        Scalar eps2 = 1.e-10;
        Scalar pHc = 7.;
        Scalar pHa = -1.;
        Scalar pHb = 15.;
        Scalar Ha,Hb,Hc;
        // Scalar CO3r,CO3l,CO3a,CO3b,CO3c;
        // Scalar NH3l,NH3r,NH3a,NH3b,NH3c;
        Scalar c=100.;
        Scalar a,b;
        Scalar error=100;
        Scalar r;
        iter_ = 0;
        int i = 0;
        // Scalar oh,hco3,nh4,co3,co2;
        // oh=hco3=nh4=co3=co2=0;
        Scalar oh,co2;
        oh=co2=0;
        while (absolute(c) > tol)
    {
      Ha =pow(10.,-pHa);
      // CO3l = 0.;
      // CO3r = cTot_;
      // error =100.;
      // i = 0;
      // CO3a = (CO3l+CO3r)/2.;
      // while (absolute(error)>1.E-10)
      //   {
      //     CO3r = CO3a + eps2;
      //     r = cTot_ - (Ha * CO3r / k2_)  - (Ha * (Ha*CO3r/k2_) / k1_) - CO3r;
      //     error = cTot_ - (Ha * CO3a / k2_)  - (Ha * (Ha*CO3a/k2_) / k1_) - CO3a;
      //     CO3a = CO3a - (eps2*error)/(r-error);
      //     i++;  if (i>1.E2) break;
      //   }
	  // 
      // hco3 = Ha * CO3a / k2_;
      // co2 = Ha * hco3 / k1_;

      // NH3l = 0.;
      // NH3r = totalnh_;
      // error =100.;
      // i=0;
      // NH3a = (NH3l+NH3r)/2.;
      // while (absolute(error)>1.E-10)
      //   {
      //     NH3r = NH3a + eps2;
      //     r = totalnh_ - Ha * NH3r / ka_ - NH3r;
      //     error = totalnh_ - Ha * NH3a / ka_ - NH3a;
      //     NH3a = NH3a - (eps2*error)/(r-error);
      //     i++;  if (i>1.E2) break;
      //   }
      // nh4 = Ha * NH3a / ka_;

      oh = kw_ / Ha;
      // a = - Ha + 2*CO3a + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;
      a = - Ha + oh - na_ + cl_ - 2.*ca_- 2.*fe2_;
	  
      Hb = pow(10.,-pHb);
      // CO3l = 0.;
      // CO3r = cTot_;
      error =100.;
      i = 0;
      // CO3b = (CO3l+CO3r)/2.;
      // while (absolute(error)>1.E-11)
      //   {
      //     CO3r = CO3b + eps2;
      //     r = cTot_ - (Hb * CO3r / k2_)  - (Hb * (Hb*CO3r/k2_) / k1_) - CO3r;
      //     error = cTot_ - (Hb * CO3b / k2_)  - (Hb * (Hb*CO3b/k2_) / k1_) - CO3b;
      //     CO3b = CO3b - (eps2*error)/(r-error);
      //     i++;  if (i>1.E2) break;
      //   }
	  // 
      // hco3 = Hb * CO3b / k2_;
      // co2 = Hb * hco3 / k1_;

      // NH3l = 0.;
      // NH3r = totalnh_;
      // error =100.;
      // i = 0;
      // NH3b = (NH3l+NH3r)/2.;
      // while (absolute(error)>1.E-11)
      //   {
      //     NH3r = NH3b + eps2;
      //     r = totalnh_ - Hb * NH3r / ka_ - NH3r;
      //     error = totalnh_ - Hb * NH3b / ka_ - NH3b;
      //     NH3b = NH3b - (eps2*error)/(r-error);
      //     i++;  if (i>1.E2) break;
      //   }
      // nh4 = Hb * NH3b / ka_;

      oh = kw_ / Hb;
      // b = - Hb + 2*CO3b + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;
      b = - Hb +  oh - na_ + cl_ - 2.*ca_- 2.*fe2_;
      pHc = (pHa + pHb)/2.;

      Hc = pow(10.,-pHc);
      // CO3l = 0.;
      // CO3r = cTot_;
      // error =100.;
      // i = 0;
      // CO3c = (CO3l+CO3r)/2.;
      // while (absolute(error)>1.E-10)
      //   {
      //     CO3r = CO3c + eps2;
      //     r = cTot_ - (Hc * CO3r / k2_)  - (Hc * (Hc*CO3r/k2_) / k1_) - CO3r;
      //     error = cTot_ - (Hc * CO3c / k2_)  - (Hc * (Hc*CO3c/k2_) / k1_) - CO3c;
      //     CO3c = CO3c - (eps2*error)/(r-error);
      //     i++; if (i>1.E2) break;
      //   }
	  // 
      // hco3 = Hc * CO3c / k2_;
      // co2 = Hc * hco3 / k1_;

    //     NH3l = 0.;
    //     NH3r = totalnh_;
    //     error =100.;
    //     i = 0;
    //     NH3c = (NH3l+NH3r)/2.;
    //     while (absolute(error)>1.E-11)
    // {
    //   NH3r = NH3c + eps2;
    //   r = totalnh_ - Hc * NH3r / ka_ - NH3r;
    //   error = totalnh_ - Hc * NH3c / ka_ - NH3c;
    //   NH3c = NH3c - (eps2*error)/(r-error);
    //   i++; if (i>1.E2) break;
    // }
    //   nh4 = Hc * NH3c / ka_;

      oh = kw_ / Hc;
      // c = - Hc + 2*CO3c + hco3 + oh - nh4 - na_ + cl_ - 2.*ca_- 2.*fe2_;
      c = - Hc + oh - na_ + cl_ - 2.*ca_- 2.*fe2_;
      if (a*c<0.) pHb = pHc;
      else pHa = pHc;
      iter_+=1;

    }

      h_ = Hc;
      oh_ = kw_ / Hc;
      // nh4_ = Hc * NH3c / ka_;
      // co3_ = CO3c;
      // hco3_ = Hc * CO3c / k2_;
      // co2_ = Hc * hco3 / k1_;

        converge = true;
        return converge;
    }



    // void H_Conly(Scalar activityH)
    //     {
	// 
    //     h_ = activityH;
    //     oh_ = kw_/h_;
	// 	hco3_ = k1_*co2aq_/(h_);
    //     co3_ = k1_*k2_*co2aq_/(h_*h_);
	// 	fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(kw_/h_, 4)), 0.25);
	// 	
	// 	
    //     //Solve the function
    //     Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ + 3.*fe3_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
    //     //Solve the derivative df/d(activityH)
    //     Scalar eps = 1e-8;
    //     Scalar xRight = h_ + eps*h_; // x + dx
    //     Scalar xLeft = h_ - eps*h_; // x - dx
    //     Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - k1_*co2aq_/xRight - 2*k1_*k2_*co2aq_/(xRight*xRight) - cl_ + 3*pow((k3_ * pow(fe2_, 4) * o2_ / pow(kw_/xRight, 4)), 0.25)+3.*al_ +2.*mg_ +k_ + 2.*mn_ +4.*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
    //     Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft - k1_*co2aq_/xLeft - 2*k1_*k2_*co2aq_/(xLeft*xLeft) - cl_ + 3*pow((k3_ * pow(fe2_, 4) * o2_ / pow(kw_/xLeft, 4)), 0.25)+3.*al_ +2.*mg_ +k_ + 2.*mn_ +4.*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx) f(x-dx)
    //     Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx
	// 
	// 
    //     fdf_[0] = f;
    //     fdf_[1] = df;
    //  }
    void H_Conly(Scalar activityCO2)
        {

        co2aq_ = activityCO2;
        oh_ = kw_/h_;
		hco3_ = k1_*co2aq_/(h_);
        co3_ = k1_*k2_*co2aq_/(h_*h_);
		fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(kw_/h_, 4)), 0.25);
		
		
        //Solve the function
        Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ + 3.*fe3_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = co2aq_ + eps*co2aq_; // x + dx
        Scalar xLeft = co2aq_ - eps*co2aq_; // x - dx
        Scalar fRight = na_ + h_ + 2*ca_ + 2*fe2_ - oh_ - k1_*xRight/h_ - 2*k1_*k2_*xRight/(h_*h_) - cl_ + 3*fe3_+3.*al_ +2.*mg_ +k_ + 2.*mn_ +4.*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = na_ + h_ + 2*ca_ + 2*fe2_ - oh_ - k1_*xLeft/h_ - 2*k1_*k2_*xLeft/(h_*h_) - cl_ + 3*fe3_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4.*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar df = (fRight - fLeft)/2/eps/co2aq_; // {f(x+dx) - f(x-dx)}/2dx


        fdf_[0] = f;
        fdf_[1] = df;
     }
    //Value of numerical derivative at xVar
    /*static*/ Scalar equationNumDeri(Scalar xVar)
    {
        Scalar eps = 1e-8;
        Scalar xRight = xVar + eps*xVar; // x + dx
        Scalar xLeft = xVar - eps*xVar; // x - dx
        Scalar fRight = equationValue(xRight); // f(x+dx)
        Scalar fLeft = equationValue(xLeft); // f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/xVar; // {f(x+dx) - f(x-dx)}/2dx
        return df;
    }



    Scalar absolute(Scalar x)
    {
        if(x<0.0)
        {
            return x*(-1);
        }
        else return x;
    }

    Scalar sign(Scalar x)
    {
        if(x > 0.0)
        {
           return 1;
        }
        else if (x < 0.0)
        {
           return -1;
        }
        else
        {
            return 0.0;
        }
    }



    int iter_; //Number of iterations the Newton solver needs until convergence
    Scalar pressure_;
    Scalar temperature_;
    Scalar salinity_;
    Scalar h2o_;
    Scalar co2aq_;
    Scalar co2aqtotal_;
	Scalar co2g_;
    Scalar o2_;
    Scalar hco3_;
    Scalar co3_;
    Scalar oh_;
    Scalar h_;
    Scalar htotal_;
    Scalar h2_;
    Scalar ca_;
    Scalar na_;
    Scalar cl_;
    Scalar fe2_;
	Scalar fe2total_;
    Scalar fe3_;
    Scalar sio2_;
    Scalar tioh4_;
    Scalar al_;
    Scalar mn_;
    Scalar mg_;
    Scalar k_;
    Scalar hpo4_;
 	
    Scalar initH_;
    Scalar initCO2_;
    Scalar initCl_;
    Scalar initCO2aq_;
    Scalar initCO3_;
    Scalar ionicStrength_;
    Scalar ctot_;
    Scalar gammaH_;
    Scalar gammaCO2_;
    Scalar gammaCa_;
    Scalar gammaOH_;
    SolVector fdf_; //Solution vector for the newtons solver every equation f solved by the newton solver for an unknown x
    // has to store f(x) in fdf_[0] and df/dx in fdf[1]
    Vector molality_;
    Vector molarity_;
    Vector charge_;
    Scalar x_;
    Scalar y_;
    Scalar k1_;
    Scalar k2_;
    Scalar k3_;
    Scalar k4_;
    Scalar k5_;
    Scalar k6_;
    Scalar k7_;
    Scalar k8_;
    Scalar k9_;
    Scalar k10_;
	Scalar k11_;
    Scalar kw_;
    Scalar ka_;
    Scalar apparentk1_;
    Scalar apparentk2_;
    Scalar apparentka_;
	
	Scalar hydrdt;
	Scalar newdt;     
    Scalar birrdt; 
	Scalar ferrdt; 
	Scalar prordt; 
	Scalar seprdt; 
	Scalar hydrprecdt;     
    Scalar birrprecdt; 
	Scalar ferrprecdt; 
	Scalar prorprecdt;
	Scalar seprprecdt; 	
	Scalar seprdissdt;  	
	Scalar hydrdissdt;     
    Scalar birrdissdt; 
	Scalar ferrdissdt; 
	Scalar prordissdt;	
	
    bool newtonOrBisection_;

    static constexpr Scalar KpHb_ = 0;//9.14e-8;//[mol/kgH2O] Kim et al. 2000 //Not implemented by Anozie!!


        Scalar glAsw_;
        Scalar glrc_;	
        Scalar glbeta_;
        Scalar glsigma_;
        Scalar glp_;


	     Scalar ferAsw_;
         Scalar ferrc_;	
         Scalar ferbeta_;
         Scalar fersigma_;
         Scalar ferp_;

	     Scalar proAsw_;
         Scalar prorc_;	
         Scalar probeta_;
         Scalar prosigma_;
         Scalar prop_;

	     Scalar birAsw_;
         Scalar birrc_;	
         Scalar birbeta_;
         Scalar birsigma_;
         Scalar birp_;
		 
	     Scalar hydAsw_;
         Scalar hydrc_;	
         Scalar hydbeta_;
         Scalar hydsigma_;
         Scalar hydp_;

	     Scalar sepAsw_;
         Scalar seprc_;	
         Scalar sepbeta_;
         Scalar sepsigma_;
         Scalar sepp_;
		 
         Scalar pKaFactor_;

public:


       Scalar glAsw()    {       return glAsw_; }
       Scalar glrc()    {       return glrc_; }
       Scalar glp()    {       return glp_; }
       Scalar glbeta()    {       return glbeta_; }
       Scalar glsigma()    {       return glsigma_; }	   



public:

  /*!
   * \brief Returns the mole fraction of NaCl \f$\mathrm{[mol \ NaCl / mol \ solution]}\f$  for a given mole fraction
   *
   * \param salinity the salinity \f$\mathrm{[kg \ NaCl / kg \ solution]}\f$
   */
  static Scalar salinityToMolFrac_(Scalar salinity) {

        const Scalar Mw = H2O::molarMass(); /* molecular weight of water [kg/mol] */
        const Scalar Ms = 58.8e-3; /* molecular weight of NaCl  [kg/mol] */

        const Scalar X_NaCl = salinity;
        /* salinity: conversion from mass fraction to mol fraction */
        const Scalar x_NaCl = -Mw * X_NaCl / ((Ms - Mw) * X_NaCl - Ms);
        return x_NaCl;
    }
};

} // end namespace

#endif










