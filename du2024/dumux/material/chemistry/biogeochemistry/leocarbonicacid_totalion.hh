
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
    // //what is the enzyme source, heat-killed cells or jack bean?
    // useHeatKilledCells_ = getParam<bool>("Problem.UseHeatKilledCells");
    // useJackBeans_ = getParam<bool>("Problem.UseJackBeans");

    // // Glass parameters
    // glac_        = getParam<Scalar>("GlassCoefficients.glac");
    // glkdiss1_    = getParam<Scalar>("GlassCoefficients.glkdiss1");
    // glkdiss2_    = getParam<Scalar>("GlassCoefficients.glkdiss2");
    // glkprec_     = getParam<Scalar>("GlassCoefficients.glkprec");
    // glndiss_     = getParam<Scalar>("GlassCoefficients.glndiss");
    // glnprec_     = getParam<Scalar>("GlassCoefficients.glnprec");
    // glAsw0_      = getParam<Scalar>("GlassCoefficients.glAsw0");

    // Glass parameters
    // glac_        = getParam<Scalar>("GlassCoefficients.glac");
    // glkdiss1_    = getParam<Scalar>("GlassCoefficients.glkdiss1");
    // glkdiss2_    = getParam<Scalar>("GlassCoefficients.glkdiss2");
    // glkprec_     = getParam<Scalar>("GlassCoefficients.glkprec");
    // glndiss_     = getParam<Scalar>("GlassCoefficients.glndiss");
    // glnprec_     = getParam<Scalar>("GlassCoefficients.glnprec");
    glAsw_      = getParam<Scalar>("GlassCoefficients.glAsw");
    glrc_      = getParam<Scalar>("GlassCoefficients.glrc");	
    glbeta_      = getParam<Scalar>("GlassCoefficients.glbeta");
    glsigma_      = getParam<Scalar>("GlassCoefficients.glsigma");
    glp_      = getParam<Scalar>("GlassCoefficients.glp");	
    // // calcite parameters
    // ac_        = getParam<Scalar>("CalciteCoefficients.ac");
    // kdiss1_    = getParam<Scalar>("CalciteCoefficients.kdiss1");
    // kdiss2_    = getParam<Scalar>("CalciteCoefficients.kdiss2");
    // kprec_     = getParam<Scalar>("CalciteCoefficients.kprec");
    // ndiss_     = getParam<Scalar>("CalciteCoefficients.ndiss");
    // nprec_     = getParam<Scalar>("CalciteCoefficients.nprec");
    // Asw0_      = getParam<Scalar>("CalciteCoefficients.Asw0");
    // calAsw_      = getParam<Scalar>("GlassCoefficients.calAsw");
    // calrc_      = getParam<Scalar>("GlassCoefficients.calrc");	
    // calbeta_      = getParam<Scalar>("GlassCoefficients.calbeta");
    // calsigma_      = getParam<Scalar>("GlassCoefficients.calsigma");
    // calp_      = getParam<Scalar>("GlassCoefficients.calp");	

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
    // // ferrohydrite parameters
    // fac_        = getParam<Scalar>("FerrohydriteCoefficients.fac");
    // fkdiss1_    = getParam<Scalar>("FerrohydriteCoefficients.fkdiss1");
    // fkdiss2_    = getParam<Scalar>("FerrohydriteCoefficients.fkdiss2");
    // fkprec_     = getParam<Scalar>("FerrohydriteCoefficients.fkprec");
    // fndiss_     = getParam<Scalar>("FerrohydriteCoefficients.fndiss");
    // fnprec_     = getParam<Scalar>("FerrohydriteCoefficients.fnprec");
    // fAsw0_      = getParam<Scalar>("FerrohydriteCoefficients.fAsw0");
	
    // // thermal ureolysis parameters
    // cu_        = getParam<Scalar>("UreolysisCoefficients.cu");
    // cuT_       = getParam<Scalar>("UreolysisCoefficients.cuT");
	// 
    // // urease parameters
    // if(useHeatKilledCells_)
    // {
    // //     ureaseInEnzymeSource_ = getParam<Scalar>("UreaseCoefficients.ureaseInEnzymeSource");
    //     cia_       = getParam<Scalar>("UreaseCoefficientsHKC.cia");
    //     ciaT_      = getParam<Scalar>("UreaseCoefficientsHKC.ciaT");
    //     cureaseT_  = getParam<Scalar>("UreaseCoefficientsHKC.cureaseT");
    //     kurease_   = getParam<Scalar>("UreaseCoefficientsHKC.kurease");
    //     ciaPrec_   = getParam<Scalar>("UreaseCoefficientsHKC.ciaPrec");
	// 
    //     //attachment and detachment parameters
    //     ka_urease_ = getParam<Scalar>("UreaseCoefficientsHKC.ka_urease");
    //     kd_urease_ = getParam<Scalar>("UreaseCoefficientsHKC.kd_urease");
    // }
    // if(useJackBeans_)
    // {
    // //     ureaseInEnzymeSource_ = getParam<Scalar>("UreaseCoefficients.ureaseInEnzymeSource");
    //     cia_       = getParam<Scalar>("UreaseCoefficientsJB.cia");
    //     ciaT_      = getParam<Scalar>("UreaseCoefficientsJB.ciaT");
    //     cureaseT_  = getParam<Scalar>("UreaseCoefficientsJB.cureaseT");
    //     kurease_   = getParam<Scalar>("UreaseCoefficientsJB.kurease");
    //     ciaPrec_   = getParam<Scalar>("UreaseCoefficientsJB.ciaPrec");
	// 
    //     //attachment and detachment parameters
    //     ka_urease_ = getParam<Scalar>("UreaseCoefficientsJB.ka_urease");
    //     kd_urease_ = getParam<Scalar>("UreaseCoefficientsJB.kd_urease");
    // }
}

    static const int wPhaseIdx    = FluidSystem::wPhaseIdx;
    static const int nPhaseIdx    = FluidSystem::nPhaseIdx;

    static const int wCompIdx     = FluidSystem::wCompIdx;
    static const int nCompIdx     = FluidSystem::nCompIdx;

    static const int H2OIdx       = FluidSystem::H2OIdx;
    static const int N2Idx        = FluidSystem::N2Idx;
    static const int O2Idx        = FluidSystem::O2Idx;
    // static const int CtotIdx       = FluidSystem::CtotIdx;
	
    static const int CaIdx        = FluidSystem::CaIdx;
    static const int NaIdx        = FluidSystem::NaIdx;
    static const int HIdx        = FluidSystem::HIdx;
    static const int HonlyIdx        = FluidSystem::HonlyIdx;
    // static const int FetotIdx       = FluidSystem::FetotIdx;
    static const int MgIdx        = FluidSystem::MgIdx;
    static const int AlIdx        = FluidSystem::AlIdx;
    static const int SiO2Idx       = FluidSystem::SiO2Idx;
    static const int HPO4Idx        = FluidSystem::HPO4Idx;
    static const int KIdx        = FluidSystem::KIdx;
    static const int MnIdx       = FluidSystem::MnIdx;
    static const int TiOH4Idx       = FluidSystem::TiOH4Idx;
	
    static const int Fe2Idx       = FluidSystem::Fe2Idx;
	// static const int Fe3Idx       = FluidSystem::Fe3Idx;
    static const int ClIdx         = FluidSystem::ClIdx;
    static const int OHIdx        = FluidSystem::OHIdx;
    // // static const int CO2gIdx       = FluidSystem::CO2gIdx;
    static const int CO2aqIdx       = FluidSystem::CO2aqIdx;
    static const int CO2aqonlyIdx       = FluidSystem::CO2aqonlyIdx;
    static const int HCO3Idx      = FluidSystem::HCO3Idx;
    static const int CO3Idx       = FluidSystem::CO3Idx;
    static const int Fe2onlyIdx       = FluidSystem::Fe2onlyIdx;
    static const int Fe3Idx       = FluidSystem::Fe3Idx;
	
    // static const int UreaIdx      = FluidSystem::UreaIdx;
    // static const int UreaseIdx    = FluidSystem::UreaseIdx;
	// 
    // static const int TNHIdx       = FluidSystem::TNHIdx;
    // static const int NH4Idx       = FluidSystem::NH4Idx;

    static const int numComponents      = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numMajorComponents;
    static const int numSecComponents   = FluidSystem::numSecComponents;
    static const int numTotComponents   = numComponents + numSecComponents;
    static const int numPhases          = FluidSystem::numPhases;

    // static const int cPhaseIdx          = SolidSystem::CalciteIdx;
    // static const int uPhaseIdx          = SolidSystem::JbmeIdx;
    static const int fPhaseIdx          = SolidSystem::FerrohydriteIdx;	
	static const int gPhaseIdx          = SolidSystem::GlassIdx;
	static const int pPhaseIdx          = SolidSystem::ProtoImogoliteIdx;
	static const int hPhaseIdx          = SolidSystem::HydroxyapatiteIdx;
	static const int bPhaseIdx          = SolidSystem::BirnessiteIdx;
	static const int sPhaseIdx          = SolidSystem::SepioliteIdx;
    static const int numSolidComponents = SolidSystem::numComponents;
    static const int numInertComponents = SolidSystem::numInertComponents;

    static const int phiGlassIdx      = numComponents + gPhaseIdx;
	// static const int phiCalciteIdx      = numComponents + cPhaseIdx;
    // static const int phiImmUreaseIdx    = numComponents + uPhaseIdx;
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
    // CompVector& calculateEquilibriumChemistry(const FluidState &fluidState, int phaseState, CompVector &variable, Scalar rhoMolar)
 {
        const VolumeVariables volVars{};
		// printf("The value of variables1 is: %.e\n", variable[2]);
        gammaCO2_ = 1.0;
        // h2o_ = 55.508; //molH2O/kgH2O
        Scalar h2o_ = moleFracToMolarity(variable[H2OIdx], rhoMolar);
        pressure_ = fluidState.pressure(wPhaseIdx);
        temperature_ = fluidState.temperature();

        // Scalar moleFracSalinity = variable[NaIdx] + variable[ClIdx] + variable[CaIdx];
        // Scalar moleFracWater = 1-variable[nCompIdx]-variable[CtotIdx]-variable[CaIdx]-variable[NaIdx]-variable[ClIdx]-variable[Fe2Idx];//
		Scalar moleFracWater = variable[H2OIdx]; //
	    // printf("The value of moleFracWater is: %.e\n", moleFracWater);
		// Scalar n2_ = moleFracToMolality(variable[N2Idx], moleFracWater);
		Scalar n2_ = moleFracToMolarity(variable[N2Idx], rhoMolar);
		// Scalar o2_ = moleFracToMolarity(variable[O2Idx], rhoMolar);
        // printf("The value of variables2 is: %.e\n", variable[2]);		
        for (int i = 0; i < numComponents + numSecComponents; ++i)
            {
                // printf("The value of variables2 is: %.10f\n", variable[i]);
				// if(variable[i]<0)
                // {   // printf("The value of variable[i] is: %.e\n", variable[i]); //added by du
                //     // printf("The value of i is: %.d\n", i); //added by du
				//     variable[i]=0;
				// }
				
				//added by du
				// std::cout << "Checking variable[" << i << "]: " << variable[i] << std::endl;
				
				
                if(std::isnan(variable[i]))
                {
                    DUNE_THROW(Dune::InvalidStateException, "Invalid component mole fraction " << variable); break;
                }
            }
		// printf("The value of variables3 is: %.4e\n", variable[2]);
        // printf("The value of co2tot_ is: %.10f\n", variable[2]);
        if (phaseState == bothPhases) //both Phases: solve an open system with co2 concentration constant
        {
            // printf("reached bothphase");// reached 
			// salinity_ = moleFracToMolality(moleFracSalinity, moleFracSalinity, variable[nCompIdx]);

            // co2_ =  (salinity_ + h2o_)/(1/variable[nCompIdx] - 1);
            // ca_ = moleFracToMolality(variable[CaIdx], moleFracSalinity, variable[nCompIdx]);
            // na_ = moleFracToMolality(variable[NaIdx], moleFracSalinity, variable[nCompIdx]);
            // cl_ = moleFracToMolality(variable[ClIdx], moleFracSalinity, variable[nCompIdx]);
            // fe2_ = moleFracToMolality(variable[Fe2Idx], moleFracSalinity, variable[nCompIdx]);
			
            // h2o_ = moleFracToMolality(variable[H2OIdx], moleFracWater);

			// co2aq_ = moleFracToMolality(variable[CO2aqIdx], moleFracWater);
			// co2g_ = moleFracToMolality(variable[CO2gIdx], moleFracWater);
			// co2aq_ = moleFracToMolality(variable[CtotIdx], moleFracWater);
			// hco3_ = moleFracToMolality(variable[CtotIdx], moleFracWater);
			// co2aq_ = moleFracToMolality(variable[CtotIdx], moleFracWater)-moleFracToMolality(variable[CO3Idx], moleFracWater)-moleFracToMolality(variable[HCO3Idx], moleFracWater);		
			// co2aqonly_ = moleFracToMolarity(variable[CO2aqIdx],rhoMolar)-moleFracToMolarity(variable[CO3Idx],rhoMolar)-moleFracToMolarity(variable[HCO3Idx],rhoMolar);		
			co2aq_ = moleFracToMolarity(variable[CO2aqIdx],rhoMolar);//-moleFracToMolarity(variable[CO3Idx],rhoMolar)-moleFracToMolarity(variable[HCO3Idx],rhoMolar);		

			// co2aq_ = moleFracToMolarity(variable[CtotIdx],rhoMolar)-moleFracToMolarity(variable[HCO3Idx],rhoMolar);		
			// ctot_ = moleFracToMolarity(variable[CtotIdx],rhoMolar);		
		    // printf("The value of co2aq_ is: %.10e\n", co2aq_);
			
            o2_ = moleFracToMolarity(variable[O2Idx],rhoMolar);			
			// co2aq_ = ctot_ / (1+1/kh_); // moleFracToMolality(variable[CtotIdx], moleFracWater);; // moleFracToMolality(variable[CO2aqIdx], moleFracWater); //issue is here
			// printf("The value of ctot_ is: %.10f\n", variable[CtotIdx]);
			// printf("The rhoMolar is: %.14f\n", rhoMolar);
			// printf("The value of ctot_ is: %.14f\n", moleFracToMolarity(variable[CtotIdx],rhoMolar));
			// printf("The value of co2aq_ is: %.14f\n", co2aq_); //NAN
            // ca_ = moleFracToMolality(variable[CaIdx], moleFracWater);
			ca_ = moleFracToMolarity(variable[CaIdx],rhoMolar);
			// printf("The value of ca_ is: %.8lf\n", ca_);
            // na_ = moleFracToMolality(variable[NaIdx], moleFracWater);
            cl_ = moleFracToMolality(variable[ClIdx], rhoMolar);//moleFracWater);
            fe2_ = moleFracToMolality(variable[Fe2Idx], rhoMolar);
			// printf("The value of variable[Fe2Idx] is: %.8e\n", fe2_);
			na_ = moleFracToMolarity(variable[NaIdx],rhoMolar);
			//[H+] – [OH-] –[HCO3]– 2[CO32-] 
            // honly_ = moleFracToMolarity(variable[HIdx],rhoMolar)+moleFracToMolarity(variable[OHIdx],rhoMolar)+moleFracToMolarity(variable[CO3Idx],rhoMolar)+2*moleFracToMolarity(variable[HCO3Idx],rhoMolar);
            h_ = moleFracToMolarity(variable[HIdx],rhoMolar);//+moleFracToMolarity(variable[OHIdx],rhoMolar)+moleFracToMolarity(variable[CO3Idx],rhoMolar)+2*moleFracToMolarity(variable[HCO3Idx],rhoMolar);
			// printf("The value of h_ is: %.e\n", h_);
			
            // fetot_ = moleFracToMolarity(variable[FetotIdx],rhoMolar);
            mg_ = moleFracToMolarity(variable[MgIdx],rhoMolar);
            al_ = moleFracToMolarity(variable[AlIdx],rhoMolar);
            sio2_ = moleFracToMolarity(variable[SiO2Idx],rhoMolar);
            hpo4_ = moleFracToMolarity(variable[HPO4Idx],rhoMolar);
            k_ = moleFracToMolarity(variable[KIdx],rhoMolar);
            mn_ = moleFracToMolarity(variable[MnIdx],rhoMolar);
			// printf("The value of variable[MnIdx] is: %.8e\n", variable[MnIdx]);
			// cl_ = moleFracToMolarity(variable[ClIdx],rhoMolar);
            tioh4_ = moleFracToMolarity(variable[TiOH4Idx],rhoMolar);
			// all zero here!!
			// n2_ = moleFracToMolality(variable[N2Idx], moleFracWater);			// fe2_ = 0; // added by du
            // totalnh_ = moleFracToMolality(variable[TNHIdx], moleFracSalinity, variable[nCompIdx]);
			// totalnh_ = 0; // added by du
           // Scalar m = na_ + ca_;
           // Scalar Temp = fluidState.temperature();
           // /* Millero et al. 2007: The dissociation of carbonic acid */
           // /* in NaCl solutions as a function of concentration and temperature */
		   // 
           // /*for pK1*/
           // Scalar a0 = 31.3616;Scalar  a1 = 0.86644;Scalar  a2 = -0.33611;Scalar  a3 = 0.05888;
           // Scalar  b0 = -1422.25317; Scalar c0 = -4.84141;
		   // 
//         //    Scalar A = a0*sqrt(m) + a1*m+ a2*pow(m,1.5) + a3*m*m ;
           // Scalar A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
//         //    Scalar B = b0 * pow (m,0.5); Scalar C = c0*pow(m,0.5);
           // Scalar B = b0 * sqrt(m);
           // Scalar C = c0*sqrt(m);
		   // 
           // Scalar dpK1 = A +  B/Temp + C*log(Temp);
           // Scalar pK1 = - 402.56788 + 11656.46/Temp + 72.173*log(Temp) - 0.161325*Temp + 7.5526E-5*Temp*Temp;
		   // 
           // pK1 = pK1 + dpK1;
		   // 
           // /*for pK2*/
           // a0 = 36.88545; a1 = 1.66599; a2 = -0.68730; a3 = 0.12070;
           // b0 = -1669.55918; c0 = -5.83555;
		   // 
           // A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
           // B = b0 * sqrt(m);
           // C = c0*sqrt(m);
		   // 
           // Scalar dpK2 = A +  B/Temp + C*log(Temp);
           // Scalar  pK2 = -122.4994 + 5811.18/Temp + 20.5263*log(Temp) - 0.0120897*Temp;
		   // 
           // pK2 = pK2 + dpK2;
		   // 
           // /*Bell et al. 2008: Ammonia/ammonium dissociation coefficient*/
           // /*in seawater: A significant numerical correction*/
           // Scalar I_f = 0.5*(na_ + 4.*ca_ + cl_);                         /*ionic strength of salt solution: here, equivalent to m, neglecting other ions*/
//         //    Scalar pKa = 10.0423 - 0.0315536*(Temp-273.15) + 0.14737*I_f;
//         //    Olofsson 1975: Thermodynamic quantities for the dissociation of the ammonium Ion and for the ionization of aqueous ammonia over a wide temperature range (up to 570K!)
           // Scalar pKa = 2533/Temp - 0.5936* log(Temp) + 4.127;
		   
           // //using the ionic strength correction by Bell et al. 2008
           // pKa += 0.14737*I_f;
           // pKa *= pKaFactor_;
		   // 
           // apparentk1_ = pow (10.,-pK1);
           // apparentk2_ = pow (10.,-pK2);
           // apparentka_ = pow (10.,-pKa);
		   // 
           // Du deleted
		   // initH_ = 1e-5; //Initial guess
           // Scalar activityH = initH_;
	       // k1_ = const1(pressure_, temperature_);
           // k2_ = const2(pressure_, temperature_);
           // k3_ = const3(pressure_, temperature_);
           // k4_ = const4(pressure_, temperature_);
           // // kh_ = constco2(pressure_, temperature_);   
           // kw_ = constW(pressure_, temperature_);
		   // // 
           // Scalar tolAbs = 1e-20;
           // Scalar tolRel = 1e-15;
           // int maxIter = 30;
		   // 
           // //Do the Newton iteration and calculate the components molalities and update the mass fraction array and
           // //the concentration of the changed primary variables
           // // if(newton1D(activityH, &ThisType::H_CO2, tolAbs, tolRel, maxIter) == false) //Alex' Newton
           // if(newton1D(activityH, &ThisType::H_CO2, tolAbs, tolRel, maxIter) == false) //Alex' Newton
		   // 
//         //   if(newton1D(activityH, tolAbs, maxIter) == false) //Anozies Newton
           // {
			//// printf("reached nonnewton"); 
           //    initH_ = 1e-5;
			//   // initCO2aq_ = 0.0132074;
			//   // printf("The value of initH2 is: %.8lf\n", initH_);
           //     activityH = initH_;
           //     // initCl_ = 1e-5;
           //     // activityCl = initCl_;				
           //     Scalar a0 = 0.0;
           //     Scalar b0 = 1e-1;
           //     Scalar tol = 1e-15;
           //     // if(bisection1D(activityH, &ThisType::H_CO2, a0, b0, tol) == false) //Alex' bisection
           //     if(bisection1D(activityH, &ThisType::H_CO2, a0, b0, tol) == false) //Alex' bisection
		   // 
//         //       if(bisection1D(tol) == false) //Anozies bisection
           //     {
           //         DUNE_THROW(Dune::InvalidStateException, "in Chemistry: Bisection did not converge!" );
           //     }
           // }
           // H_CO2(activityH); //update component molalities
           // h2_ = h_;
		   // 
		   // printf("The value of co2aq_ is: %.8e\n", co2aq_);
			//// printf("The value of hco3_ is: %.8e\n", hco3_);
			//// printf("The value of co3_ is: %.8e\n", co3_);
		   // printf("The value of h_ is: %.8e\n", h_);		   
			//
           // // h_ = moleFracToMolarity(variable[HIdx],rhoMolar);		   
		   // // 
		   // initCO2aq_ = 0.0132074; //Initial guess
		   // printf("The value of initH is: %.8lf\n", initH_);
           // Scalar activityCO2aq = initCO2aq_;
		   initH_ = 8.07096e-13;//h_; //Initial guess
		   // // printf("The value of initH is: %.8lf\n", initH_);
           Scalar activityH = initH_;	
		   // initCO3_ = co2aq_; //Initial guess
		   // printf("The value of initH is: %.8lf\n", initCO3_);
           // Scalar activityCO3 = initCO3_;		   
           // //Anozies apparent constants
           // k1_ = apparentk1_;
           // k2_ = apparentk2_;
           // ka_ = apparentka_;
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
		   // 
           Scalar tolAbs = 1e-20; //1e-20;
           Scalar tolRel = 1e-20;// 1e-15;
           int maxIter = 40;
		   
		   // printf("The value of co2aq1 is: %.8e\n", co2aq_);
		   // printf("The value of co31 is: %.8e\n", co3_);
		   // printf("The value of hco31 is: %.8e\n", hco3_);
		   // printf("The value of oh1 is: %.8e\n", oh_);
		   // printf("The value of h1 is: %.8e\n", h_);
		   // h_ = h_ + oh_ + 2*co3_ + hco3_;
		   // co2aq_ = co2aq_ - co3_ - hco3_;
		   // printf("The value of co2aq2 is: %.8e\n", co2aq_);
		   // printf("The value of co32 is: %.8e\n", co3_);
		   // printf("The value of hco32 is: %.8e\n", hco3_);
		   // printf("The value of oh2 is: %.8e\n", oh_);
		   // printf("The value of h2 is: %.8e\n", h_);
           // oh_ = kw_/honly_;
		   // //  // printf("The value of oh_ is: %.8e\n", oh_);
		   // //  // co2aq_ = ctot_;
		   // //  // co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + (1/kh_) + 1);
		   // //  // co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + 1);
           // hco3_ = k1_*co2aqonly_/honly_;
		   // //  // printf("The value of co2aq is: %.8e\n", co2aq_);
		   // //  // printf("The value of k1_ is: %.e\n", k1_);
//         // //    co3_ = k1_*k2_*co2_/pow(h_, 2.);
           // co3_ = k1_*k2_*co2aqonly_/(honly_*honly_);
		   // // 
		   // printf("The value of co2aq3 is: %.8e\n", co2aq_);
		   // printf("The value of co33 is: %.8e\n", co3_);
		   // printf("The value of hco33 is: %.8e\n", hco3_);
		   // printf("The value of oh3 is: %.8e\n", oh_);
		   // printf("The value of h3 is: %.8e\n", h_);
		   // h_ = honly_ - oh_ - 2*co3_ - hco3_;
		   // co2aq_ = co2aqonly_ + co3_ + hco3_;
		   // printf("The value of co2aq4 is: %.8e\n", co2aq_);
		   // printf("The value of co34 is: %.8e\n", co3_);
		   // printf("The value of hco34 is: %.8e\n", hco3_);
		   // printf("The value of oh4 is: %.8e\n", oh_);
		   // printf("The value of h4 is: %.8e\n", h_);		 
           //Do the Newton iteration and calculate the components molalities and update the mass fraction array and
           //the concentration of the changed primary variables
           // if(newton1D(activityH, &ThisType::H_CO2, tolAbs, tolRel, maxIter) == false) //Alex' Newton
           // if(newton1D(activityCl, &ThisType::Cl_CO2, tolAbs, tolRel, maxIter) == false) //Alex' Newton
		    // if(newton1D(activityCO2aq, &ThisType::H_varyCO2, tolAbs, tolRel, maxIter) == false)
		    // if(newton1D(activityCO3, &ThisType::H_CO3, tolAbs, tolRel, maxIter) == false)
			if(newton1D(activityH, &ThisType::H_Ctot, tolAbs, tolRel, maxIter) == false) //Alex' Newton
//            if(newton1D(activityH, tolAbs, maxIter) == false) //Anozies Newton
            {
		        // printf("reached nonnewton"); 
                initH_ = 8.07096e-13 ;//h_;
		        // initCl_ = cl_;
		        // initCO2aq_ = co2aq_;
				// initCO3_ = co2aq_;
		        // printf("The value of initH2 is: %.8lf\n", initH_);
                // activityCO2aq = initCO2aq_;
				activityH = initH_;
				// activityCO3 = initCO3_;
                // initCl_ = 1e-5;
                // activityCl = initCl_;				
                Scalar a0 = 0.0;
                Scalar b0 = 1e-1;
                Scalar tol = 1e-15;//1e-15;
                // if(bisection1D(activityH, &ThisType::H_CO2, a0, b0, tol) == false) //Alex' bisection
				if(bisection1D(activityH, &ThisType::H_Ctot, a0, b0, tol) == false) //Alex' bisection
                // if(bisection1D(activityCO2aq, &ThisType::H_varyCO2, a0, b0, tol) == false) //Alex' bisection
                // if(bisection1D(activityCO3, &ThisType::H_CO3, a0, b0, tol) == false) //Alex' bisection
		    
//                if(bisection1D(tol) == false) //Anozies bisection
                {
                    DUNE_THROW(Dune::InvalidStateException, "in Chemistry: Bisection did not converge!" );
                }
            }
            // // H_varyCO2(activityCO2aq); //update component molalities
			// H_CO3(activityCO3);
			// oh_ = kw_/h_;
            // hco3_ = k1_*co2aq_/h_;
            // co3_ = k1_*k2_*co2aq_/(h_*h_);
			// printf("The value of fe21_ is: %.8e\n", fe2_);
		    H_Ctot(activityH);
		    // co2aq_CO2(activityco2aq); //update component molalities
            // printf("The value of co2g_ is: %.e\n", co2g_);
			// printf("The value of co2aqonly_ is: %.8e\n", co2aqonly_);
			// printf("The value of honly_ is: %.8e\n", honly_);
			// printf("The value of co2aq_ is: %.8e\n", co2aq_);
			// printf("The value of hco3_ is: %.8e\n", hco3_);
			// printf("The value of co3_ is: %.8e\n", co3_);
			// printf("The value of h_ is: %.8e\n", h_);
			// printf("The value of oh_ is: %.8e\n", oh_);
            //update mole fractions in the variable vector for the open system

            // co2aq_ = hco3_ + co3_ + co2aq_; //calculate the molarity of cTot from the other c-components
		    h_ = honly_ - oh_ - 2*co3_ - hco3_ + fe3_ ;
		    co2aq_ = co2aqonly_ + co3_ + hco3_;
			fe2_ = fe2only_ + fe3_ ; 
            // printf("The value of honly_ is: %.8e\n", honly_);
			// printf("The value of o2_ is: %.8e\n", o2_);
			// printf("The value of co2aq_ is: %.e\n", co2aq_);
			// printf("The value of hco3_ is: %.e\n", hco3_);
			// printf("The value of co3_ is: %.e\n", co3_);  
            // printf("The value of fe2_ is: %.8e\n", fe2_);	
            // printf("The value of fe2only_ is: %.8e\n", fe2only_);			
            // printf("The value of fe3_ is: %.8e\n", fe3_);						
			// cTot_ = co2_ ; // added by du
			
            // Scalar moleFracCTot = molalityToMoleFrac(ctot_, moleFracWater);
            // Scalar moleFracCTot = molarityToMoleFrac(ctot_,rhoMolar);  
//             Scalar urea = moleFracToMolality(variable[UreaIdx], moleFracWater);
            // Scalar urea = moleFracToMolality(variable[UreaIdx], moleFracSalinity, variable[nCompIdx]);
            // Scalar urease = moleFracToMolality(variable[UreaseIdx], moleFracSalinity, variable[nCompIdx]);

            // Scalar totalMolality = h2o_ + cTot_ + na_ + cl_ + ca_ + totalnh_ + urea + urease + fe2_;
            // Scalar totalMolarity = h2o_ + n2_ +  h_ + o2_ + ctot_  + na_ + cl_ + ca_ + fe2_ + al_ + mg_ +k_ + mn_ + hpo4_ + tioh4_ + sio2_; // du do not know
            Scalar totalMolarity = h2o_ + n2_ + o2_ +  co2aq_  + na_ + cl_ + ca_ + fe2_ + al_ + mg_ +k_ + mn_ + hpo4_ + tioh4_ + sio2_ + h_; // du do not know
            // Scalar totalMolarity = h2o_ + n2_ + o2_ + na_ + cl_ + ca_ + fe2_ + al_ + mg_ +k_ + mn_ + hpo4_ + tioh4_ + sio2_ + h_; // du do not know

			// // printf("The value of totalMolality] is: %.8lf\n", totalMolality); //Nan
            // variable[CtotIdx] = ctot_/totalMolarity; //calculate the mole fraction of cTot in terms of mol CO2 / mol solution
            // // variable[CO2gIdx] = co2g_/totalMolality;
			variable[CO2aqIdx] = co2aq_/totalMolarity;
			variable[CO2aqonlyIdx] = co2aqonly_/totalMolarity;
            variable[HCO3Idx] = hco3_/totalMolarity;
            variable[CO3Idx] = co3_/totalMolarity;
            // // variable[NH4Idx] = nh4_/totalMolality;
            variable[OHIdx] = oh_/totalMolarity;
            variable[HIdx] = h_/totalMolarity;
            variable[HonlyIdx] = honly_/totalMolarity;
			variable[Fe2onlyIdx] = fe2only_/totalMolarity;
			variable[Fe2Idx] = fe2_/totalMolarity;
			variable[Fe3Idx] = fe3_/totalMolarity;
            // variable[Fe3Idx] = fe3_/totalMolarity;
		    // printf("The value of moleFracforchem0 is: %.5e\n", variable[CO2aqIdx]);	
		    // printf("The value of totalMolarity is: %.5e\n", totalMolarity);			
			// same as the output mole fraction
			// printf("The value of variable[CtotIdx] is: %.8lf\n", variable[CtotIdx]);
			// printf("The value of variable[CO2aqonlyIdx] is: %.8e\n", variable[CO2aqonlyIdx]);
			// printf("The value of variable[CO2aqIdx] is: %.8e\n", variable[CO2aqIdx]);
			// printf("The value of variable[HCO3Idx] is: %.8lf\n", variable[HCO3Idx]);
			// printf("The value of variable[CO3Idx] is: %.8lf\n", variable[CO3Idx]);
			// printf("The value of variable[H2Idx] is: %.8e\n", variable[HIdx]);
			// printf("The value of variable[OHIdx] is: %.8lf\n", variable[OHIdx]); // Nan
     
			// Scalar pH = -log10(h_); // added by du
			// printf("The value of pH is: %.e\n", pH); 
        }

        else if (phaseState == wPhaseOnly) //wPhaseOnly: solve a closed system with cTot concentration constant
        {
			// printf("reached wphaseonly"); 
            // salinity_ = moleFracToMolality(moleFracSalinity,  moleFracSalinity, variable[nCompIdx]);
            // n2_ =  moleFracToMolality(variable[nCompIdx], moleFracSalinity, variable[nCompIdx]);
            // ca_ = moleFracToMolality(variable[CaIdx], moleFracSalinity, variable[nCompIdx]);
            // na_ = moleFracToMolality(variable[NaIdx], moleFracSalinity, variable[nCompIdx]);
            // cl_ = moleFracToMolality(variable[ClIdx],  moleFracSalinity, variable[nCompIdx]);
            // fe2_ = moleFracToMolality(variable[Fe2Idx], moleFracSalinity, variable[nCompIdx]);  
            // co2_ = moleFracToMolality(variable[CO2Idx], moleFracSalinity, variable[nCompIdx]);
			// co2aq_ = moleFracToMolality(variable[CO2aqIdx], moleFracWater);
            // ca_ = moleFracToMolality(variable[CaIdx], moleFracWater);
            // na_ = moleFracToMolality(variable[NaIdx], moleFracWater);
            // cl_ = moleFracToMolality(variable[ClIdx], moleFracWater);
            // fe2_ = moleFracToMolality(variable[Fe2Idx], moleFracWater);		
			// co2aq_ = moleFracToMolarity(variable[CtotIdx],rhoMolar)-moleFracToMolarity(variable[CO3Idx],rhoMolar)-moleFracToMolarity(variable[HCO3Idx],rhoMolar);		
			co2aq_ = moleFracToMolarity(variable[CO2aqIdx],rhoMolar);				
			ca_ = moleFracToMolarity(variable[CaIdx],rhoMolar);
			na_ = moleFracToMolarity(variable[NaIdx],rhoMolar);
            cl_ = moleFracToMolarity(variable[ClIdx],rhoMolar);
            fe2_ = moleFracToMolarity(variable[Fe2Idx],rhoMolar);			
            // fe2_ = 0;   // added by du
		    // totalnh_ = moleFracToMolality(variable[TNHIdx], moleFracSalinity, variable[nCompIdx]);
            // totalnh_ = 0;   // added by du
            // Scalar m = na_ + ca_;
            // Scalar Temp = fluidState.temperature();
            // /* Millero et al. 2007: The dissociation of carbonic acid */
            // /* in NaCl solutions as a function of concentration and temperature */
			// 
            // /*for pK1*/
            // Scalar a0 = 31.3616;Scalar  a1 = 0.86644;Scalar  a2 = -0.33611;Scalar  a3 = 0.05888;
            // Scalar  b0 = -1422.25317; Scalar c0 = -4.84141;
			// 
            // Scalar A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
            // Scalar B = b0 * sqrt(m);
            // Scalar C = c0*sqrt(m);
			// 
            // Scalar dpK1 = A +  B/Temp + C*log(Temp);
            // Scalar pK1 = - 402.56788 + 11656.46/Temp + 72.173*log(Temp) - 0.161325*Temp + 7.5526E-5*Temp*Temp;
			// 
            // pK1 = pK1 + dpK1;
			// 
            // /*for pK2*/
            // a0 = 36.88545; a1 = 1.66599; a2 = -0.68730; a3 = 0.12070;
            // b0 = -1669.55918; c0 = -5.83555;
			// 
//          //    A = a0*pow(m,0.5) + a1*m+ a2*pow(m,1.5) + a3*m*m;
            // A = a0*sqrt(m) + a1*m+ a2*sqrt(m*m*m) + a3*m*m ;
//          //    B = b0 * pow (m,0.5); C = c0*pow(m,0.5);
            // B = b0 * sqrt(m);
            // C = c0*sqrt(m);
			// 
            // Scalar dpK2 = A +  B/Temp + C*log(Temp);
            // Scalar  pK2 = -122.4994 + 5811.18/Temp + 20.5263*log(Temp) - 0.0120897*Temp;
			// 
            // pK2 = pK2 + dpK2;
			// 
            // /*Bell et al. 2008: Ammonia/ammonium dissociation coefficient*/
            // /*in seawater: A significant numerical correction*/
            // Scalar I_f = 0.5*(na_ + 4.*ca_ + cl_);                         /*ionic strength of salt solution: here, equivalent to m, neglecting other ions*/
//          //    Scalar pKa = 10.0423 - 0.0315536*(Temp-273.15) + 0.14737*I_f;
//          //    Olofsson 1975: Thermodynamic quantities for the dissociation of the ammonium Ion and for the ionization of aqueous ammonia over a wide temperature range (up to 570K!)
            // Scalar pKa = 2533/Temp - 0.5936* log(Temp) + 4.127;
            // //using the ionic strength correction by Bell et al. 2008
            // pKa += 0.14737*I_f;
            // pKa *= pKaFactor_;
			// 
            // apparentk1_ = pow (10.,-pK1);
            // apparentk2_ = pow (10.,-pK2);
            // apparentka_ = pow (10.,-pKa);
			// 
            // //Anozies apparent constants
            // k1_ = apparentk1_;
            // k2_ = apparentk2_;
            // ka_ = apparentka_;
            kw_ = constW(pressure_, temperature_);

            //Parameters for the newton solver
//            Scalar CTotLow = 1e-4;
//            Scalar CTotHigh = 1e-2;
            Scalar tolAbs = 1e-20; //1e-11;old
            Scalar tolRel = 1e-20; //1e-11;old
            int maxIter = 40;
            initH_ = 1e-7;
            Scalar activityH = initH_;
            // if(newton1D(activityH, &ThisType::H_Ctot, tolAbs, tolRel, maxIter) == true) //Alex' Newton
//          // if(newton1D(activityH, tolAbs, maxIter) == true)  //Anozies Newton
            // {
            //     //update all component molalities
            //     H_Ctot(activityH);
            // }
            // else //solve with the bisection method and hco3 as primary variable
            // {
            //     Scalar a0 = 0;
            //     Scalar b0 = 1e-3;//1e-4; old
//          //       initH_ = 1e-5;
            //     Scalar tol = 1e-12;//1e-8; old
            //     Scalar activityH = b0;
            //     if(bisection1D(activityH, &ThisType::H_Ctot, a0, b0, tol) == true) //Alex' Bisection
//          //       if(bisection1D(tol) == true) //Anozies Bisection
            //     {
            //         H_Ctot(activityH);//   CTot_HCO3(activityHCO3); //update all component molalities
            //     }
            //     else
            //     {
            //         DUNE_THROW(Dune::InvalidStateException, "in Chemistry: Bisection did not converge!" );
            //     }
			// 
            // }
			// 
            // // Scalar urea = moleFracToMolality(variable[UreaIdx], moleFracSalinity, variable[nCompIdx]);
            // // Scalar urease = moleFracToMolality(variable[UreaseIdx], moleFracSalinity, variable[nCompIdx]);
			// 
            // // Scalar totalMolality = h2o_ + cTot_ + na_ + cl_ + ca_ + totalnh_ + urea + urease + fe2_;
            // Scalar totalMolarity = h2o_ + ctot_ + n2_ + o2_ + na_ + cl_ + ca_ + fe2_;
            // // calculate the secondary component mole fractions
            // // variable[CO2Idx] = co2_/totalMolality;
            // // variable[HCO3Idx] = hco3_/totalMolality;
            // // variable[CO3Idx] = co3_/totalMolality;
            // // variable[NH4Idx] = nh4_/totalMolality; added by du
			// // variable[NH4Idx] = 0;
            // variable[OHIdx] = oh_/totalMolarity;
            // variable[HIdx] = h_/totalMolarity;
			// 
            // // Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_ + 2*fe2_;
            // // Scalar fmolfrac = 2*variable[CaIdx] + variable[NaIdx] + variable[NH4Idx] + variable[HIdx] - variable[ClIdx] - variable[HCO3Idx] - 2*variable[CO3Idx] - variable[OHIdx] + 2*variable[Fe2Idx];
            // Scalar f = na_ + h_ + 2*ca_ - oh_ - cl_  + 2*fe2_ + 3*fe3_;
            // Scalar fmolfrac = 2*variable[CaIdx] + variable[NaIdx] + variable[HIdx] - variable[ClIdx] - variable[OHIdx] + 2*variable[Fe2Idx]  + 3*variable[Fe3Idx];

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
    // //Return equlibrium constant for chemical equation:
    // // NH4 <--> H + NH3
    // /*static*/ Scalar consta(const Scalar pw, const Scalar T)
    // {
    //     return 5.12861e-10;//return(pow(10,-9.29)); //pow(10,-9.25)
    // }


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
        // + mNH4  * FluidSystem::charge(NH4Idx) * FluidSystem::charge(NH4Idx)
		+ mFe2  * FluidSystem::charge(Fe2Idx) * FluidSystem::charge(Fe2Idx));
        return ionicStrength;
    }
    // static Scalar ionicStrength(Scalar mNa, Scalar mCl, Scalar mCa, Scalar mNH4 )
    // {
    //     Scalar ionicStrength = 0.5*( mNa    * FluidSystem::charge(NaIdx) * FluidSystem::charge(NaIdx)
    //     + mCl   * FluidSystem::charge(ClIdx) * FluidSystem::charge(ClIdx)
    //     + mCa   * FluidSystem::charge(CaIdx) * FluidSystem::charge(CaIdx)
    //     // + mNH4  * FluidSystem::charge(NH4Idx) * FluidSystem::charge(NH4Idx)
	// 	);
	// 
    //     return ionicStrength;
    // }
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
      // Scalar beta_nahco3_0, beta_nahco3_1, C_nahco3_phi;
      // Scalar beta_naco3_0, beta_naco3_1, C_naco3_phi;
      // Scalar psi_canacl, psi_co3nacl, theta_naca, theta_clco3;
      // Scalar B1_nacl, C_nacl, B1_nahco3, C_nahco3, B1_naco3, C_naco3, B_naco3;
      // Scalar A_phi, a[6], T,x_clco3,x_clcl, x_co3co3,x_cana,x_caca,x_nana;
      // Scalar E_theta_cana, E_theta_clco3, E1_theta_cana, E1_theta_clco3;
      Scalar psi_canacl, theta_naca;
      Scalar B1_nacl, C_nacl;
      Scalar A_phi, a[6], T,x_clcl, x_cana,x_caca,x_nana;
      Scalar E_theta_cana,  E1_theta_cana;
      // Scalar beta_nh4cl_0, beta_nh4cl_1, beta_nh4co3_0, beta_nh4co3_1, beta_nh4hco3_0, beta_nh4hco3_1; /*new*/
      // Scalar B_nh4cl, B_nh4co3, B_nh4hco3, B1_nh4cl, B1_nh4co3, B1_nh4hco3, C_nh4cl, C_nh4co3, C_nh4cl_phi, C_nh4co3_phi; /*new*/

      // I = 0.5*( mNa + 4.*mCa + mNH4 + mHCO3 + 4*mCO3 + mCl) + 1.E-20;
      I = 0.5*( mNa + 4.*mCa + mCl) + 1.E-20;
      sqrt_I = sqrt(I);

      T = temp;
      a[0]=-8.1765300E-1; a[1]=-8.6852760E-1; a[2]=1.9251000E+4; a[3]=5.2514840E-3; a[4]=-7.1493970E-6; a[5]=9.3385590E-12;

      A_phi = a[0] + a[1]/(T-222.) + a[2]/(T*T) + a[3]*T + a[4]*T*T + a[5]*T*T*T*T;
      /*MODELING AND NUMERICAL SIMULATION OF SALT TRANSPORT AND PHASE TRANSITIONS IN UNSATURATED POROUS BUILDING MATERIALS By Andreas Nicolai*/

      beta_cacl_0 = 0.3159;  beta_cacl_1 = 1.614; C_cacl_phi = -0.00034;
      beta_nacl_0 = 0.0765; beta_nacl_1 = 0.2664; C_nacl_phi = 0.00127;
      // beta_nahco3_0 = 0.0277; beta_nahco3_1 = 0.0411; C_nahco3_phi = 0.0;
      // beta_naco3_0 = 0.1898; beta_naco3_1 = 0.846; C_naco3_phi = -0.048;
      psi_canacl = -0.014; // psi_co3nacl = 0.016;
      theta_naca = 0.07; // theta_clco3 = -0.053;


      // beta_nh4co3_0 = 0.1288; beta_nh4co3_1 = 1.433; C_nh4co3_phi = 0.0005; /*new*/
      // beta_nh4hco3_0 = -0.038; beta_nh4hco3_1 = 0.07;                       /*new*/
      // beta_nh4cl_0 = 0.0522; beta_nh4cl_1 = 0.1918; C_nh4cl_phi = 0.003;    /*new*/


      // x_clco3 = 6.*(-1.)*(-2.)*A_phi*sqrt_I;
      x_clcl = 6.*(-1.)*(-1.)*A_phi*sqrt_I;
      // x_co3co3 = 6.*(-2.)*(-2.)*A_phi*sqrt_I;
      x_cana = 6.*(+2.)*(+1.)*A_phi*sqrt_I;
      x_caca = 6.*(+2.)*(+2.)*A_phi*sqrt_I;
      x_nana = 6.*(+1.)*(+1.)*A_phi*sqrt_I;


      E_theta_cana = ((+2.)*(+1.)/(4.*I))*( J(x_cana) - 0.5*J(x_caca) - 0.5*J(x_nana) );
      // E_theta_clco3 = ((-1.)*(-2.)/(4.*I))*( J(x_clco3) - 0.5*J(x_clcl) - 0.5*J(x_co3co3) );

      E1_theta_cana = -(E_theta_cana/I) + ((+2)*(+1)/(8*I*I))*( x_cana*Jprime(x_cana) - 0.5*x_caca*Jprime(x_caca) - 0.5*x_nana*Jprime(x_nana) );
      // E1_theta_clco3 = -(E_theta_clco3/I) + ((-1)*(-2)/(8*I*I))*( x_clco3*Jprime(x_clco3) - 0.5*x_clcl*Jprime(x_clcl) - 0.5*x_co3co3*Jprime(x_co3co3) );

        f = -A_phi * ( sqrt_I/(1. + 1.2*sqrt_I) + (2./1.2)*log(1. + 1.2*sqrt_I) );
        B_cacl = beta_cacl_0 + (beta_cacl_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));
        B1_cacl = (beta_cacl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_cacl = C_cacl_phi  / (2.*pow ( (2.*1.) , 0.5 ));
        C_cacl = C_cacl_phi  / (2.*sqrt(2.*1.));

        B1_nacl = (beta_nacl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_nacl = C_nacl_phi  / (2.*pow ( (1.*1.) , 0.5 ));
        C_nacl = C_nacl_phi  / (2.*sqrt(1.*1.));

        // B1_nahco3 = (beta_nahco3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_nahco3 = C_nahco3_phi  / (2.*pow ( (1.*1.) , 0.5 ));
        // C_nahco3 = C_nahco3_phi  / (2.*sqrt(1.*1.));

        // B_naco3 = beta_naco3_0 + (beta_naco3_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));
        // B1_naco3 = (beta_naco3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));
//         C_naco3 = C_naco3_phi  / (2.*pow ( (1.*2.) , 0.5 ));
        // C_naco3 = C_naco3_phi  / (2.*sqrt(1.*2.));


        // B_nh4cl = beta_nh4cl_0 + (beta_nh4cl_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));  /*new*/
        // B1_nh4cl = (beta_nh4cl_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));    /*new*/
//      //    C_nh4cl = C_nh4cl_phi  / (2.*pow ( (2.*1.) , 0.5 ));       /*new*/
        // C_nh4cl = C_nh4cl_phi  / (2.*sqrt(2.*1.));                                         /*new*/
		// 
        // B_nh4co3 = beta_nh4co3_0 + (beta_nh4co3_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));  /*new*/
        // B1_nh4co3 = (beta_nh4co3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));      /*new*/
//      //    C_nh4co3 = C_nh4co3_phi  / (2.*pow ( (1.*2.) , 0.5 ));
        // C_nh4co3 = C_nh4co3_phi  / (2.*sqrt(1.*2.));                                          /*new*/
		// 
        // B_nh4hco3 = beta_nh4hco3_0 + (beta_nh4hco3_1 / (2.*I)) * (1. - exp(-2.*sqrt_I) * (1. + 2.*sqrt_I));  /*new*/
        // B1_nh4hco3 = (beta_nh4hco3_1 / (2.*I*I)) * ( -1. + exp(-2.*sqrt_I) * ( 1. + 2.*sqrt_I + 2.*I));      /*new*/



        gamma_Ca = exp (
                4.*f
                + mCl*(2.*B_cacl + mCl*C_cacl)
                + mNa*mCl*(4.*B1_nacl + 2.*C_nacl)
                // + mNa*mHCO3*(4.*B1_nahco3 + 2.*C_nahco3)
                // + mNa*mCO3*(4.*B1_naco3 + 2.*C_naco3)
                + mCa*mCl*(4.*B1_cacl + 2.*C_cacl)

                        // + mNH4*mCl*(4.*B1_nh4cl + 2.*C_nh4cl)        /*new*/
                        // + mNH4*mHCO3*(4.*B1_nh4hco3)  /*new*/
                        // + mNH4*mCO3*(4.*B1_nh4co3 + 2.*C_nh4co3)     /*new*/


                + mNa*(2.*theta_naca + 2.*E_theta_cana + mCl*psi_canacl)
                + 4.*mNa*mCa*E1_theta_cana // + 4.*mCl*mCO3*E1_theta_clco3
                );
        // gamma_CO3 = exp (
        //          4.*f
        //          + mNa*(2.*B_naco3 + mNa*C_naco3)
        //          + mNa*mCl*(4.*B1_nacl + 2.*C_nacl)
        //          + mNa*mHCO3*(4.*B1_nahco3 + 2.*C_nahco3)
        //          + mNa*mCO3*(4.*B1_naco3 + 2.*C_naco3)
		// 
        //                 // + mNH4*(2.*B_nh4co3 + mNH4*C_nh4co3)        /*new*/
        //                 // + mNH4*mCl*(4.*B1_nh4cl + 2.*C_nh4cl)        /*new*/
        //                 // + mNH4*mHCO3*(4.*B1_nh4hco3)                /*new*/
        //                 // + mNH4*mCO3*(4.*B1_nh4co3 + 2.*C_nh4co3)       /*new*/
		// 
        //          + mCa*mCl*(4.*B1_cacl + 2.*C_cacl)
        //          + mCl*(2.*theta_clco3 + 2.*E_theta_clco3 + mNa*psi_co3nacl)
        //          + 4.*mNa*mCa*E1_theta_cana + 4.*mCl*mCO3*E1_theta_clco3
        //          );

//        Ksp = pow (10.,-8.48)/(gamma_Ca*gamma_CO3);
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
   //  Scalar Omega(const Scalar mNa,
   //          const Scalar mCa,
   //          // const Scalar mNH4,
   //          const Scalar mHCO3,
   //          const Scalar mCO3,
   //          const Scalar mCl,
   //          const Scalar temperature)
   //  {
   //   // Scalar Ksp = Appa_Ksp( mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, temperature);
   //   Scalar Ksp = Appa_Ksp( mNa,  mCa,  mHCO3,  mCO3,  mCl, temperature);
   //   // Scalar Omega_ = mCa * mCO3 / Ksp;
   //   Scalar Omega_ = mCa * mCO3 / 100; //adde by du 
   //   return Omega_;
   //  }
   //  Scalar OmegaApprox(const Scalar mCa,
   //          const Scalar mCO3)
   //  {
   //   Scalar Omega_ = mCa * mCO3 / 100; // pow (10.,-8.48); adde by du // = 3.3e-9= Ksp(Standard) = pow (10.,-8.48);
   //   return Omega_;
   //  }
   //  Scalar rdiss(const Scalar initialPorosity,
   //          const Scalar volFracCalcite,
   //          const Scalar mNa,
   //          const Scalar mCa,
   //          // const Scalar mNH4,
   //          const Scalar mHCO3,
   //          const Scalar mCO3,
   //          const Scalar mCl,
   //          const Scalar temperature,
   //          const Scalar mH)
   //  {
   //      Scalar Asw = Asw0_ * cbrt((1-volFracCalcite/initialPorosity)*(1-volFracCalcite/initialPorosity));   // TODO Asw should be a function of Sw, too!
   //      if (Asw < 1e-8 || isnan(Asw))
   //      {
   //          std::cout<< "Asw = "<<Asw<<std::endl;
   //          Asw = 0;
   //          std::cout<< "Asw, corrected = "<<Asw<<std::endl;
   //      }
   //      Scalar Acw = ac_ * volFracCalcite;
   //      if (ac_ * volFracCalcite > Asw)
   //          Acw = Asw;
   //      Scalar rdiss_ = 0;
   //      // Scalar Omega_ = Omega(mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, temperature);
   //      Scalar Omega_ = Omega(mNa,  mCa,  mHCO3,  mCO3,  mCl, temperature);
   //      if (Omega_ <= 1)
   //      {
   //            rdiss_ = (kdiss1_ * mH + kdiss2_) * Acw * pow((1 - Omega_),ndiss_); //[mol/dm³s]
   //            rdiss_ *= 1000; // rdiss [mol/m³s]
   //      }
   // 
   //      return rdiss_;
   //  }
   //  Scalar rprec(const Scalar initialPorosity,
   //          const Scalar volFracCalcite,
   //          const Scalar mNa,
   //          const Scalar mCa,
   //          // const Scalar mNH4,
   //          const Scalar mHCO3,
   //          const Scalar mCO3,
   //          const Scalar mCl,
   //          const Scalar temperature)
   //  {
   //      Scalar Asw = Asw0_ * cbrt((1-volFracCalcite/initialPorosity)*(1-volFracCalcite/initialPorosity));   // TODO Asw should be a function of Sw, too!
   //           if (Asw < 1e-8 || isnan(Asw))
   //           {
   //               std::cout<< "Asw = "<<Asw<<std::endl;
   //               Asw = 0;
   //               std::cout<< "Asw, corrected = "<<Asw<<std::endl;
   //           }
   //           Scalar rprec_ = 0;
   //           // Scalar Omega_ = Omega(mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, temperature);
   //           Scalar Omega_ = Omega(mNa,  mCa,  mHCO3,  mCO3,  mCl, temperature);
   //           if (Omega_ >= 1)
   //           {
   //               rprec_ = kprec_ * Asw * pow(Omega_ - 1 , nprec_);//[mol/dm³s]
   //               rprec_ *= 1000; // rprec [mol/m³s]
   //           }
   // 
   // 
   //      return rprec_;
   //  }
   //  // Scalar Fe2Omega, consider the activity, added by DU
   //  Scalar Fe2OmegaApprox(const Scalar mFe2,
   //          const Scalar mOH)
   //  {
   //   // Scalar Fe2Omega_ = mFe2 * mOH * mOH / pow (10.,-4.89); // = 3.3e-9= Ksp(Standard) = pow (10.,-8.48);
   //   Scalar Fe2Omega_ = mFe2 * mOH * mOH / 100; // added by du = 3.3e-9= Ksp(Standard) = pow (10.,-8.48);
   //   return Fe2Omega_;
   //  }
   //  Scalar frdiss(const Scalar initialPorosity,
   //          // const Scalar volFracCalcite,
   //          const Scalar volFracFerrohydrite,
   //          const Scalar mNa,
   //          const Scalar mCa,
   //          // const Scalar mNH4,
   //          const Scalar mHCO3,
   //          const Scalar mCO3,
   //          const Scalar mCl,
   //          const Scalar mFe2,
   //          const Scalar mOH,
   //          const Scalar temperature,
   //          const Scalar mH)
   //  {
   //      // volFracFerrohydrite=0;// added by du
	//	Scalar fAsw = fAsw0_ * cbrt((1-volFracFerrohydrite/initialPorosity)*(1-volFracFerrohydrite/initialPorosity));   // TODO Asw should be a function of Sw, too!
   //      if (fAsw < 1e-8 || isnan(fAsw))
   //      {
   //          std::cout<< "fAsw = "<<fAsw<<std::endl;
   //          fAsw = 0;
   //          std::cout<< "fAsw, corrected = "<<fAsw<<std::endl;
   //      }
   //      Scalar fAcw = fac_ * volFracFerrohydrite;
   //      if (fac_ * volFracFerrohydrite > fAsw)
   //          fAcw = fAsw;
   //      Scalar frdiss_ = 0;
   //      Scalar Fe2OmegaApprox_ = Fe2OmegaApprox(mFe2,  mOH);
   //      if (Fe2OmegaApprox_ <= 1)
   //      {
   //            frdiss_ = (fkdiss1_ * mH + fkdiss2_) * fAcw * pow((1 - Fe2OmegaApprox_),fndiss_); //[mol/dm³s]
   //            frdiss_ *= 1000; // rdiss [mol/m³s]
	//		  
   //      }
   //      // frdiss_ = 0.0; // added by Du, delete Fe precipitation
   //      return frdiss_;
   //  }
   //  Scalar frprec(const Scalar initialPorosity,
   //          // const Scalar volFracCalcite,
   //          const Scalar volFracFerrohydrite,
   //          const Scalar mNa,
   //          const Scalar mCa,
   //          // const Scalar mNH4,
   //          const Scalar mHCO3,
   //          const Scalar mCO3,
   //          const Scalar mCl,
   //          const Scalar mFe2,
   //          const Scalar mOH,
   //          const Scalar temperature)
   //  {
   //      // volFracFerrohydrite=0; // added by du
	//	Scalar fAsw = fAsw0_ * cbrt((1-volFracFerrohydrite/initialPorosity)*(1-volFracFerrohydrite/initialPorosity));   // TODO Asw should be a function of Sw, too!
   //           if (fAsw < 1e-8 || isnan(fAsw))
   //           {
   //               std::cout<< "fAsw = "<<fAsw<<std::endl;
   //               fAsw = 0;
   //               std::cout<< "fAsw, corrected = "<<fAsw<<std::endl;
   //           }
   //           Scalar frprec_ = 0;
   //           Scalar Fe2OmegaApprox_ = Fe2OmegaApprox(mFe2,  mOH);
   //           if (Fe2OmegaApprox_ >= 1)
   //           {
   //               frprec_ = fkprec_ * fAsw * pow(Fe2OmegaApprox_ - 1 , fnprec_);//[mol/dm³s]
   //               frprec_ *= 1000; // rprec [mol/m³s]
   //           }
   // 
   //      // frprec_ = 0.0; // added by Du, delete Fe precipitation
   //      return frprec_;
   //  }

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
        // Scalar initialPorosity = 1.0;
        // for (int i=numSolidComponents-numInertComponents; i<numSolidComponents ; ++i)
        // for (int i=0; i<numSolidComponents ; ++i)
        // {
        //     initialPorosity   -= volVars.solidVolumeFraction(i);
		//     
        // }
		// Scalar initialvolFracGlass = volVars.solidVolumeFraction(gPhaseIdx);
        
		const Scalar initialvolFracGlass = 6.1e-1; // volVars.solidVolumeFraction(gPhaseIdx);
		// const Scalar initialvolFracCalcite = 1.0e-20; // volVars.solidVolumeFraction(gPhaseIdx);
        const Scalar initialvolFracFerrohydrite = 1e-30;
		const Scalar initialvolFracProtoImogolite = 1e-30;
		const Scalar initialvolFracBirnessite = 1e-30;
		const Scalar initialvolFracHydroxyapatite = 1e-30;
		const Scalar initialvolFracSepiolite = 1e-30;
		// Scalar tempPorosity = 1.0;  // Assuming an initial value of 1.0
        // for (int i = 0; i < numSolidComponents; ++i)
        // {
        //     tempPorosity -= volVars.solidVolumeFraction(i);
        // }
        const Scalar initialPorosity = 1.0-initialvolFracGlass-initialvolFracFerrohydrite-initialvolFracProtoImogolite- initialvolFracBirnessite-initialvolFracHydroxyapatite-initialvolFracSepiolite; // tempPorosity;
		
        // for (int i=0; i<numSolidComponents-numInertComponents ; ++i)
        // {
        //     initialPorosity   -= volVars.solidVolumeFraction(i);
	    //  
        // }
        // printf("The value of initialPorosity is: %.5e\n", initialPorosity);   
        Scalar Sw   =  volVars.saturation(wPhaseIdx);
        Scalar xWWater = volVars.moleFraction(wPhaseIdx,wCompIdx);
		
        Scalar volFracGlass = volVars.solidVolumeFraction(gPhaseIdx);
        if (volFracGlass < 1e-30) // Added by du
        {// printf("The value of volFracGlass is: %.12e\n", volFracGlass); 
		volFracGlass = 1e-30;}
        
        // Scalar volFracCalcite = volVars.solidVolumeFraction(cPhaseIdx);
        // if (volFracCalcite < 1e-20)
        // volFracCalcite = 1e-20;
		// 
        Scalar volFracFerrohydrite = volVars.solidVolumeFraction(fPhaseIdx);
        if (volFracFerrohydrite < 1e-30) // Added by du
        {
		// printf("The value of volFracFerrohydrite is: %.12e\n", volFracFerrohydrite); 
		volFracFerrohydrite = 1e-30;}
		
        Scalar volFracProtoImogolite = volVars.solidVolumeFraction(pPhaseIdx);
        if (volFracProtoImogolite < 1e-30) // Added by du
        {
		// printf("The value of volFracProtoImogolite is: %.12e\n", volFracProtoImogolite); 
		volFracProtoImogolite = 1e-30;}		
		
        Scalar volFracBirnessite = volVars.solidVolumeFraction(bPhaseIdx);
		if (volFracBirnessite < 1e-30) // Added by du
		{
		// printf("The value of birrdiss is: %.20e\n", birrdiss);
		// printf("The value of birrprec is: %.20e\n", birrprec);
		// printf("The value of volFracBirnessite is: %.20e\n", volFracBirnessite);
        volFracBirnessite = 1e-30;
	    }
		
		Scalar volFracHydroxyapatite = volVars.solidVolumeFraction(hPhaseIdx);
        if (volFracHydroxyapatite < 1e-30) // Added by du
        {// printf("The value of volFracHydroxyapatite is: %.20e\n", volFracHydroxyapatite);
		volFracHydroxyapatite = 1e-30;}
        
		Scalar volFracSepiolite = volVars.solidVolumeFraction(sPhaseIdx);
        if (volFracSepiolite < 1e-30) // Added by du
        {// printf("The value of volFracSepiolite is: %.20e\n", volFracSepiolite);
		volFracSepiolite = 1e-30;}		
        // Scalar massImmUrease = volVars.solidVolumeFraction(uPhaseIdx)*volVars.solidComponentDensity(uPhaseIdx);
        // if (massImmUrease < 0) added by du
        // massImmUrease = 0;

        // Scalar cUrea = volVars.moleFraction(wPhaseIdx, UreaIdx) * rhoMolar * FluidSystem::molarMass(UreaIdx);
        // Scalar cUrea = 0; // added by du
		// Scalar cUrease = volVars.moleFraction(wPhaseIdx, UreaseIdx) * rhoMolar * FluidSystem::molarMass(UreaseIdx);
        // Scalar cUrease = 0; // added by du
        // Scalar xlSalinity = volVars.moleFraction(wPhaseIdx,NaIdx)
                            // + volVars.moleFraction(wPhaseIdx,CaIdx)
                            // + volVars.moleFraction(wPhaseIdx,ClIdx);
		Scalar moleFracWater = volVars.moleFraction(wPhaseIdx,H2OIdx);
        // Scalar mH = moleFracToMolality(volVars.moleFraction(wPhaseIdx,HIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_H/kg_H2O]
        // Scalar mNa = moleFracToMolality(volVars.moleFraction(wPhaseIdx,NaIdx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_sodium/kg_H2O]
        // Scalar mH = moleFracToMolality(volVars.moleFraction(wPhaseIdx,HIdx), moleFracWater);  //[mol_H/kg_H2O]
        // Scalar mNa = moleFracToMolality(volVars.moleFraction(wPhaseIdx,NaIdx), moleFracWater);  //[mol_sodium/kg_H2O]
        Scalar mH = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,HonlyIdx),volVars.molarDensity(wPhaseIdx));  //[mol_H/kg_H2O]
        if (mH < 0 )//|| std::isnan(mH))
		{
			mH = 0;
			// printf("The value of mH is: %.e\n", mH);
		}
        Scalar mNa = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,NaIdx),volVars.molarDensity(wPhaseIdx));  //[mol_sodium/kg_H2O]
        if (mNa < 0)// || std::isnan(mNa))
		{
			mNa = 0;
			// printf("The value of mNa is: %.e\n", mNa);
		}
        // Scalar mCl = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,ClIdx),volVars.molarDensity(wPhaseIdx));  //[mol_chloride/kg_H2O]
        // if (mCl < 0)
        //     mCl = 0;
        // Scalar mNH4 = moleFracToMolality(volVars.moleFraction(wPhaseIdx,NH4Idx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_NH4/kg_H2O]
        // if (mNH4 < 0) added by du
        //     mNH4 = 0;
        Scalar mCa = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,CaIdx),volVars.molarDensity(wPhaseIdx));  //[mol_calcium/kg_H2O]
        if (mCa < 0 )//|| std::isnan(mCa))
		{
			mCa = 0;
			// printf("The value of mCa is: %.e\n", mCa);
		}
        Scalar mCO2aq = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,CO2aqonlyIdx), volVars.molarDensity(wPhaseIdx));  //[mol_CO3/kg_H2O]
        if (mCO2aq < 0)
        {
			mCO2aq = 0;
		}
        Scalar mO2 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,O2Idx), volVars.molarDensity(wPhaseIdx));  //[mol_CO3/kg_H2O]
        if (mO2 < 0)
        {
			mO2 = 0;
		}
        // Scalar mCO3 = moleFracToMolality(volVars.moleFraction(wPhaseIdx,CO3Idx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_CO3/kg_H2O]
        // if (mCO3 < 0)
        //     mCO3 = 0;
        // Scalar mHCO3 = moleFracToMolality(volVars.moleFraction(wPhaseIdx,HCO3Idx), xlSalinity, volVars.moleFraction(wPhaseIdx,nCompIdx));  //[mol_HCO3/kg_H2O]
        // if (mHCO3 < 0)
        //     mHCO3 = 0;
        Scalar mFe2 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,Fe2onlyIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mFe2 < 0 )//|| std::isnan(mFe2)) // added by DU
		{
			mFe2 = 0;
			// printf("The value of mFe2 is: %.e\n", mFe2);
		}
        Scalar mFe3 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,Fe3Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mFe3 < 0 )//|| std::isnan(mFe2)) // added by DU
		{
			mFe3 = 0;
			// printf("The value of mFe2 is: %.e\n", mFe2);
		}
        Scalar mSiO2 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,SiO2Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mSiO2 < 0 )//|| std::isnan(mSiO2)) // added by DU
		{
			mSiO2 = 0;
			// printf("The value of mSiO2 is: %.e\n", mSiO2);
		}
        Scalar mTiOH4 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,TiOH4Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mTiOH4 < 0 )//|| std::isnan(mTiOH4)) // added by DU
		{
			mTiOH4 = 0;
			// printf("The value of mTiOH4 is: %.e\n", mTiOH4);
		}
        Scalar mAl = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,AlIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mAl < 0 )//|| std::isnan(mAl)) // added by DU
		{
			mAl = 0;
			// printf("The value of mAl is: %.e\n", mAl);
		}
        Scalar mMn = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,MnIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mMn < 0 )//|| std::isnan(mMn)) // added by DU
		{
			mMn = 0;
			// printf("The value of mMn is: %.e\n", mMn);
		}
        Scalar mMg = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,MgIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mMg < 0 )// || std::isnan(mMg)) // added by DU
		{
			mMg = 0;
			// printf("The value of mMg is: %.e\n", mMg);
		}
        Scalar mK = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,KIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mK < 0 )//|| std::isnan(mK)) // added by DU
		{
			mK = 0;
			// printf("The value of mK is: %.e\n", mK);
		}
        Scalar mHPO4 = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,HPO4Idx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mHPO4 < 0)// || std::isnan(mHPO4)) // added by DU
		{
			mHPO4 = 0;
			// printf("The value of mHPO4 is: %.e\n", mHPO4); 
		}
        Scalar mOH = moleFracToMolarity(volVars.moleFraction(wPhaseIdx,OHIdx),volVars.molarDensity(wPhaseIdx));  //[mol_HCO3/kg_H2O]
        if (mOH < 0)
		{
			mOH = 0;
		}
        // compute dissolution and precipitation rate of basalt glass
        Scalar glsp = const4(pressure_, temperature_);
        Scalar glOmegaApprox_ = pow(mSiO2,1.8) * pow(mTiOH4,0.07) * pow(mAl,0.62) * pow(mFe2,0.38) * pow(mMn,0.01)* pow(mMg,0.3)
		* pow(mCa,0.44)* pow(mNa,0.26)* pow(mK,0.06)* pow(mHPO4,0.03)/ pow (mH,4.38) /glsp;
       // if (std::isnan(glOmegaApprox_)) // added by DU
		//{
		//	glOmegaApprox_ = 0;
		//}
		// printf("The value of glsp is: %.8e\n", glsp); 
		// printf("The value of glOmegaApprox1_ is: %.8e\n", glOmegaApprox_);    
		// printf("The value of mHPO4 is: %.7e\n", mHPO4); 
		// printf("The value of mAl is: %.7e\n", mAl); 
		// printf("The value of mH is: %.7e\n", mH); //0.0
		// printf("The value of mSiO2 is: %.7e\n", mSiO2); 
		// printf("The value of mFe is: %.7e\n", mFe2); 
		// printf("The value of mNa is: %.7e\n", mNa); 
		// printf("The value of mK is: %.7e\n", mK); 
		// printf("The value of mMg is: %.7e\n", mMg); 
		// printf("The value of mCa is: %.7e\n", mCa);
		// printf("The value of mMn is: %.7e\n", mMn); 			
		// printf("The value of mTiOH4 is: %.7e\n", mTiOH4); // -Nan		
		// Scalar glAw = glAsw_ * initialvolFracGlass * volVars.solidComponentDensity(gPhaseIdx);
		Scalar glAw0 = glAsw_ * initialvolFracGlass * volVars.solidComponentDensity(gPhaseIdx);
		// Scalar glAsw = glAsw0_ * cbrt((1-volFracGlass/initialPorosity)*(1-volFracGlass/initialPorosity));   // TODO Asw should be a function of Sw, too!
 		Scalar glAwcd = glAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracGlass/initialvolFracGlass)*(volFracGlass/initialvolFracGlass));   // TODO Asw should be a function of Sw, too!
 		Scalar glAwcp = glAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!

		// Scalar glAwf = glAsw_ * volFracGlass * volVars.solidComponentDensity(gPhaseIdx);
		// Scalar glAw = glAwc;
		// if (glAw > glAwf)
		// {glAw = glAwf;}
		// Scalar glAsw = glAsw0_ * cbrt((1-volFracGlass/initialPorosity)*(1-volFracGlass/initialPorosity));   // TODO Asw should be a function of Sw, too!
 		// Scalar glAw = glAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracGlass/initialvolFracGlass)*(volFracGlass/initialvolFracGlass));   // TODO Asw should be a function of Sw, too!
		// printf("The value of volFracGlass is: %.e\n", volFracGlass); 
		// printf("The value of initialPorosity is: %.e\n", initialPorosity); 
		// printf("The value of Porosityratio is: %.e\n", 1-volFracGlass/initialPorosity); 
		// printf("The value of volFracGlassratio is: %.e\n", volFracGlass/initialvolFracGlass); 
		// printf("The value of glAw0_ is: %.10e\n", glAw0); 
		// printf("The value of glAw is: %.10e\n", glAw); 		
        // if (glAsw < 1e-8 || std::isnan(glAsw)) // changed by du
        if (glAwcd < 0 || std::isnan(glAwcd))
        {
		// printf("The value of volFracGlass is: %.e\n", volFracGlass); 
		// printf("The value of initialPorosity is: %.e\n", initialPorosity); 
		// printf("The value of initialvolFracGlass is: %.e\n", initialvolFracGlass); 		
		// printf("The value of Porosityratio is: %.e\n", 1-volFracGlass/initialPorosity); 
		// printf("The value of volFracGlassratio is: %.e\n", volFracGlass/initialvolFracGlass); 
		// printf("The value of glAsw0_ is: %.e\n", glAsw0_); 
		// printf("The value of glAsw is: %.e\n", glAsw); 	
        std::cout<< "glAwcd = "<<glAwcd<<std::endl;
        glAwcd = 0;
        std::cout<< "glAwcd, corrected = "<<glAwcd<<std::endl;
        }

        if (glAwcp < 0 || std::isnan(glAwcp))
        {
		// printf("The value of volFracGlass is: %.e\n", volFracGlass); 
		// printf("The value of initialPorosity is: %.e\n", initialPorosity); 
		// printf("The value of initialvolFracGlass is: %.e\n", initialvolFracGlass); 		
		// printf("The value of Porosityratio is: %.e\n", 1-volFracGlass/initialPorosity); 
		// printf("The value of volFracGlassratio is: %.e\n", volFracGlass/initialvolFracGlass); 
		// printf("The value of glAsw0_ is: %.e\n", glAsw0_); 
		// printf("The value of glAsw is: %.e\n", glAsw); 	
        std::cout<< "glAwcp = "<<glAwcp<<std::endl;
        glAwcp = 0;
        std::cout<< "glAwcp, corrected = "<<glAwcp<<std::endl;
        }
		
		Scalar glr = 0;
		// glr = - glAsw * pow(10,glrc_) * glp_ * pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_);
			if ( glOmegaApprox_ < 1 && glr > 1e-6)
			{
			// printf("The value of glr is: %.e\n", glr); 
			// printf("The value of glAsw is: %.e\n", glAsw); 
		    // printf("The value of Porosityratio is: %.e\n", (1-volFracGlass)/initialPorosity); 
		    // printf("The value of volFracGlassratio is: %.e\n", volFracGlass/initialvolFracGlass);
            // printf("The value of 1-glOmegaApprox_ is: %.e\n",pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_));	// Nan	
		    // printf("The value of glAsw0_ is: %.e\n", glAsw0_);
			// printf("The value of glOmegaApprox2_ is: %.7e\n", glOmegaApprox_); //-nan
			// printf("The value of mHPO4 is: %.7e\n", mHPO4); 
			// printf("The value of mAl is: %.7e\n", mAl); 
			// printf("The value of mH is: %.7e\n", mH); //0.0
			// printf("The value of mSiO2 is: %.7e\n", mSiO2); 
			// printf("The value of mFe is: %.7e\n", mFe2); 
			// printf("The value of mNa is: %.7e\n", mNa); 
			// printf("The value of mK is: %.7e\n", mK); 
			// printf("The value of mMg is: %.7e\n", mMg); 
			// printf("The value of mCa is: %.7e\n", mCa);
			// printf("The value of mMn is: %.7e\n", mMn); 			
			// printf("The value of mTiOH4 is: %.7e\n", mTiOH4); // -Nan
			}		
        // Scalar glAcw = glac_ * volFracGlass;
        // if (glac_ * volFracGlass > glAsw)
        //     glAcw = glAsw;
		// 
        Scalar glrdiss = 0;
        Scalar glrprec = 0;
		
        if (glOmegaApprox_ > 1)
        {
		    glr = glAwcp * pow(10,glrc_) * glp_ * pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_);
            // printf("Precipitation");   
            glrdiss = 0;
            glrprec = glr;//[mol/m³s] // [mol/kg s]
			// glrprec *= 0.00001;
            // glrprec *= 1000; // rprec [mol/m³s]
		    // grprec = 0; // added by du
			// printf("The value of glOmegaApprox2_ is: %.7e\n", glOmegaApprox_);
			// printf("The value of mHPO4 is: %.7e\n", mHPO4); 
			// printf("The value of mAl is: %.7e\n", mAl); 
			// printf("The value of mH is: %.7e\n", mH); 
			// printf("The value of mSiO2 is: %.7e\n", mSiO2); 
			// printf("The value of mFe is: %.7e\n", mFe2); 
			// printf("The value of mNa is: %.7e\n", mNa); 
			// printf("The value of mK is: %.7e\n", mK); 
			// printf("The value of mMg is: %.7e\n", mMg); 
			// printf("The value of mCa is: %.7e\n", mCa);
			// printf("The value of mMn is: %.7e\n", mMn); 			
			// printf("The value of mTiOH4 is: %.7e\n", mTiOH4); 
		    // printf("The value of Porosityratio is: %.e\n", (1-volFracGlass)/initialPorosity); 
		    // printf("The value of volFracGlassratio is: %.e\n", volFracGlass/initialvolFracGlass);
            // printf("The value of 1-glOmegaApprox_ is: %.e\n",pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_));		
		    // printf("The value of glAsw0_ is: %.e\n", glAsw0_); 
		    // printf("The value of glAsw is: %.e\n", glAsw); 	
	
        }
        else
        {
		    glr = glAwcd * pow(10,glrc_) * glp_ * pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_);
		    // printf("Dissolution");   
			glrdiss = glr; //[mol/kg s]
			// glrdiss *= 0.00001;
			// glrdiss *= 1000; // rdiss [mol/m³s]
            glrprec = 0;
			// printf("The value of glAwcd is: %.7e\n", glAwcd);
		    // printf("The value of glr is: %.7e\n", glr); 
			if (std::isnan(glr))
			{
			// printf("The value of glr is: %.e\n", glr); 
			// printf("The value of glAsw is: %.e\n", glAsw); 
		    // printf("The value of Porosityratio is: %.e\n", (1-volFracGlass)/initialPorosity); 
		    // printf("The value of volFracGlassratio is: %.e\n", volFracGlass/initialvolFracGlass);
            // printf("The value of 1-glOmegaApprox_ is: %.e\n",pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_));	// Nan	
		    // printf("The value of glAsw0_ is: %.e\n", glAsw0_);
			// printf("The value of glOmegaApprox2_ is: %.7e\n", glOmegaApprox_); //-nan
			// printf("The value of mHPO4 is: %.7e\n", mHPO4); 
			// printf("The value of mAl is: %.7e\n", mAl); 
			// printf("The value of mH is: %.7e\n", mH); //0.0
			// printf("The value of mSiO2 is: %.7e\n", mSiO2); 
			// printf("The value of mFe is: %.7e\n", mFe2); 
			// printf("The value of mNa is: %.7e\n", mNa); 
			// printf("The value of mK is: %.7e\n", mK); 
			// printf("The value of mMg is: %.7e\n", mMg); 
			// printf("The value of mCa is: %.7e\n", mCa);
			// printf("The value of mMn is: %.7e\n", mMn); 			
			// printf("The value of mTiOH4 is: %.7e\n", mTiOH4); // -Nan

			}
			// printf("The value of mNa is: %.7e\n", mNa);
        }

	    // printf("The value of glr is: %.e\n", glr); 
	    // printf("The value of glAsw is: %.e\n", glAsw); 
	    // printf("The value of Porosityratio is: %.e\n", (1-volFracGlass)/initialPorosity); 
	    // printf("The value of volFracGlassratio is: %.e\n", volFracGlass/initialvolFracGlass);
        // printf("The value of 1-glOmegaApprox_ is: %.e\n",pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_));	// Nan	
	    // printf("The value of glAsw0_ is: %.e\n", glAsw0_);
	    // printf("The value of glOmegaApprox2_ is: %.7e\n", glOmegaApprox_); //-nan
		
        // // calcite
        // Scalar calsp = const5(pressure_, temperature_);
        // Scalar calOmegaApprox_ = pow(mCO2aq,1) * pow(mCa,1) / pow (mH,2) /calsp;
		// Scalar calAw0 = calAsw_ * volFracCalcite * volVars.solidComponentDensity(cPhaseIdx);
 		// Scalar calAw = calAw0 * cbrt(((1-volFracCalcite-volFracGlass)/initialPorosity)*((1-volFracCalcite-volFracGlass)/initialPorosity))* cbrt((volFracCalcite/initialvolFracCalcite)*(volFracCalcite/initialvolFracCalcite));   // TODO Asw should be a function of Sw, too!
		// 
        // if (calAw < 0 || std::isnan(calAw))
        // {
        // std::cout<< "calAw = "<<calAw<<std::endl;
        // calAw = 0;
        // std::cout<< "calAw, corrected = "<<calAw<<std::endl;
        // }
		// 
		// Scalar calr = 0;
		// calr = calAw * pow(10,calrc_) * calp_ * pow(abs(1-pow(calOmegaApprox_,1/calsigma_)),calbeta_);
		// 	if ( calOmegaApprox_ < 1 && calr > 1e-6)
		// 	{
		// 	}		
		// 
        // Scalar calrdiss = 0;
        // Scalar calrprec = 0;
		// 
        // if (calOmegaApprox_ > 1)
        // {
		//     // printf("Precipitation");   
        //     calrdiss = 0;
        //     calrprec = calr;//[mol/kg s]
		// 	// calrprec *= 0.00001;
        //     // calrprec *= 1000; // rprec [mol/m³s]
        // }
        // else
        // {
		//     // printf("Dissolution");   
		// 	calrdiss = calr; //[mol/kg s]
		// 	// calrdiss *= 0.00001;
		// 	// calrdiss *= 1000; // rdiss [mol/m³s]
        //     calrprec = 0;
		// 	if (std::isnan(calr))
		// 	{
		// 	}
        // }
		
        // ferrohydrite
        Scalar fersp = const6(pressure_, temperature_);
        Scalar ferOmegaApprox_ = pow(mFe3,1) / pow (mH,3) /fersp;
        Scalar ferAw0 = ferAsw_ * initialvolFracFerrohydrite * volVars.solidComponentDensity(fPhaseIdx);
 		Scalar ferAwcp = ferAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!
 		Scalar ferAwcd = ferAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracFerrohydrite*volFracFerrohydrite));   // TODO Asw should be a function of Sw, too!
		Scalar ferAwf = ferAsw_ * volFracFerrohydrite * volVars.solidComponentDensity(fPhaseIdx);
		
		Scalar ferAwp = ferAwcp;
		Scalar ferAwd = ferAwcd;
		
		// if (ferAw > ferAwf)
		// {ferAw = ferAwf;}		

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
		// if ( ferOmegaApprox_ < 1 && ferr > 1e-6)
		// {
		// }		

        Scalar ferrdiss = 0;
        Scalar ferrprec = 0;
		
        if (ferOmegaApprox_ > 1)
        {
		    Scalar ferr = ferAwp * pow(10,ferrc_) * ferp_ * pow(abs(1-pow(ferOmegaApprox_,1/fersigma_)),ferbeta_);
			// printf("Precipitation");   
            ferrdiss = 0;
            ferrprec = ferr;//[mol/kg s]
        }
        else
        {
		    Scalar ferr = ferAwd * pow(10,ferrc_) * ferp_ * pow(abs(1-pow(ferOmegaApprox_,1/fersigma_)),ferbeta_);
			// printf("Dissolution");   
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
		// printf("The value of ferOmegaApprox2_ is: %.7e\n", ferOmegaApprox_); //-nan
		
	    // printf("The value of mFe is: %.7e\n", mFe2); 
	    // printf("The value of mH is: %.7e\n", mH);
		// printf("The value of ferAw is: %.7e\n", ferAw);
		
		// if (0.02 < dt && dt < 0.03)
		// {
		// printf("The value of ferOmegaApprox_ is: %.12e\n", ferOmegaApprox_);
		// printf("The value of ferrdt is: %.12e\n", ferrdt);
		// printf("The value of ferrprec is: %.12e\n", ferrprec);
		// printf("The value of conc is: %.12e\n",volVars.moleFraction(wPhaseIdx,Fe3Idx) * Sw * porosity * volVars.molarDensity(wPhaseIdx));
        // }
		
        // ProtoImogolite
        Scalar prosp = const7(pressure_, temperature_);
        Scalar proOmegaApprox_ = pow(mAl,2) * pow(mSiO2,1) / pow (mH,6) /prosp;

        Scalar proAw0 = proAsw_ * initialvolFracProtoImogolite * volVars.solidComponentDensity(pPhaseIdx);
		Scalar proAwf = proAsw_ * volFracProtoImogolite * volVars.solidComponentDensity(pPhaseIdx);
 		Scalar proAwcd = proAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracProtoImogolite*volFracProtoImogolite));   // TODO Asw should be a function of Sw, too!
 		Scalar proAwcp = proAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!

		Scalar proAwd = proAwcd;
		Scalar proAwp = proAwcp;		
		// if (proAw > proAwf)
		// {proAw = proAwf;}
		// printf("The value of proAw1 is: %.10e\n", proAw);	
	    // proAw = 6 * cbrt(volFracProtoImogolite*volFracProtoImogolite) * cbrt(volVars.solidComponentDensity(pPhaseIdx));
		// printf("The value of proAw0_ is: %.10e\n", proAw0); 
		// printf("The value of proAwc_ is: %.10e\n", proAwc); 
		// printf("The value of proAwd_ is: %.10e\n", proAwd); 
		// printf("The value of proAwp is: %.10e\n", proAwp);
	
	
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
		    // printf("Precipitation");   
            prordiss = 0;
            prorprec = pror;//[mol/kg s]
			// printf("The value of proAwp is: %.7e\n", proAwp);
		    // printf("The value of pror is: %.7e\n",  pror); 
        }
        else
        {
			pror = proAwd * pow(10,prorc_) * prop_ * pow(abs(1-pow(proOmegaApprox_,1/prosigma_)),probeta_);

		    // printf("Dissolution");   
			prordiss = pror; //[mol/kg s]
            prorprec = 0;
			if (std::isnan(pror))
			{
			}
        }
		
		// if (dt == 0.01)
		// {
		// printf("The value of prorprec1 is: %.12e\n", prorprec); 
		// }
		// if (dt == 0.015)
		// {
		// printf("The value of prorprec2 is: %.12e\n", prorprec); 
		// }
		
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

		
		// if (dt == 0.015)
		// {
		// printf("The value of prorprec2 is: %.12e\n", prorprec); 
		// }

		// if (dt > 0.0168)
		// {
		// printf("The value of prorprecdt is: %.12e\n", prorprecdt); 
		// printf("The value of prorprec is: %.12e\n", prorprec);
		// printf("The value of conc is: %.12e\n",volVars.moleFraction(wPhaseIdx,AlIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx));
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
		
		// if (0.02 < dt && dt < 0.03)
		// {
		// printf("The value of prorprecdt is: %.12e\n", prorprecdt); 
		// printf("The value of prordissdt is: %.12e\n", prordissdt); 
		// printf("The value of proOmegaApprox_ is: %.12e\n", proOmegaApprox_);
		// printf("The value of prordt is: %.12e\n", prordt);
		// printf("The value of prorprec is: %.12e\n", prorprec);
		// printf("The value of conc is: %.12e\n",volVars.moleFraction(wPhaseIdx,AlIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx));
        // }
		// if (dt < 0.000038)
		// {
		// printf("The value of prordt is: %.12e\n", prordt); 
		// }		
		
		// FILE *file; // Declare a file pointer
		// 
        // // Open the file in write mode ("w")
        // file = fopen("output.txt", "a");
	    // fprintf(file, "The value of prolr is: %.e\n", pror); 
		// fprintf(file, "The value of proAw is: %.e\n", proAw); 
        // fprintf(file, "The value of 1-proOmegaApprox_ is: %.e\n",pow(abs(1-pow(proOmegaApprox_,1/prosigma_)),probeta_));	// Nan	
		// fprintf(file, "The value of proAw0_ is: %.e\n", proAw0);
		// fprintf(file, "The value of proOmegaApprox2_ is: %.7e\n", proOmegaApprox_); //-nan
		// fprintf(file, "The value of mAl is: %.7e\n", mAl); 
		// fprintf(file, "The value of mSiO2 is: %.7e\n", mSiO2);
		// fprintf(file, "The value of mH is: %.7e\n", mH);
		// fprintf(file, "The value of q[phiProtoImogoliteIdx] is: %.7e\n", prorprec - prordiss);
		// // Write unformatted text (string) to the file
        // fputs("This is another line written using fputs.\n", file);
		// 
        // // Close the file
        // fclose(file);
		
        // Birnessite
        Scalar birsp = const8(pressure_, temperature_);
        Scalar birOmegaApprox_ = pow(mMn,8) * pow(mO2,3) / pow (mH,16) /birsp;
		Scalar birAw0 = birAsw_ * initialvolFracBirnessite * volVars.solidComponentDensity(bPhaseIdx);
 		Scalar birAwf = birAsw_ * volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx);
 		Scalar birAwcp = birAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   // TODO Asw should be a function of Sw, too!
 		Scalar birAwcd = birAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracBirnessite)*(volFracBirnessite));   // TODO Asw should be a function of Sw, too!

		Scalar birAwd = birAwcd;
		Scalar birAwp = birAwcp;		
		// if (birAw > birAwf)
		// {birAw = birAwf;}

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
		// if ( birOmegaApprox_ < 1 && birr > 1e-6)
		// {
		// }		

        Scalar birrdiss = 0;
        Scalar birrprec = 0;
		
        if (birOmegaApprox_ > 1)
        {
			Scalar birr = birAwp * pow(10,birrc_) * birp_ * pow(abs(1-pow(birOmegaApprox_,1/birsigma_)),birbeta_);
		    // printf("Precipitation");   
            birrdiss = 0;
            birrprec = birr;//[mol/kg s]
        }
        else
        {
		    Scalar birr = birAwd * pow(10,birrc_) * birp_ * pow(abs(1-pow(birOmegaApprox_,1/birsigma_)),birbeta_);
			// printf("Dissolution");   
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

		// if (0.02 < dt && dt < 0.03)
		// {
		// printf("The value of birOmegaApprox_ is: %.12e\n", birOmegaApprox_);
		// printf("The value of birrdt is: %.12e\n", birrdt);
		// printf("The value of birrprec is: %.12e\n", birrprec);
		// printf("The value of conc is: %.12e\n",volVars.moleFraction(wPhaseIdx,MnIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx));
        // }
		// if(birrdiss > volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / dt)
        // {
        //     printf("The value of birrdiss is: %.12e\n", birrdiss);
		// 	printf("The value of volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / dt is: %.12e\n", volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / dt);
		// 	printf("The value of volFracBirnessite is: %.7e\n", volFracBirnessite);
		// 	birrdiss =  volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / dt;
		// 	 
        // }	
		// if (q[phiBirnessiteIdx] < 0)
        // {
		// printf("The value of birrdiss is: %.8e\n", birrdiss); 
		// printf("The value of birrdiss2 is: %.8e\n", volFracBirnessite * volVars.solidComponentDensity(bPhaseIdx) / dt); 
		// printf("The value of birrprec is: %.8e\n", birrprec); 
		// printf("The value of volFracBirnessite is: %.7e\n", volFracBirnessite);
		// // printf("The value of birlr is: %.e\n", birr); 
		// //printf("The value of birAw is: %.e\n", birAw); 
        // //printf("The value of 1-birOmegaApprox_ is: %.e\n",pow(abs(1-pow(birOmegaApprox_,1/birsigma_)),birbeta_));	// Nan	
		// //printf("The value of birAw0_ is: %.e\n", birAw0);
		// printf("The value of birOmegaApprox2_ is: %.7e\n", birOmegaApprox_); //-nan
		// }
		// printf("The value of volFracBirnessite is: %.7e\n", volFracBirnessite);
		// printf("The value of q[phiBirnessiteIdx] is: %.7e\n", q[phiBirnessiteIdx]);
		// printf("The value of birrdiss is: %.8e\n", birrdiss); 
		// printf("The value of birrprec is: %.8e\n", birrprec); 
		// // 
		// printf("The value of mO2 is: %.7e\n", mO2); 
		// printf("The value of mH is: %.7e\n", mH);
		// printf("The value of mMn is: %.7e\n", mMn);
		// printf("The value of birAw is: %.e\n", birAw);
		
        // Hydroxyapatite
        Scalar hydsp = const9(pressure_, temperature_);
        Scalar hydOmegaApprox_ = pow(mCa,5) * pow(mHPO4,3) / pow (mH,4) /hydsp;
		Scalar hydAw0 = hydAsw_ * initialvolFracHydroxyapatite * volVars.solidComponentDensity(hPhaseIdx);
 		Scalar hydAwf = hydAsw_ * volFracHydroxyapatite * volVars.solidComponentDensity(hPhaseIdx);
 		Scalar hydAwcd = hydAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracHydroxyapatite)*(volFracHydroxyapatite));   
 		Scalar hydAwcp = hydAw0 * cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity));   

		Scalar hydAwd = hydAwcd;
		Scalar hydAwp = hydAwcp;
		
		// if (hydAw > hydAwf)
		// {hydAw = hydAwf;}
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
		
		
		// if ( hydOmegaApprox_ < 1 && hydr > 1e-6)
		// {
		// }		

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
		
		// if (0.02 < dt && dt < 0.03)
		// {
		// printf("The value of hydOmegaApprox_ is: %.12e\n", hydOmegaApprox_);
		// printf("The value of hydrdt is: %.12e\n", hydrdt);
		// printf("The value of hydrprec is: %.12e\n", hydrprec);
		// printf("The value of conc is: %.12e\n",volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx));
        // }
		// printf("The value of hydlr is: %.e\n", hydr); 
		// printf("The value of hydAw is: %.e\n", hydAw); 
        // printf("The value of 1-hydOmegaApprox_ is: %.e\n",pow(abs(1-pow(hydOmegaApprox_,1/hydsigma_)),hydbeta_));	// Nan	
		// printf("The value of hydAw0_ is: %.e\n", hydAw0);
		// printf("The value of hydOmegaApprox2_ is: %.7e\n", hydOmegaApprox_); //-nan
		// printf("The value of mCa is: %.7e\n", mCa); 
		// printf("The value of mH is: %.7e\n", mH);
		// printf("The value of mHPO4 is: %.7e\n", mHPO4);
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
		
		// if ( sepOmegaApprox_ < 1 && sepr > 1e-6)
		// {
		// }		

        Scalar seprdiss = 0;
        Scalar seprprec = 0;
		
        if (sepOmegaApprox_ > 1)
        {
		    Scalar sepr = sepAwp * pow(10,seprc_) * sepp_ * pow(abs(1-pow(sepOmegaApprox_,1/sepsigma_)),sepbeta_);
			// printf("Precipitation");   
            seprdiss = 0;
            seprprec = sepr;//[mol/kg s]
        }
        else
        {
		    Scalar sepr = sepAwd * pow(10,seprc_) * sepp_ * pow(abs(1-pow(sepOmegaApprox_,1/sepsigma_)),sepbeta_);
			// printf("Dissolution");   
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
	
		// if (0.02 < dt && dt < 0.03)
		// {
		// printf("The value of sepOmegaApprox_ is: %.12e\n", sepOmegaApprox_);
		// printf("The value of seprdt is: %.12e\n", seprdt);
		// printf("The value of seprprec is: %.12e\n", seprprec);
		// printf("The value of conc is: %.12e\n",volVars.moleFraction(wPhaseIdx,MgIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx));
        // }
	    // printf("The value of sepOmegaApprox2_ is: %.7e\n", sepOmegaApprox_); //-nan
	    // printf("The value of mMg is: %.7e\n", mMg); 
	    // printf("The value of mH is: %.7e\n", mH);
	    // printf("The value of mSiO2 is: %.7e\n", mSiO2);		
		// 
		// // need to figure out this part
        // if (glOmegaApprox_ >= 1)
        // {
		//     // printf("Precipitation");   
        //     glrdiss = 0;
        //     glrprec = glkprec_ * glAsw * pow(glOmegaApprox_ - 1 , glnprec_);//[mol/dm³s]
        //     glrprec *= 1000; // rprec [mol/m³s]
		//     // grprec = 0; // added by du
        // }
        // else
        // {
		//     // printf("Dissolution");   
		// 	glrdiss = (glkdiss1_ * mH + glkdiss2_) * glAcw * pow((1 - glOmegaApprox_),glndiss_); //[mol/dm³s]
		// 	glrdiss *= 1000; // rdiss [mol/m³s]
        //     glrprec = 0;
        // }
		
		// printf("The value of glrprec is: %.e\n", glrprec); 
		// printf("The value of glrdiss is: %.e\n", glrdiss); 
		// question?
        // if(glrprec >
        //     volVars.moleFraction(wPhaseIdx,Fe2Idx) * Sw * porosity * rhoMolar / dt)
        // {
        //     glrprec =  volVars.moleFraction(wPhaseIdx,Fe2Idx) * Sw * porosity * rhoMolar / dt;
        // }
		
		
        // // compute dissolution and precipitation rate of ferrohydrite
        // Scalar Fe2OmegaApprox_ = mFe2 * mOH * mOH / pow (10.,-4.89);
        // Scalar fAsw = fAsw0_ * cbrt((1-volFracFerrohydrite/initialPorosity)*(1-volFracFerrohydrite/initialPorosity));   // TODO Asw should be a function of Sw, too!
        // if (fAsw < 1e-8 || std::isnan(fAsw))
        // {
        // std::cout<< "fAsw = "<<fAsw<<std::endl;
        // fAsw = 0;
        // std::cout<< "fAsw, corrected = "<<fAsw<<std::endl;
        // }
        // Scalar fAcw = fac_ * volFracFerrohydrite;
        // if (fac_ * volFracFerrohydrite > fAsw)
        //     fAcw = fAsw;
		// 
        // Scalar frdiss = 0;
        // Scalar frprec = 0;
        // if (Fe2OmegaApprox_ >= 1)
        // {
        // frdiss = 0;
        // frprec = fkprec_ * fAsw * pow(Fe2OmegaApprox_ - 1 , fnprec_);//[mol/dm³s]
        // frprec *= 1000; // rprec [mol/m³s]
		// frprec = 0; // added by du
        // }
        // else
        // {
        //     frprec = 0;
        // }
        // if(frprec >
        //     volVars.moleFraction(wPhaseIdx,Fe2Idx) * Sw * porosity * rhoMolar / dt)
        // {
        //     frprec =  volVars.moleFraction(wPhaseIdx,Fe2Idx) * Sw * porosity * rhoMolar / dt;
        // }
		
		// 
        // // compute dissolution and precipitation rate of calcite
        // // Scalar Ksp = Appa_Ksp( mNa,  mCa,  mNH4,  mHCO3,  mCO3,  mCl, volVars.temperature());
        // Scalar Ksp = Appa_Ksp( mNa,  mCa,  mHCO3,  mCO3,  mCl, volVars.temperature());
        // Scalar Omega = mCa * mCO3 / Ksp;
        // Scalar Asw = Asw0_ * cbrt((1-volFracCalcite/initialPorosity)*(1-volFracCalcite/initialPorosity));   // TODO Asw should be a function of Sw, too!
        // if (Asw < 1e-8 || std::isnan(Asw))
        // {
        // std::cout<< "Asw = "<<Asw<<std::endl;
        // Asw = 0;
        // std::cout<< "Asw, corrected = "<<Asw<<std::endl;
        // }
        // Scalar Acw = ac_ * volFracCalcite;
        // if (ac_ * volFracCalcite > Asw)
        //     Acw = Asw;
		// 
        // Scalar rdiss = 0;
        // Scalar rprec = 0;
        // if (Omega >= 1)
        // {
        // rdiss = 0;
        // rprec = kprec_ * Asw * pow(Omega - 1 , fnprec_);//[mol/dm³s]
        // rprec *= 1000; // rprec [mol/m³s]
        // }
        // else
        // {
//      //        rdiss = (kdiss1_ * mH + kdiss2_) * Acw * pow((1 - Omega),ndiss_); //[mol/dm³s]
//      //        rdiss *= 1000; // rdiss [mol/m³s]
        //     rprec = 0;
        // }
        // if(rprec >
        //     volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * rhoMolar / dt)
        // {
        //     rprec =  volVars.moleFraction(wPhaseIdx,CaIdx) * Sw * porosity * rhoMolar / dt;
        // }

        // // compute first-order rate of ureolysis:
        // // calculate the temperature-dependent rate coefficient
        // Scalar kureaseT = kurease_ * exp(cureaseT_ / temperature);
		// 
        // Scalar rurea_urease = kureaseT * cUrease * porosity * Sw * cUrea;
//      //    Scalar rurea_immUrease = kureaseT * massImmUrease * ureaseInEnzymeSource_ * cUrea;
        // Scalar rurea_immUrease = kureaseT * massImmUrease * cUrea;
		// 
        // //[mol_urea/m³s]
        // rurea_urease /=FluidSystem::molarMass(UreaIdx);
        // rurea_immUrease /=FluidSystem::molarMass(UreaIdx);
		// 
        // //compute rate of enzyme inactivation
        // Scalar kia = cia_ * exp(ciaT_ / temperature); //[1/s]
        // //kg_urease/m³s]
        // Scalar ria_urease = kia * cUrease * porosity * Sw;
        // Scalar ria_immUrease = (kia  +
        //     pow(rprec * SolidSystem::molarMass(cPhaseIdx) /
        //     (volVars.solidComponentDensity(cPhaseIdx) * (initialPorosity - volFracCalcite)), ciaPrec_)
        //                        )
        //     * massImmUrease;
		// 
		// 
        // // compute rate of ureolysis due to temperature:
        // //using actual parameters fitted to experimental data:
        // // calculate the temperature-dependent rate coefficient October 2019 data fitted
        // Scalar ku = cu_ * exp(cuT_ / temperature - 0.5*mCa);
        // Scalar rurea_temp = ku * cUrea/FluidSystem::molarMass(UreaIdx);
		// 
        // //compute the combined rate
        // Scalar rurea = rurea_urease + rurea_immUrease + rurea_temp;
		// 
        // // compute attachment rates:
        // Scalar ra_urease = ka_urease_ * cUrease * porosity * Sw;          //[kg/m³s]
		// 
        // // compute detachment rates:
        // Scalar rd_urease = kd_urease_ * massImmUrease;                      //[kg/m³s]

        // rdiss+rprec[mol/m³s]
        // rurea[mol/m³s]
        // q[kg/m³s]
        // newdt = prordt; // Assume seprdt is the smallest initially
	    
        // if (hydrdt < newdt) {
        //     newdt = hydrdt;
        // }
        // if (birrdt < newdt) {
        //     newdt = birrdt;
        // }
        // if (ferrdt < newdt) {
        //     newdt = ferrdt;
        // }
        // if (prordt < newdt) {
        //     newdt = prordt;
        // }
		
		if (1e-3 < newdt && newdt < 2e-3)
		{
		// printf("The value of newdt is: %.12e\n", newdt); 
		
        // printf("The value of seprdt is: %.12e\n", seprdt); 
		// printf("The value of hydrdt is: %.12e\n", hydrdt); 
		// printf("The value of ferrdt is: %.12e\n", ferrdt); 
		// printf("The value of prordt is: %.12e\n", prordt); 
		
		// printf("The value of prorprecdt is: %.12e\n", prorprecdt); 
		// printf("The value of prordissdt is: %.12e\n", prordissdt); 
		// printf("The value of proOmegaApprox_ is: %.12e\n", proOmegaApprox_);
		// printf("The value of prordt is: %.12e\n", prordt);
		// printf("The value of prorprec is: %.12e\n", prorprec);
		// printf("The value of conc is: %.12e\n",volVars.moleFraction(wPhaseIdx,AlIdx) * Sw * porosity * volVars.molarDensity(wPhaseIdx));
        // }
	    }
		
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
        // q[nCompIdx] += rurea - rprec + rdiss;
        // q[nCompIdx] += - rprec + rdiss; // added by du
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
        // q[HIdx]+= 0;
		q[HIdx]+= 4.38*(+ glrprec - glrdiss)+16*(birrprec - birrdiss)+4*(hydrprec - hydrdiss)+8*(seprprec - seprdiss)+ 2* (ferrprec - ferrdiss) +6*(prorprec - prordiss) ; // added by du +2*(calrprec - calrdiss) 
        q[phiGlassIdx] += (glrprec - glrdiss);
		q[phiFerrohydriteIdx] += (ferrprec - ferrdiss) ;
		q[phiProtoImogoliteIdx] += (prorprec - prordiss) ;
		q[phiBirnessiteIdx] += (birrprec - birrdiss);
		q[phiHydroxyapatiteIdx] += (hydrprec - hydrdiss);
		q[phiSepioliteIdx] += (seprprec - seprdiss);
        // q[phiCalciteIdx] += (calrprec - calrdiss);
		// if (q[phiGlassIdx] = -0.0)
		// {
        // FILE *file = fopen("output.txt", "w");

        // Write to the file
        // printf("The value of volume fraction of Glass is: %.7e\n", volFracGlass);
        // printf("The value of volume fraction of Ferrohydrite is: %.7e\n", volFracFerrohydrite);
        // printf("The value of volume fraction of ProtoImogolite is: %.7e\n", volFracProtoImogolite);
        // printf("The value of q[phiGlassIdx] is: %.7e\n", q[phiGlassIdx]);
		// printf("The value of q[phiFerrohydriteIdx] is: %.7e\n", q[phiFerrohydriteIdx]);
		// if (q[phiProtoImogoliteIdx] < 0)
		// {
        // printf("The value of q[phiProtoImogoliteIdx] is: %.7e\n", q[phiProtoImogoliteIdx]);
		// }
		
		Scalar Hdt = 0;
		// if (0 > q[AlIdx] && mAl < abs(q[AlIdx]))
		
	    if (0 > q[AlIdx])
		{
		prorprecdt = mAl / abs(q[AlIdx]);
		newdt = prorprecdt;
		

		// if (0 > q[HIdx])
		// {
		// Hdt = mH / abs(q[HIdx]);
		// if (newdt > Hdt)
		// {newdt = Hdt;}
		// }
		// printf("The value of q[AlIdx] is: %.7e\n", q[AlIdx]);
		// printf("The value of mAl is: %.7e\n", mAl); 
		// printf("The value of q[SiIdx] is: %.7e\n", q[SiO2Idx]);
		// printf("The value of mSiO2 is: %.7e\n", mSiO2);
		// printf("The value of glAwcd is: %.7e\n", glAwcd);
		// printf("The value of glr is: %.7e\n", glr);
		// printf("The value of proAwp is: %.7e\n", proAwp);
		// printf("The value of pror is: %.7e\n",  pror); 

		}
		// else if (0 > q[HIdx])
		// {
		// Hdt = mH / abs(q[HIdx]);
		// newdt = Hdt;
		// }
		else 
		{newdt = 0;}
		// printf("The value of q[HIdx] is: %.7e\n", q[HIdx]);
		// printf("The value of mH is: %.7e\n", mH);
		// printf("The value of volFracProtoImogolite is: %.7e\n", volFracProtoImogolite);
		// printf("The value of volume fraction of Glass is: %.7e\n", volFracGlass);
		// printf("The value of volFracProtoImogolite1 is: %.12e\n", volFracProtoImogolite); 
		// printf("The value of q[phiProtoImogoliteIdx1] is: %.7e\n", q[phiProtoImogoliteIdx]);
		
		// if (q[phiProtoImogoliteIdx] > 10)
		// {
        // printf("The value of q[phiProtoImogoliteIdx] is: %.7e\n", q[phiProtoImogoliteIdx]);
		// printf("The value of q[phiGlassIdx] is: %.7e\n", q[phiGlassIdx]);
		// printf("The value of q[AlIdx] is: %.7e\n", q[AlIdx]);
		// printf("The value of mAl is: %.7e\n", mAl); 
		// printf("The value of q[SiIdx] is: %.7e\n", q[SiO2Idx]);
		// printf("The value of mSiO2 is: %.7e\n", mSiO2);
		// printf("The value of q[HIdx] is: %.7e\n", q[HIdx]);
		// printf("The value of mH is: %.7e\n", mH);
	    // // printf("The value of prolr is: %.e\n", pror); 
		// printf("The value of proAw is: %.e\n", proAw); 
		// // printf("The value of proAsw is: %.e\n", proAsw_);
		// // printf("The value of volFracProtoImogolite is: %.7e\n", volFracProtoImogolite);
		// // printf("The value of volVars.solidComponentDensity(pPhaseIdx) is: %.7e\n", volVars.solidComponentDensity(pPhaseIdx));
        // // printf("The value of 1-proOmegaApprox_ is: %.e\n",pow(abs(1-pow(proOmegaApprox_,1/prosigma_)),probeta_));	// Nan	
		// printf("The value of proOmegaApprox_ is: %.7e\n", proOmegaApprox_); //-nan
		// // printf("The value of prorc_ is: %.7e\n", prorc_); //-nan
		// printf("The value of proAw0_ is: %.e\n", proAw0);
        // printf("The value of ratio is: %.7e\n",cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity))* cbrt((volFracProtoImogolite/initialvolFracProtoImogolite)*(volFracProtoImogolite/initialvolFracProtoImogolite)));   // TODO Asw shou		
        // // printf("The value of ratio1 is: %.7e\n",cbrt(((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)*((1-volFracFerrohydrite-volFracGlass-volFracHydroxyapatite-volFracProtoImogolite-volFracSepiolite-volFracBirnessite)/initialPorosity)));	
        // // printf("The value of ratio2 is: %.7e\n",cbrt((volFracProtoImogolite/initialvolFracProtoImogolite)*(volFracProtoImogolite/initialvolFracProtoImogolite)));   // TODO Asw shou		
        // printf("The value of volFracProtoImogolite is: %.7e\n",volFracProtoImogolite);
		// printf("The value of volFracFerrohydrite is: %.7e\n", volFracFerrohydrite);
		// printf("The value of volFracHydroxyapatite is: %.7e\n", volFracHydroxyapatite);
		// printf("The value of volFracBirnessite is: %.7e\n", volFracBirnessite);
		// printf("The value of volFracSepiolite is: %.7e\n", volFracSepiolite);
        // printf("The value of initialvolFracProtoImogolite is: %.7e\n",initialvolFracProtoImogolite);
		// }
		
		
		// if (q[phiFerrohydriteIdx] > 1)
		// {
        // printf("The value of q[phiProtoImogoliteIdx] is: %.7e\n", q[phiProtoImogoliteIdx]);
		// printf("The value of q[phiFerrohydriteIdx] is: %.7e\n", q[phiFerrohydriteIdx]);
		// printf("The value of volume fraction of Ferrohydrite is: %.7e\n", volFracFerrohydrite);
	    // printf("The value of mFe is: %.7e\n", mFe2); 
	    // printf("The value of mH is: %.7e\n", mH);		
		// }
        // printf("The value of q[phiBirnessiteIdx] is: %.7e\n", q[phiBirnessiteIdx]);
        // printf("The value of q[phiHydroxyapatiteIdx] is: %.7e\n", q[phiHydroxyapatiteIdx]);
        // printf("The value of q[phiSepioliteIdx] is: %.7e\n", q[phiSepioliteIdx]);
		// printf("The value of volFracBirnessite is: %.7e\n", volFracBirnessite);
		// printf("The value of volFracSepiolite is: %.20e\n", volFracSepiolite);
        // printf("The value of glOmega is: %.7e\n", pow(abs(1-pow(glOmegaApprox_,1/glsigma_)),glbeta_));
        // printf("The value of glAsw is: %.7e\n", glAsw);
		// printf("The value of volFracGlass is: %.7e\n", volFracGlass);
        // printf("The value of glOmegaApprox_ is: %.7e\n", glOmegaApprox_);
		// printf("The value of glrdiss is: %.10e\n", glrdiss);
		// printf("The value of glrprec is: %.10e\n", glrprec);     
		// printf("The value of calrdiss is: %.10e\n", calrdiss);
		// printf("The value of calrprec is: %.10e\n", calrprec); 
        // printf("The value of q[CaIdx] is: %.10f\n", q[CaIdx]);  
        // printf("The value of mCa is: %.7e\n", mCa);		
		
		// if (- glrprec + glrdiss < 0.0)
        // { printf("The value of - glrprec + glrdiss is: %.7e\n", - glrprec + glrdiss);}
	
        // Close the file
        // fclose(file);		
		// }
        // q[Fe2Idx] += - frprec + frdiss; 
        // q[UreaIdx] += 0;
        // q[TNHIdx] += 0;
        // q[UreaseIdx] += 0;
		// // q[UreaIdx] += - rurea;
        // // q[TNHIdx] += 2 * rurea;
        // // q[UreaseIdx] += (- ra_urease + rd_urease - ria_urease)/FluidSystem::molarMass(UreaseIdx);
//      // q[JBMIdx] += (- ra_iJBM + rd_iJBM)/FluidSystem::molarMass(JBMIdx);
        // q[phiCalciteIdx] += + rprec - rdiss;//added by du. key
		// q[phiCalciteIdx] += 0;
		// q[phiFerrohydriteIdx] += 0;
        // q[phiFerrohydriteIdx] += + frprec - frdiss; // #added by du. key
        // q[phiImmUreaseIdx] += (ra_urease - rd_urease - ria_immUrease)/SolidSystem::molarMass(uPhaseIdx);
		// printf("The value of q1 is: %.10f\n", q[0]); 
		// printf("The value of q is: %.10f\n", q[1]); 
		// printf("The value of q is: %.10f\n", q[2]); 
		// printf("The value of q is: %.10f\n", q[3]); 
		// printf("The value of q is: %.10f\n", q[4]); 
		// printf("The value of q is: %.10f\n", q[5]); 
		// printf("The value of q is: %.10f\n", q[6]); 
		// printf("The value of q is: %.10f\n", q[7]); 
		// printf("The value of q is: %.10f\n", q[8]); 
		// printf("The value of q is: %.10f\n", q[9]); 
		// printf("The value of q is: %.10f\n", q[10]); 
		// printf("The value of q is: %.10f\n", q[11]); 
		// printf("The value of q is: %.10f\n", q[12]); 
		// printf("The value of q is: %.10f\n", q[13]); 
		// printf("The value of q is: %.10f\n", q[14]); 
		// printf("The value of q is: %.10f\n", q[15]); 
		// printf("The value of q is: %.10f\n", q[16]); 
		// printf("The value of q is: %.10f\n", q[17]); 
		// printf("The value of q is: %.10f\n", q[18]); 
        
		//return q;
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

    //Function solves electro neutrality equation f and derivative df/dH for H with constant vary CO2
    void Cl_CO2(Scalar activityCl)
    {
        // printf("reached H_CO2"); reached
        cl_ = activityCl;
		// printf("The value of co2aq_ is: %.8e\n", co2aq_);
        oh_ = kw_/h_;
		// printf("The value of oh_ is: %.8e\n", oh_);
		// co2aq_ = ctot_;
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + (1/kh_) + 1);
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + 1);
        hco3_ = k1_*co2aq_/h_;
		// printf("The value of co2aq is: %.8e\n", co2aq_);
		// printf("The value of k1_ is: %.e\n", k1_);
//         co3_ = k1_*k2_*co2_/pow(h_, 2.);
        co3_ = k2_*hco3_/h_;
		// fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(oh_, 4)), 0.25);
		
        // nh4_ = totalnh_/(1+ka_/h_);
        // co2g_ = co2aq_/kh_;
		// printf("The value of co2g_ is: %.e\n", co2g_);
		// printf("The value of co2aq_ is: %.e\n", co2aq_);
		// printf("The value of hco3_ is: %.8e\n", hco3_);
		// printf("The value of co3_ is: %.8e\n", co3_);
        //Solve the function
        // Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
        Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = cl_ + eps*cl_; // x + dx
        Scalar xLeft = cl_ - eps*cl_; // x - dx
        // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
        // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft - cl_ ; //  f(x-dx)
        Scalar fRight = na_ + h_ + 2*ca_ + 2*fe2_ - oh_ -  hco3_ - 2*co3_ - xRight +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = na_ + h_ + 2*ca_ + 2*fe2_ - oh_ -  hco3_ - 2*co3_ - xLeft +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xLeft); //  f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/cl_; // {f(x+dx) - f(x-dx)}/2dx


        fdf_[0] = f;
        fdf_[1] = df;
     }
	 
    //Function solves electro neutrality equation f and derivative df/dH for H with constant vary CO2
    void H_varyCO2(Scalar activityCO2aq)
    {
        // printf("reached H_CO2"); reached
        co2aq_ = activityCO2aq;
		// printf("The value of co2aq_ is: %.8e\n", co2aq_);
        oh_ = kw_/h_;
		// printf("The value of oh_ is: %.8e\n", oh_);
		// co2aq_ = ctot_;
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + (1/kh_) + 1);
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + 1);
        hco3_ = k1_*co2aq_/h_;
		// printf("The value of co2aq is: %.8e\n", co2aq_);
		// printf("The value of k1_ is: %.e\n", k1_);
//         co3_ = k1_*k2_*co2_/pow(h_, 2.);
        co3_ = k1_*k2_*co2aq_/(h_*h_);
		// fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(oh_, 4)), 0.25);
		
        // nh4_ = totalnh_/(1+ka_/h_);
        // co2g_ = co2aq_/kh_;
		// printf("The value of co2g_ is: %.e\n", co2g_);
		// printf("The value of co2aq_ is: %.e\n", co2aq_);
		// printf("The value of hco3_ is: %.8e\n", hco3_);
		// printf("The value of co3_ is: %.8e\n", co3_);
        //Solve the function
        // Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
        Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = co2aq_ + eps*co2aq_; // x + dx
        Scalar xLeft = co2aq_ - eps*co2aq_; // x - dx
        // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
        // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft - cl_ ; //  f(x-dx)
        Scalar fRight = na_ + h_ + 2*ca_ + 2*fe2_ - kw_/h_ - k1_*xRight/h_ - 2*k1_*k2_*xRight/(h_*h_) - cl_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = na_ + h_ + 2*ca_ + 2*fe2_ - kw_/h_ -  k1_*xLeft/h_ - 2*k1_*k2_*xLeft/(h_*h_) - cl_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xLeft); //  f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/co2aq_; // {f(x+dx) - f(x-dx)}/2dx


        fdf_[0] = f;
        fdf_[1] = df;
     }

    void H_CO3(Scalar activityCO3)
    {
        // printf("reached H_CO2"); reached
        hco3_ = activityCO3; //co2aq_ * k1_ / h_;
		// printf("The value of co2aq_ is: %.8e\n", co2aq_);
        oh_ = kw_/h_;
		// printf("The value of oh_ is: %.8e\n", oh_);
		// co2aq_ = ctot_;
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + (1/kh_) + 1);
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + 1);
        // hco3_ = k1_*co2aq_/h_;
		// h+co3=hco3
		co3_ = hco3_*k2_/h_;
		// printf("The value of co2aq is: %.8e\n", co2aq_);
		// printf("The value of k1_ is: %.e\n", k1_);
//         co3_ = k1_*k2_*co2_/pow(h_, 2.);
        // co3_ = k2_*hco3_/(h_);
		
		// fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(oh_, 4)), 0.25);
		
        // nh4_ = totalnh_/(1+ka_/h_);
        // co2g_ = co2aq_/kh_;
		// printf("The value of co2g_ is: %.e\n", co2g_);
		// printf("The value of co2aq_ is: %.e\n", co2aq_);
		// printf("The value of hco3_ is: %.8e\n", hco3_);
		// printf("The value of co3_ is: %.8e\n", co3_);
        //Solve the function
        // Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
        Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = hco3_ + eps*hco3_; // x + dx
        Scalar xLeft = hco3_ - eps*hco3_; // x - dx
        // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
        // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft - cl_ ; //  f(x-dx)
        Scalar fRight = na_ + h_ + 2*ca_ + 2*fe2_ - kw_/h_ - xRight - 2*xRight*k2_/h_  - cl_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = na_ + h_ + 2*ca_ + 2*fe2_ - kw_/h_ -  xLeft - 2*xLeft*k2_/h_ - cl_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xLeft); //  f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/hco3_; // {f(x+dx) - f(x-dx)}/2dx


        fdf_[0] = f;
        fdf_[1] = df;
     }

	  //Function solves electro neutrality equation f and derivative df/dH for H with constant CO2
     void H_CO2(Scalar activityH)
     {
         // printf("reached H_CO2"); reached
         h_ = activityH;
	 	// printf("The value of h_ is: %.8e\n", h_);
         oh_ = kw_/h_;
	 	// printf("The value of oh_ is: %.8e\n", oh_);
	 	// co2aq_ = ctot_;
	 	// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + (1/kh_) + 1);
	 	// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + 1);
         hco3_ = k1_*co2aq_/h_;
	 	// printf("The value of co2aq is: %.8e\n", co2aq_);
	 	// printf("The value of k1_ is: %.e\n", k1_);
//          co3_ = k1_*k2_*co2_/pow(h_, 2.);
         co3_ = k1_*k2_*co2aq_/(h_*h_);
	 	// fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(oh_, 4)), 0.25);
	 	
         // nh4_ = totalnh_/(1+ka_/h_);
         // co2g_ = co2aq_/kh_;
	 	// printf("The value of co2g_ is: %.e\n", co2g_);
	 	// printf("The value of co2aq_ is: %.e\n", co2aq_);
	 	// printf("The value of hco3_ is: %.8e\n", hco3_);
	 	// printf("The value of co3_ is: %.8e\n", co3_);
         //Solve the function
         // Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
         Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
         //Solve the derivative df/d(activityH)
         Scalar eps = 1e-8;
         Scalar xRight = h_ + eps*h_; // x + dx
         Scalar xLeft = h_ - eps*h_; // x - dx
         // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
         // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft - cl_ ; //  f(x-dx)
         Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - k1_*co2aq_/xRight - 2*k1_*k2_*co2aq_/(xRight*xRight) - cl_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
         Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft -  k1_*co2aq_/xLeft - 2*k1_*k2_*co2aq_/(xLeft*xLeft) - cl_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xLeft); //  f(x-dx)
         Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx
	 
	 
         fdf_[0] = f;
         fdf_[1] = df;
      }
	 
    // need to changed
	// void H_CO2(Scalar activityH)
    // {
    //     // printf("reached H_CO2"); reached
    //     h_ = activityH;
    //     oh_ = kw_/h_;
    //     hco3_ = k1_*co2aqonly_/h_;
    //     co3_ = k1_*k2_*co2aqonly_/(h_*h_);
    //     Scalar f = honly_ - oh_ - 2*co3_ - hco3_ - h_;
    //     //Solve the derivative df/d(activityH)
    //     Scalar eps = 1e-8;
    //     Scalar xRight = honly_ + eps*honly_; // x + dx
    //     Scalar xLeft = honly_ - eps*honly_; // x - dx        Scalar fRight = xRight - kw_/xRight - co2aq_/(xRight/k1_ + 1 + k2_/xRight ) - 2*co2aq_/((xRight*xRight)/k1_/k2_ + xRight/k2_ + 1 ) - h_;// + totalnh_/(1+ka_/xRight); // f(x+dx)
    //     Scalar fLeft = xLeft - kw_/xLeft - co2aq_/(xLeft/k1_ + 1 + k2_/xLeft) - 2*co2aq_/((xLeft*xLeft)/k1_/k2_ + xLeft/k2_ + 1 ) - h_;// + totalnh_/(1+ka_/xLeft); //  f(x-dx)
    //     Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx
	// 
	// 
    //     fdf_[0] = f;
    //     fdf_[1] = df;
    //  }


    void co2aq_CO2(Scalar activityco2aq)
    {
        // printf("reached H_CO2"); reached
        co2aq_ = activityco2aq;
		// printf("The value of h_ is: %.8e\n", h_);
        oh_ = kw_/h_;
		// printf("The value of oh_ is: %.8e\n", oh_);
		// co2aq_ = ctot_;
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + (1/kh_) + 1);
		// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + 1);
        hco3_ = k1_*co2aq_/(h_);
		// printf("The value of co2aq is: %.8e\n", co2aq_);
		// printf("The value of k1_ is: %.e\n", k1_);
//         co3_ = k1_*k2_*co2_/pow(h_, 2.);
        co3_ = k1_*k2_*co2aq_/(h_*h_);
		fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(oh_, 4)), 0.25);
		
        // nh4_ = totalnh_/(1+ka_/h_);
        // co2g_ = co2aq_/kh_;
		// printf("The value of co2g_ is: %.e\n", co2g_);
		// printf("The value of co2aq_ is: %.e\n", co2aq_);
		// printf("The value of hco3_ is: %.8e\n", hco3_);
		// printf("The value of co3_ is: %.8e\n", co3_);
        //Solve the function
        // Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
        Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ + 3.*fe3_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = co2aq_ + eps*co2aq_; // x + dx
        Scalar xLeft = co2aq_ - eps*co2aq_; // x - dx
        // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
        // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft - cl_ ; //  f(x-dx)
        Scalar fRight = na_ + h_ + 2*ca_ + 2*fe2_ - oh_ - k1_*xRight/h_ - 2*k1_*k2_*xRight/(h_*h_) - cl_ + 3*fe3_+3.*al_ +2.*mg_ +k_ + 2.*mn_ +4.*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = na_ + h_ + 2*ca_ + 2*fe2_ - oh_ -  k1_*xLeft/h_ - 2*k1_*k2_*xLeft/(h_*h_) - cl_ +  3*fe3_+3.*al_ +2.*mg_ +k_ + 2.*mn_ +4.*tioh4_-2.*hpo4_; // + totalnh_/(1+ka_/xLeft); //  f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/co2aq_; // {f(x+dx) - f(x-dx)}/2dx


        fdf_[0] = f;
        fdf_[1] = df;
     }

    // void Cl_CO2(Scalar activityCl)
    // {
    //     // printf("reached H_CO2"); reached
    //     cl_ = activityCl;
	// 	// printf("The value of h_ is: %.8e\n", h_);
    //     oh_ = kw_/h_;
	// 	// printf("The value of oh_ is: %.8e\n", oh_);
	// 	// co2aq_ = ctot_;
	// 	// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + (1/kh_) + 1);
	// 	// co2aq_ = ctot_ / (k1_/h_ + k1_*k2_/(h_*h_) + 1);
    //     hco3_ = k1_*co2aq_/h_;
	// 	// printf("The value of co2aq is: %.8e\n", co2aq_);
	// 	// printf("The value of k1_ is: %.e\n", k1_);
//  //        co3_ = k1_*k2_*co2_/pow(h_, 2.);
    //     co3_ = k1_*k2_*co2aq_/(h_*h_);
	// 	fe3_ = pow((k3_ * pow(fe2_, 4) * o2_ / pow(oh_, 4)), 0.25);
	// 	
    //     // nh4_ = totalnh_/(1+ka_/h_);
    //     // co2g_ = co2aq_/kh_;
	// 	// printf("The value of co2g_ is: %.e\n", co2g_);
	// 	// printf("The value of co2aq_ is: %.e\n", co2aq_);
	// 	// printf("The value of hco3_ is: %.8e\n", hco3_);
	// 	// printf("The value of co3_ is: %.8e\n", co3_);
    //     //Solve the function
    //     // Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ + nh4_+ 2.*fe2_;
    //     Scalar f = na_ + h_ + 2*ca_ - oh_ -  hco3_ - 2*co3_ - cl_ + 2.*fe2_ + 3.*fe3_ +3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_ -2.*hpo4_ ;
    //     //Solve the derivative df/d(activityH)
    //     Scalar eps = 1e-8;
    //     Scalar xRight = h_ + eps*h_; // x + dx
    //     Scalar xLeft = h_ - eps*h_; // x - dx
    //     // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
    //     // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_ - kw_/xLeft - cl_ ; //  f(x-dx)
    //     Scalar fRight = na_ - xRight + 2*ca_ + 2*fe2_ - kw_/h_ - k1_*co2aq_/h_ - 2*k1_*k2_*co2aq_/(h_*h_) + h_ + 3*pow((k3_ * pow(fe2_, 4) * o2_ / pow((kw_/h_), 4)), 0.25+3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_); // + totalnh_/(1+ka_/xRight); // f(x+dx)
    //     Scalar fLeft = na_ - xLeft + 2*ca_ + 2*fe2_ - kw_/h_ -  k1_*co2aq_/h_ - 2*k1_*k2_*co2aq_/(h_*h_) + h_ +  3*pow((k3_ * pow(fe2_, 4) * o2_ / pow((kw_/h_), 4)), 0.25+3.*al_ +2.*mg_ +k_ + 2.*mn_ +4*tioh4_-2.*hpo4_); // + totalnh_/(1+ka_/xLeft); //  f(x-dx)
    //     Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx
	// 
	// 
    //     fdf_[0] = f;
    //     fdf_[1] = df;
    //  }
    // void H_Ctot(Scalar activityH)
    //     {
	// 
    //     h_ = activityH;
    //     oh_ = kw_/h_;
    //     // hco3_ = ctot_/(h_/k1_ + 1 + k2_/h_ + h_/(kh_ * k1_));
    //     // co3_ = ctot_/((h_*h_)/k1_/k2_ + h_/k2_ + 1 +(h_*h_)/(kh_ * k1_ *k2_));
    //     // co2aq_ = (ctot_-co3_-hco3_)/(1+1/kh_);//h_*hco3_/k1_;
    //     hco3_ = ctot_/(h_/k1_ + 1 + k2_/h_ );
    //     co3_ = ctot_/((h_*h_)/k1_/k2_ + h_/k2_ + 1 );
	// 	co2aq_ = ctot_-co3_-hco3_;
	// 	// co2g_ = co2aq_ / kh_;
    //     // nh4_ = totalnh_/(1+ka_/h_);
	// 
    //     //Solve the function
    //     Scalar f = na_ + h_ + 2*ca_ - oh_ - hco3_ - 2*co3_ - cl_ +  2.*fe2_;
    //     // Scalar f = na_ + h_ + 2*ca_ - oh_ - cl_ + 2.*fe2_;
    //     //Solve the derivative df/d(activityH)
    //     Scalar eps = 1e-8;
    //     Scalar xRight = h_ + eps*h_; // x + dx
    //     Scalar xLeft = h_ - eps*h_; // x - dx
    //     Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - ctot_/(xRight/k1_ + 1 + k2_/xRight) - 2*ctot_/((xRight*xRight)/k1_/k2_ + xRight/k2_ + 1) - cl_;// + totalnh_/(1+ka_/xRight); // f(x+dx)
    //     Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_  - kw_/xLeft - ctot_/(xLeft/k1_ + 1 + k2_/xLeft) - 2*ctot_/((xLeft*xLeft)/k1_/k2_ + xLeft/k2_ + 1) - cl_;// + totalnh_/(1+ka_/xLeft); //  f(x-dx)
    //     // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
    //     // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_  - kw_/xLeft - cl_ ; //  f(x-dx)
    //     Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx
	// 
	// 
    //     fdf_[0] = f;
    //     fdf_[1] = df;
    // }
    void H_Ctot(Scalar activityH)
        {

        honly_ = activityH;
        oh_ = kw_/honly_;
        hco3_ = co2aq_/(honly_/k1_ + 1 + k2_/honly_ );
        co3_ = co2aq_/((honly_*honly_)/k1_/k2_ + honly_/k2_ + 1 );
		co2aqonly_ = co2aq_-co3_-hco3_;
		fe3_ = fe2_ / (1+k11_/(honly_* pow(o2_,0.25)));
		fe2only_ = fe2_ - fe3_;
		
		// co2g_ = co2aq_ / kh_;
        // nh4_ = totalnh_/(1+ka_/h_);

        //Solve the function
        Scalar f = honly_ - oh_ - 2*co3_ - hco3_ + fe3_ - h_;
        // Scalar f = na_ + h_ + 2*ca_ - oh_ - cl_ + 2.*fe2_;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = honly_ + eps*honly_; // x + dx
        Scalar xLeft = honly_ - eps*honly_; // x - dx
        // Scalar f = honly_ - kw_/honly_ - co2aq_/(honly_/k1_ + 1 + k2_/honly_ )_ - 2*co2aq_/((honly_*honly_)/k1_/k2_ + honly_/k2_ + 1 )_ - h_;// + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fRight = xRight - kw_/xRight - co2aq_/(xRight/k1_ + 1 + k2_/xRight ) - 2*co2aq_/((xRight*xRight)/k1_/k2_ + xRight/k2_ + 1 ) + fe2_ / (1+k11_/(xRight * pow(o2_,0.25))) - h_;// + totalnh_/(1+ka_/xRight); // f(x+dx)
        Scalar fLeft = xLeft - kw_/xLeft - co2aq_/(xLeft/k1_ + 1 + k2_/xLeft) - 2*co2aq_/((xLeft*xLeft)/k1_/k2_ + xLeft/k2_ + 1 ) + fe2_ / (1+k11_/(xLeft * pow(o2_,0.25))) - h_;// + totalnh_/(1+ka_/xLeft); //  f(x-dx)
        // Scalar fRight = na_ + xRight + 2*ca_ + 2*fe2_ - kw_/xRight - cl_ ; // f(x+dx)
        // Scalar fLeft = na_ + xLeft + 2*ca_ + 2*fe2_  - kw_/xLeft - cl_ ; //  f(x-dx)
        Scalar df = (fRight - fLeft)/2/eps/honly_; // {f(x+dx) - f(x-dx)}/2dx


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
    Scalar co2aqonly_;
	Scalar co2g_;
    Scalar o2_;
    Scalar hco3_;
    Scalar co3_;
    Scalar oh_;
    Scalar h_;
    Scalar honly_;
    Scalar h2_;
    Scalar ca_;
    Scalar na_;
    Scalar cl_;
    Scalar fe2_;
	Scalar fe2only_;
    Scalar fe3_;
    Scalar sio2_;
    Scalar tioh4_;
    Scalar al_;
    Scalar mn_;
    Scalar mg_;
    Scalar k_;
    Scalar hpo4_;
 	
    // Scalar totalnh_;
    // Scalar nh4_;
    Scalar initH_;
    Scalar initCl_;
    Scalar initCO2aq_;
    Scalar initCO3_;
    Scalar ionicStrength_;
    Scalar ctot_;
    Scalar gammaH_;
    Scalar gammaCO2_;
    Scalar gammaCa_;
    Scalar gammaOH_;
    // Scalar gammaHCO3_;
    // Scalar gammaCO3_;
    // Scalar gammaNH3_;
    // Scalar gammaNH4_;
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
	// Scalar kh_;
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

    // //attachment and detachment parameters
    //     Scalar ka_urease_;
    //     Scalar kd_urease_;
    //     Scalar ka_iJBM_;
    //     Scalar kd_iJBM_;

    // glass parameters
        // Scalar glac_;
        // Scalar glkdiss1_;
        // Scalar glkdiss2_;
        // Scalar glkprec_;
        // Scalar glndiss_;
        // Scalar glnprec_;
        Scalar glAsw_;
        Scalar glrc_;	
        Scalar glbeta_;
        Scalar glsigma_;
        Scalar glp_;
    // // calcite parameters
    //     Scalar ac_;
    //     Scalar kdiss1_;
    //     Scalar kdiss2_;
    //     Scalar kprec_;
    //     Scalar ndiss_;
    //     Scalar nprec_;
    //     Scalar Asw0_;
	// 
	     // Scalar calAsw_;
         // Scalar calrc_;	
         // Scalar calbeta_;
         // Scalar calsigma_;
         // Scalar calp_;

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
		 
    // // ferrohydrite parameters
    //     Scalar fac_;
    //     Scalar fkdiss1_;
    //     Scalar fkdiss2_;
    //     Scalar fkprec_;
    //     Scalar fndiss_;
    //     Scalar fnprec_;
    //     Scalar fAsw0_;

    // // urease parameters
    //     bool useHeatKilledCells_;
    //     bool useJackBeans_;
    //     Scalar kub_;
    //     Scalar kurease_;
    //     Scalar nub_;
    //     Scalar Keu1_;
    //     Scalar Keu2_;
    //     Scalar cia_;
    //     Scalar ciaT_;
    //     Scalar cureaseT_;
//  //        Scalar sorptionCoeffUrease_;
//  //        Scalar maxSorptionUrease_;
//  //        Scalar pmDensity_;
    //     Scalar cu_,cuT_;
    //     Scalar ciaPrec_;
        Scalar pKaFactor_;

public:

    // calcite parameters
       // Scalar glac()    {       return glac_; }
       // Scalar glkdiss1()    {    return glkdiss1_; }
       // Scalar glkdiss2()    {    return glkdiss2_; }
       // Scalar glkprec()    {       return glkprec_; }
       // Scalar glndiss()    {       return glndiss_; }
       // Scalar glnprec()    {       return glnprec_; }
       Scalar glAsw()    {       return glAsw_; }
       Scalar glrc()    {       return glrc_; }
       Scalar glp()    {       return glp_; }
       Scalar glbeta()    {       return glbeta_; }
       Scalar glsigma()    {       return glsigma_; }	   
    // // calcite parameters
    //    Scalar ac()    {       return ac_; }
    //    Scalar kdiss1()    {    return kdiss1_; }
    //    Scalar kdiss2()    {    return kdiss2_; }
    //    Scalar kprec()    {       return kprec_; }
    //    Scalar ndiss()    {       return ndiss_; }
    //    Scalar nprec()    {       return nprec_; }
    //    Scalar Asw0()    {       return Asw0_; }
	// 
    // // ferrohydrite parameters
    //    Scalar fac()    {       return fac_; }
    //    Scalar fkdiss1()    {    return fkdiss1_; }
    //    Scalar fkdiss2()    {    return fkdiss2_; }
    //    Scalar fkprec()    {       return fkprec_; }
    //    Scalar fndiss()    {       return fndiss_; }
    //    Scalar fnprec()    {       return fnprec_; }
    //    Scalar fAsw0()    {       return fAsw0_; }

    // // urease parameters
    //     Scalar kub()    {       return kub_; }
    //     Scalar kurease()    {   return kurease_; }
    //     Scalar nub()    {       return nub_; }
    //     Scalar Keu1()    {       return Keu1_; }
    //     Scalar Keu2()    {       return Keu2_; }
    //     Scalar cia()    {       return cia_; }
    //     Scalar ciaT()    {       return ciaT_; }
    //     Scalar cureaseT()    {       return cureaseT_; }


public:
   // Scalar kprec() const
   // {   return kprec_;}
    // Scalar kub() const
    // {   return kub_;}
    // Scalar kurease() const
    // {   return kurease_;}
    // Scalar nprec() const
    // {   return nprec_;}
    // Scalar Asw0() const
    // {   return Asw0_;}
    // Scalar fnprec() const
    // {   return fnprec_;}
    // Scalar fAsw0() const
    // {   return fAsw0_;}
    // Scalar Keu1() const
    // {   return Keu1_;}
    // Scalar Keu2() const
    // {   return Keu2_;}
    // Scalar cia() const
    // {   return cia_;}
    // Scalar ciaT() const
    // {   return ciaT_;}
    // Scalar cureaseT() const
    // {   return cureaseT_;}

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










