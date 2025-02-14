
#ifndef LEO_CARBONIC_ACID_HH_
#define LEO_CARBONIC_ACID_HH_

#include <dumux/common/exceptions.hh>
#include "leomin.hh"

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
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ThisType = LeoMinCarbonicAcid<TypeTag, CO2Tables, ModelTraits>;
    using H2O = Components::H2O<Scalar>;

    enum
    {

    };

public:

    static const int H2OIdx       = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::H2OIdx);
    static const int N2Idx        = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::N2Idx);
    static const int HIdx        = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::HIdx);
    static const int HtotalIdx        = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::HtotalIdx);
    static const int OHIdx        = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::OHIdx);
    static const int CO2aqIdx       = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::CO2aqIdx );
    static const int CO2aqtotalIdx       = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::CO2aqtotalIdx);
    static const int HCO3Idx      = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::HCO3Idx );
    static const int CO3Idx       = FluidSystem::compIdx(FluidSystem::MultiPhaseFluidSystem::CO3Idx );
	
    static const int numComponents      = FluidSystem::numComponents;
    static const int numSecComponents   = FluidSystem::numSecComponents;
    static const int numTotComponents   = numComponents + numSecComponents;

    typedef Dune::FieldVector<Scalar, 4> Vector;   // Ionic Strength with NH4/totalnh
    typedef Dune::FieldVector<Scalar, 2> SolVector;
	typedef Dune::FieldVector<Scalar, numTotComponents> CompVector;
	
    template <class FluidState>
    void calculateEquilibriumChemistry(const FluidState &fluidState, CompVector &variable, Scalar rhoMolar)
 {
        const VolumeVariables volVars{};
        gammaCO2_ = 1.0;
        Scalar h2o_ = moleFracToMolarity(variable[H2OIdx], rhoMolar);
        pressure_ = fluidState.pressure(0);
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

	co2aq_ = moleFracToMolarity(variable[CO2aqIdx],rhoMolar)-moleFracToMolarity(variable[CO3Idx],rhoMolar)-moleFracToMolarity(variable[HCO3Idx],rhoMolar);		
        h_ = moleFracToMolarity(variable[HIdx],rhoMolar)+moleFracToMolarity(variable[OHIdx],rhoMolar)+moleFracToMolarity(variable[CO3Idx],rhoMolar)+moleFracToMolarity(variable[HCO3Idx],rhoMolar);	
	// co2aq_ = moleFracToMolarity(variable[CO2aqIdx],rhoMolar);
        // h_ = moleFracToMolarity(variable[HIdx],rhoMolar);
		
	initH_ = 8.07096e-12;//h_; //Initial guess
        Scalar activityH = initH_;		
	k1_ = const1(pressure_, temperature_);
        k2_ = const2(pressure_, temperature_);
        kw_ = constW(pressure_, temperature_);
		
        Scalar tolAbs = 1e-20;
        Scalar tolRel = 1e-20;
        int maxIter = 500;
		
	if(newton1D(activityH, &ThisType::H_Conly, tolAbs, tolRel, maxIter) == false) //Alex' Newton
         {
             initH_ = 8.07096e-12; //8.07096e-13 ;//h_;
			 activityH = initH_;		
             Scalar a0 = 0.0;
             Scalar b0 = 1e-1;
             Scalar tol = 1e-15; //1e-20 makes it unable to run, the H becomes 0
			if(bisection1D(activityH, &ThisType::H_Conly, a0, b0, tol) == false) //Alex' bisection
		
             {
                 DUNE_THROW(Dune::InvalidStateException, "in Chemistry: Bisection did not converge!" );
             }
         }

	H_Conly(activityH);
	htotal_ = h_ - oh_ - co3_ - hco3_ ;
	co2aqtotal_ = co2aq_ + co3_ + hco3_;
	 
        Scalar totalMolarity = h2o_ + n2_ + co2aqtotal_ + htotal_; // du do not know
   
	variable[CO2aqIdx] = co2aq_/totalMolarity;
	variable[CO2aqtotalIdx] = co2aqtotal_/totalMolarity;
        variable[HCO3Idx] = hco3_/totalMolarity;
        variable[CO3Idx] = co3_/totalMolarity;
        variable[OHIdx] = oh_/totalMolarity;
        variable[HIdx] = h_/totalMolarity;
        variable[HtotalIdx] = htotal_/totalMolarity;

 }

    //Return equlibrium constant for chemical equation:
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
    // H2O <--> H + OH
    static Scalar constW(const Scalar pw, const Scalar T)
    {
        return 1e-14;
    }


    // added by du
    static Scalar moleFracToMolarity(Scalar moleFracX, Scalar molarDensity)
    {
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

   Scalar pH(const VolumeVariables &volVars)
   {
      Scalar mH = moleFracToMolarity(volVars.moleFraction(0,HIdx));  //[mol_H/kg_H2O]
      Scalar pH = -log10(mH);
         return pH;
   }


private:
    //Newton Solver which returns true if convergence is reached and false if not.
    // x(i+1) = x(i) - f(i)/df(i) = x(i) + h(i)

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
        int maxIter = 500;
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

    void H_Conly(Scalar activityH)
        {

        h_ = activityH;
        oh_ = kw_/h_;
		hco3_ = k1_*co2aq_/(h_);
        co3_ = k1_*k2_*co2aq_/(h_*h_);
		
		
        //Solve the function
        Scalar f = h_ - oh_ -  hco3_ - 2*co3_;
        //Solve the derivative df/d(activityH)
        Scalar eps = 1e-8;
        Scalar xRight = h_ + eps*h_; // x + dx
        Scalar xLeft = h_ - eps*h_; // x - dx
        Scalar fRight =  xRight  - kw_/xRight - k1_*co2aq_/xRight - 2*k1_*k2_*co2aq_/(xRight*xRight); // f(x+dx)
        Scalar fLeft =  xLeft - kw_/xLeft - k1_*co2aq_/xLeft - 2*k1_*k2_*co2aq_/(xLeft*xLeft); // f(x+dx)
        Scalar df = (fRight - fLeft)/2/eps/h_; // {f(x+dx) - f(x-dx)}/2dx

        fdf_[0] = f;
        fdf_[1] = df;
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
    Scalar hco3_;
    Scalar co3_;
    Scalar oh_;
    Scalar h_;
    Scalar htotal_;
 	
    Scalar initH_;
    Scalar initCO2_;
    Scalar initCl_;
    Scalar initCO2aq_;
    Scalar initCO3_;

    Scalar ctot_;
    Scalar gammaH_;
    Scalar gammaCO2_;
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
    Scalar kw_;
    Scalar ka_;
    Scalar apparentk1_;
    Scalar apparentk2_;
    Scalar apparentka_;
	
    bool newtonOrBisection_;

    static constexpr Scalar KpHb_ = 0;//9.14e-8;//[mol/kgH2O] Kim et al. 2000 //Not implemented by Anozie!!

};

} // end namespace

#endif










