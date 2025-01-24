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
 * \ingroup ConstraintSolvers
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 */
#ifndef DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH
#define DUMUX_MISCIBLE_MULTIPHASE_COMPOSITION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/valgrind.hh>

namespace Dumux {
/*!
 * \ingroup ConstraintSolvers
 * \brief Computes the composition of all phases of a N-phase,
 *        N-component fluid system assuming that all N phases are
 *        present
 *
 * The constraint solver assumes the following quantities to be set:
 *
 * - temperatures of *all* phases
 * - saturations of *all* phases
 * - pressures of *all* phases
 *
 * It also assumes that the mole/mass fractions of all phases sum up
 * to 1. After calling the solve() method the following quantities
 * are calculated in addition:
 *
 * - temperature of *all* phases
 * - density, molar density, molar volume of *all* phases
 * - composition in mole and mass fractions and molarities of *all* phases
 * - mean molar masses of *all* phases
 * - fugacity coefficients of *all* components in *all* phases
 */
template <class Scalar, class FluidSystem>
class MiscibleMultiPhaseComposition
{
    static constexpr int numPhases = FluidSystem::numPhases;
    static constexpr int numComponents = FluidSystem::numComponents;
    static const int numMajorComponents = FluidSystem::numPhases;

public:
    /*!
     * \brief @copybrief Dumux::MiscibleMultiPhaseComposition
     *
     * This function additionally considers a lowering of the saturation vapor pressure
     * of the wetting phase by the Kelvin equation:
     * \f[
     * p^\textrm{w}_\textrm{sat,Kelvin}
     * = p^\textrm{w}_\textrm{sat}
     *   \exp \left( -\frac{p_\textrm{c}}{\varrho_\textrm{w} R_\textrm{w} T} \right)
     * \f]
     *
     * \param fluidState A container with the current (physical) state of the fluid
     * \param paramCache A container for iterative calculation of fluid composition
     * \param knownPhaseIdx The index of the phase with known properties
     */
    template <class FluidState, class ParameterCache>
    static void solve(FluidState &fluidState,
                      ParameterCache &paramCache,
                      int knownPhaseIdx = 0) //, int ii = 0) //, double co2aqn = 0.00039) //, double ratio = 1)
    {
#ifndef NDEBUG
        // currently this solver can only handle fluid systems which
        // assume ideal mixtures of all fluids. TODO: relax this
        // (requires solving a non-linear system of equations, i.e. using
        // Newton method.)
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            assert(FluidSystem::isIdealMixture(phaseIdx));

        }
#endif

        //get the known mole fractions from the fluidState
        //in a 2pnc system the n>2 mole fractions are primary variables and are already set in the fluidstate
        Dune::FieldVector<Scalar, numComponents-numMajorComponents> xKnown(0.0);
        for (int knownCompIdx = 0; knownCompIdx < numComponents-numMajorComponents; ++knownCompIdx)
        {
            xKnown[knownCompIdx] = fluidState.moleFraction(knownPhaseIdx, knownCompIdx + numMajorComponents);
			// printf("reached here");
        }

        // compute all fugacity coefficients
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            paramCache.updatePhase(fluidState, phaseIdx);

		    // // added by Du
            // Scalar CO2aqIdx = FluidSystem::CO2aqIdx;
			// Scalar CO2aqonlyIdx = FluidSystem::CO2aqonlyIdx;
			// // double ratio = 1.0;
		    // // printf("The value of moleFrac[CO2aqonlyIdx1] is: %.8e\n", fluidState.moleFraction(0, CO2aqonlyIdx));
		    // // printf("The value of moleFrac[CO2aqIdx] is: %.8e\n", fluidState.moleFraction(0,CO2aqIdx));
		    // double ratio = fluidState.moleFraction(0, CO2aqonlyIdx) / fluidState.moleFraction(0,CO2aqIdx);
            // // printf("ratio1: %.5e\n",ratio);
		    // if (ratio < 0.0001)
		    // {
		    // 	ratio = 2.38640000e-07 / 2.65687305e-06;
		    // 	// printf("ratio1: %.5e\n",ratio);
		    // }
            // ratio = 1.0;
			// Scalar fugCoeff = 0.0;
            // since we assume ideal mixtures, the fugacity
            // coefficients of the components cannot depend on
            // composition, i.e. the parameters in the cache are valid
			
			// TRUE
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
				// if (compIdx == CO2aqIdx) {
				// fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, CO2aqIdx) * ratio;
				// }
                // else
				// {
				Scalar fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
				// }
                fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
			}

            // // ADDED BY DU, Dec 5
			// Scalar CO2aqIdx = FluidSystem::CO2aqIdx;
			// Scalar O2Idx = FluidSystem::O2Idx;
			// Scalar fugCoeff = 0.0;
	        // for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
			// if (compIdx == CO2aqIdx) {
			// fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, CO2aqIdx) * 1E-3;
			// fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
			// }
			// else if (compIdx == O2Idx) {
			// fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, O2Idx) * 1E-3;
			// fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
			// }		
            // else
			// {
			// fugCoeff = FluidSystem::fugacityCoefficient(fluidState, paramCache, phaseIdx, compIdx);
			// fluidState.setFugacityCoefficient(phaseIdx, compIdx, fugCoeff);
			// }
            // 
            // }		
			
		//}

			// printf("fugCoeff1: %.5e\n",fluidState.fugacityCoefficient(0, CO2aqIdx));

            // fluidState.setFugacityCoefficient(phaseIdx, CO2aqIdx, fugCoeffCO2);	
            // printf("reach here");
            // printf("ratio2: %.5e\n",ratio);	
			// printf("fugCoeff2: %.5e\n",fluidState.fugacityCoefficient(0, CO2aqIdx));
        }


        // create the linear system of equations which defines the
        // mole fractions
        Dune::FieldMatrix<Scalar, numComponents*numPhases, numComponents*numPhases> M(0.0);
        Dune::FieldVector<Scalar, numComponents*numPhases> x(0.0);
        Dune::FieldVector<Scalar, numComponents*numPhases> b(0.0);

        // assemble the equations expressing the assumption that the
        // sum of all mole fractions in each phase must be 1
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            int rowIdx = numComponents*(numPhases - 1) + phaseIdx;
            b[rowIdx] = 1.0;

            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int colIdx = phaseIdx*numComponents + compIdx;

                M[rowIdx][colIdx] = 1.0;
            }
        }

        // set the additional equations for the numComponents-numMajorComponents
        // this is only relevant if numphases != numcomponents e.g. in a 2pnc system
        // Components, of which the mole fractions are known, set to molefraction(knownCompIdx)=xKnown
        for(int knownCompIdx = 0; knownCompIdx < numComponents-numMajorComponents; ++knownCompIdx)
        {
            int rowIdx = numComponents + numPhases + knownCompIdx;
            int colIdx = knownPhaseIdx*numComponents + knownCompIdx + numMajorComponents;
            M[rowIdx][colIdx] = 1.0;
            b[rowIdx] = xKnown[knownCompIdx];
        }

        // assemble the equations expressing the fact that the
        // fugacities of each component are equal in all phases
        for (int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
		   //// added by Du
           //Scalar CO2aqIdx = FluidSystem::CO2aqIdx;
			//Scalar CO2aqonlyIdx = FluidSystem::CO2aqonlyIdx;
			//// double ratio = 1.0;
		   //// printf("The value of moleFrac[CO2aqonlyIdx1] is: %.8e\n", fluidState.moleFraction(0, CO2aqonlyIdx));
		   //// printf("The value of moleFrac[CO2aqIdx] is: %.8e\n", fluidState.moleFraction(0,CO2aqIdx));
		   //double ratio = fluidState.moleFraction(0, CO2aqonlyIdx) / fluidState.moleFraction(0,CO2aqIdx);
           //// printf("ratio1: %.5e\n",ratio);
		   //if (ratio < 0.0001)
		   //{
		   //	ratio = 2.38640000e-07 / 2.65687305e-06;
		   //	// printf("ratio1: %.5e\n",ratio);
		   //}

            int col1Idx = compIdx;
            const auto entryPhase0 = fluidState.fugacityCoefficient(0, compIdx)*fluidState.pressure(0);// *ratio;

            for (unsigned int phaseIdx = 1; phaseIdx < numPhases; ++phaseIdx)
            {
                int rowIdx = (phaseIdx - 1)*numComponents + compIdx;
                int col2Idx = phaseIdx*numComponents + compIdx;
                M[rowIdx][col1Idx] = entryPhase0;
                M[rowIdx][col2Idx] = -fluidState.fugacityCoefficient(phaseIdx, compIdx)*101325; // fluidState.pressure(phaseIdx);
            }
        }

        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            // Multiply row of main component (Raoult's Law) with 10e-5 (order of magn. of pressure)
            if (compIdx < numMajorComponents)
                M[compIdx] *= 10e-5;

            // Multiply row of sec. components (Henry's Law) with 10e-9 (order of magn. of Henry constant)
            else
                M[compIdx] *= 10e-9;

        }
		
		// added by du
        // Index of the equation/variable you want to exclude from being recalculated
        Scalar CO2aqIdx = FluidSystem::CO2aqIdx;
        Scalar CO2aqonlyIdx = FluidSystem::CO2aqonlyIdx;
		// if (fluidState.moleFraction(0, CO2aqonlyIdx) == 0){
		// // co2aqn = 2.38640000e-07 * fluidState.fugacityCoefficient(0, CO2aqIdx);
		// fluidState.setMoleFractionSecComp(0, CO2aqonlyIdx, 2.38640000e-07);
		// // fluidState.moleFraction(0, CO2aqonlyIdx) = 2.38640000e-07;	
		// }
        // Scalar co2aqn = 2.38640000e-07 * fluidState.fugacityCoefficient(0, CO2aqIdx);
        Scalar co2aqn = fluidState.moleFraction(0, CO2aqonlyIdx) * fluidState.fugacityCoefficient(0, CO2aqIdx);

        // Scalar co2aqn = variable[CO2aqonlyIdx] * fluidState.fugacityCoefficient(0, CO2aqIdx);
		// printf("fluidState.moleFraction(0, CO2aqonlyIdx): %.8e\n", variable[CO2aqonlyIdx]);
		
		// if (co2aqn == 0){
		// // co2aqn = 2.38640000e-07 * fluidState.fugacityCoefficient(0, CO2aqIdx);
		// co2aqn = fluidState.moleFraction(0, CO2aqonlyIdx) * fluidState.fugacityCoefficient(0, CO2aqIdx);
		// // printf("reach here");
		// }
		
		// printf("co2aqn: %.8e\n", co2aqn);
		
		Scalar fixedIdx = CO2aqIdx;  // Set the index you want to exclude
        Scalar fixedValue = co2aqn;  // Set the value for this index
        for (int j = 0; j < numComponents * numPhases; ++j) {
        // printf("The value of co2aqn is: %.8e\n", co2aqn);
        // Before solving, modify the row corresponding to fixedIdx
        M[fixedIdx][j] = (j == fixedIdx+numComponents) ? 1.0 : 0.0;
        // b[fixedIdx] = fixedValue * fluidState.pressure(0) / fluidState.pressure(1) /fluidState.fugacityCoefficient(1, CO2aqIdx); //weird here
        b[fixedIdx] = fixedValue * fluidState.pressure(0) / 101325 /fluidState.fugacityCoefficient(1, CO2aqIdx);
        }
                // Set b[fixedIdx] to the fixed value
        // b[fixedIdx] = fixedValue * 1e-8; //weird here

        // if (fluidState.moleFraction(0, CO2aqonlyIdx) != 0 && fluidState.moleFraction(0, CO2aqonlyIdx) != 2.38640000e-07){
		// printf("fluidState.moleFraction(0, CO2aqonlyIdx)4: %.8e\n", fluidState.moleFraction(0, CO2aqonlyIdx));
		// printf("fluidState.fugacityCoefficient(0, CO2aqIdx)4: %.8e\n", fluidState.fugacityCoefficient(0, CO2aqIdx));
		// printf("fluidState.fugacityCoefficient(1, CO2aqIdx)4: %.8e\n", fluidState.fugacityCoefficient(1, CO2aqIdx));
		// printf("fluidState.pressure(0)4: %.8e\n", fluidState.pressure(0));
		// printf("fluidState.pressure(1)4: %.8e\n", fluidState.pressure(1));	
        // printf("b[fixedIdx]: %.8e\n", b[fixedIdx]);		
		// }        
        // for (int i = 0; i < numComponents*numPhases; i++) {
        //     for (int j = 0; j < numComponents*numPhases; j++) {
        //         printf("%.2f ", M[i][j]);  // Adjust the format specifier based on your data type (e.g., %d for integers).
        //     }
        //     printf("\n");  // Newline after each row.
        // }
		// 
        // for (int i = 0; i < numComponents*numPhases; i++) {
        //         printf("%.2e ", b[i]);  // Adjust the format specifier based on your data type (e.g., %d for integers).
        //     }
        // printf("\n");  // Newline after each row.

		
		// printf ("reach here");
        // preconditioning of M to reduce condition number

        


        // solve for all mole fractions
        try { M.solve(x, b); }
        catch (Dune::FMatrixError & e) {
            DUNE_THROW(NumericalProblem,
                    "Matrix for composition of phases could not be solved. \n"
                    "Throwing NumericalProblem for trying with smaller timestep.");
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }

		
        // for (int i = 0; i < numComponents*numPhases; i++) {
        //         printf("%.8e ", x[i]);  // Adjust the format specifier based on your data type (e.g., %d for integers).
        //     }
        // printf("\n");  // Newline after each row.
        // set all mole fractions and the the additional quantities in
        // the fluid state
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                int rowIdx = phaseIdx*numComponents + compIdx;
                fluidState.setMoleFraction(phaseIdx, compIdx, x[rowIdx]);
            }
            paramCache.updateComposition(fluidState, phaseIdx);

            Scalar value = FluidSystem::density(fluidState, paramCache, phaseIdx);
            fluidState.setDensity(phaseIdx, value);

            value = FluidSystem::molarDensity(fluidState, paramCache, phaseIdx);
            fluidState.setMolarDensity(phaseIdx, value);
        }
    }
};

} // end namespace Dumux

#endif
