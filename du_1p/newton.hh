// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © Timo Koch
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_LOGARITHM_TRANSFORMATION_NEWTON_SOLVER_HH
#define DUMUX_LOGARITHM_TRANSFORMATION_NEWTON_SOLVER_HH

#include <algorithm>

#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/assembly/partialreassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/discretization/elementsolution.hh>

namespace Dumux {

template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator> >
class LogarithmTransformationNewtonSolver : public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using Scalar = typename Assembler::Scalar;
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using Indices = typename Assembler::GridVariables::VolumeVariables::Indices;

    using typename ParentType::Backend;
    using typename ParentType::SolutionVector;
    using typename ParentType::ResidualVector;

public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

private:

    void newtonBeginStep(const Variables& vars) override
    {
        ParentType::newtonBeginStep(vars);
        uLastIter_ = Backend::dofs(vars);
    }

    /*!
     * \brief Solve the linear system of equations \f$\mathbf{A}x - b = 0\f$.
     * Throws Dumux::NumericalProblem if the linear solver didn't converge.
     */
    bool solveLinearSystem_(ResidualVector& deltaU) override
    {
        // auto jacobian = this->assembler().jacobian();  // Get the Jacobian matrix
        // // Print custom text
        // std::cout << "Printing the JC matrix:\n";
        // // Iterate over rows and columns to print the matrix
        // for (const auto& row : jacobian) {
        //     for (const auto& elem : row) {
        //         std::cout << elem << " ";  // Print each element in the row
        //     }
        //     std::cout << std::endl;  // Print a newline after each row
        // }


        static const bool useLogTranformation = getParam<bool>("Newton.UseLogTransformation", false);
        if (!useLogTranformation){
            printf("JC\n");
            Dune::printmatrix(std::cout, this->assembler().jacobian(), "", "", 10/*width*/, 2/*precision*/);
            printf("deltaU\n");
            Dune::printvector(std::cout, deltaU, "", "", 10/*width*/, 2/*precision*/);       
            printf("residual\n");
            Dune::printvector(std::cout, this->assembler().residual(), "", "", 10/*width*/, 2/*precision*/);

            return this->linearSolver().solve(
                this->assembler().jacobian(),deltaU, this->assembler().residual());
				
        }

        // Solve the linear system in a modified way
        auto JLogC = this->assembler().jacobian(); // copy
        auto deltaLogC = deltaU; // copy
        auto residual = this->assembler().residual(); // copy
        // modify the Jacobian by multiplying each molefraction derivative by the molefraction itself
        // JLogC = dR/dlnC = J * dC/dlnC
        // dx/dlnx = (dlnx/dx)^-1 = (1/x)^-1 = x
        // for (auto rowIt = JLogC.begin(); rowIt != JLogC.end(); ++rowIt)
        //     for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
        //         for (int i = 0; i < colIt->M(); ++i)
        //             for (int j = 1; j < colIt->N(); ++j)
        //                 (*colIt)[i][j] *= uLastIter_[colIt.index()][j];

        for (auto i = 0UL; i < residual.size(); ++i)
            for (int j = 1; j < residual[i].size(); ++j)
                residual[i][j] /= uLastIter_[i][j];

        bool converged = this->linearSolver().solve(
            JLogC,
            deltaLogC,
            residual
        );
        printf("JLogC\n");
        Dune::printmatrix(std::cout, JLogC, "", "", 10/*width*/, 2/*precision*/);
        printf("deltaLogC\n");
        Dune::printvector(std::cout, deltaLogC, "", "", 10/*width*/, 2/*precision*/);       
        printf("residual\n");
        Dune::printvector(std::cout, residual, "", "", 10/*width*/, 2/*precision*/);

// 
        // bool converged = this->linearSolver().solve(
        //     JLogC,
        //     deltaLogC,
        //     this->assembler().residual()
        // );
// 
        if (!converged)
            DUNE_THROW(Dumux::NumericalProblem, "Linear solver did not converge.");

        // Transform the solution back to the original space
        for (int i = 0; i < deltaU.size(); ++i)
        {
            deltaU[i] = deltaLogC[i]; // copy back

            // overwrite the solution with the transformed solution for the mole fractions
            // lnCNew = lnCOld - deltaLnC
            // => CNew = COld * exp(-deltaLnC) = COld - deltaC
            // => deltaC = COld - COld * exp(-deltaLnC)
            for (int j = 1; j < deltaU[i].size(); ++j){
                // deltaLogC[i][j] = std::copysign(std::min(5.0, std::abs(deltaLogC[i][j])), deltaLogC[i][j]);
                deltaU[i][j] = uLastIter_[i][j] -  uLastIter_[i][j]*std::exp(-deltaLogC[i][j]);
        }
        }
        printf("deltaU\n");
        Dune::printvector(std::cout, deltaU, "", "", 10/*width*/, 2/*precision*/); 
        return converged;
    }
private:
    SolutionVector uLastIter_;

};

} // end namespace Dumux

#endif
