//@HEADER
// ************************************************************************
//
//                     		       Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Chris Wentland (crwentl@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef PRESSIODEMOAPPS_SCHWARZ_CUSTOMBCS_HPP_
#define PRESSIODEMOAPPS_SCHWARZ_CUSTOMBCS_HPP_

#include "pressiodemoapps/impl/ghost_relative_locations.hpp"
#include "pressiodemoapps/euler2d.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "pressiodemoapps/advection_diffusion2d.hpp"


namespace pschwarz{

namespace pda = pressiodemoapps;

enum class BCType {
    HomogNeumannVert,
    HomogNeumannHoriz,
    HomogDirichletVert,
    HomogDirichletHoriz,
    SlipWallVert,
    SlipWallHoriz,
    SchwarzDirichlet,
};

template<class mesh_t>
struct BCFunctor
{
    using graph_t  = typename mesh_t::graph_t;
    using scalar_t = typename mesh_t::scalar_t;
    // TODO: not sure if there's a way to template state_t, since app type is templated on BCFunctor (circular?)
    using state_t  = Eigen::Matrix<scalar_t,-1,1>;

    BCType m_bcSwitch;

    // m_stateBcs is shared by EVERY BCFunctor
    // m_graphBcs is a unique mask on m_stateBcs, for each BCFunctor
    // This is because BCFunctor has no internal left/right/front/back reference, and
    //  only understand position in graph from GLOBAL (subdomain) cell index
    state_t* m_stateBcs = nullptr;
    graph_t* m_graphBcs = nullptr;

    BCFunctor(BCType bcSwitch) : m_bcSwitch(bcSwitch){}

    void setInternalPtr(state_t* stateBcs){
        m_stateBcs = stateBcs;
    }

    void setInternalPtr(graph_t* graphBcs){
        m_graphBcs = graphBcs;
    }

    template<class ...Args>
    void operator()(Args && ... args) const
    {
        switch(m_bcSwitch)
        {
            case BCType::HomogNeumannVert:
                HomogNeumannVertBC(std::forward<Args>(args)...);
                break;
            case BCType::HomogNeumannHoriz:
                HomogNeumannHorizBC(std::forward<Args>(args)...);
                break;
            case BCType::HomogDirichletVert:
                HomogDirichletVertBC(std::forward<Args>(args)...);
                break;
            case BCType::HomogDirichletHoriz:
                HomogDirichletHorizBC(std::forward<Args>(args)...);
                break;
            case BCType::SlipWallVert:
                SlipWallVertBC(std::forward<Args>(args)...);
                break;
            case BCType::SlipWallHoriz:
                SlipWallHorizBC(std::forward<Args>(args)...);
                break;
            case BCType::SchwarzDirichlet:
                SchwarzDirichletBC(std::forward<Args>(args)...);
                break;
            default:
          throw std::runtime_error("Invalid probId for getPhysBCs()");
        };
    }

private:

    /*=========================
        PHYSICAL BOUNDARIES
    =========================*/

    template<class ConnecRowType, class StateT, class T>
    void HomogNeumannVertBC(
        const int /*unused*/, ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        const StateT & currentState, int numDofPerCell,
        const double cellWidth, T & ghostValues) const
    {

        // TODO: generalize to 1D/3D

        // this operates under the assumption that this cell does not have ghost cells in two parallel walls
        int stencilSize1D = ghostValues.cols() / numDofPerCell;
        const int cellGID = connectivityRow[0];
        const auto uIndex  = cellGID * numDofPerCell;

        const auto left0  = connectivityRow[1];
        const auto right0  = connectivityRow[3];
        if ((left0 == -1) && (right0 == -1)) {
            throw std::runtime_error("Should not have walls to left and right of same cell");
        }

        if ((left0 == -1) || (right0 == -1)) {
            for (int i = 0; i < numDofPerCell; ++i) {
                ghostValues[i] = currentState(uIndex+i);
            }
        }

        // TODO: extend to WENO5
        if (stencilSize1D > 1) {
            const auto left1  = connectivityRow[5];
            const auto right1  = connectivityRow[7];
            if ((left1 == -1) && (right1 == -1)) {
                throw std::runtime_error("Should not have walls to left and right of same cell");
            }

            // TODO: I don't think this is actually valid for cells that are more than 1 cell away from the boundary?
            if (left1 == -1) {
                const auto ind = right0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
            }
            if (right1 == -1) {
                const auto ind = left0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
            }
        }

    }

    template<class ConnecRowType, class FactorsType>
    void HomogNeumannVertBC(
        ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        int numDofPerCell, FactorsType & factorsForBCJac) const
    {
        for (int i = 0; i < numDofPerCell; ++i) {
            factorsForBCJac[i] = 1.0;
        }
    }

    template<class ConnecRowType, class StateT, class T>
    void HomogNeumannHorizBC(
        const int /*unused*/, ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        const StateT & currentState, int numDofPerCell,
        const double cellWidth, T & ghostValues) const
    {
        // TODO: generalize to 3D

        // this operates under the assumption that this cell does not have ghost cells in two parallel walls
        int stencilSize1D = ghostValues.cols() / numDofPerCell;
        const int cellGID = connectivityRow[0];
        const auto uIndex  = cellGID * numDofPerCell;

        const auto front0  = connectivityRow[2];
        const auto back0  = connectivityRow[4];
        if ((front0 == -1) && (back0 == -1)) {
            throw std::runtime_error("Should not have walls to left and right of same cell");
        }

        if ((front0 == -1) || (back0 == -1)) {
            for (int i = 0; i < numDofPerCell; ++i) {
                ghostValues[i] = currentState(uIndex+i);
            }
        }

        // TODO: extend to WENO5
        if (stencilSize1D > 1) {
            const auto front1  = connectivityRow[6];
            const auto back1  = connectivityRow[8];
            if ((front1 == -1) && (back1 == -1)) {
                throw std::runtime_error("Should not have walls to left and right of same cell");
            }

            if (front1 == -1) {
                const auto ind = back0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
            }
            if (back1 == -1) {
                const auto ind = front0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
            }
        }
    }

    template<class ConnecRowType, class FactorsType>
    void HomogNeumannHorizBC(
        ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        int numDofPerCell, FactorsType & factorsForBCJac) const
    {
        for (int i = 0; i < numDofPerCell; ++i) {
            factorsForBCJac[i] = 1.0;
        }
    }

    template<class ConnecRowType, class StateT, class T>
    void HomogDirichletVertBC(
        const int /*unused*/, ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        const StateT & currentState, int numDofPerCell,
        const double cellWidth, T & ghostValues) const
    {
        // this operates under the assumption that this cell does not have ghost cells in two parallel walls
        int stencilSize1D = ghostValues.cols() / numDofPerCell;

        const auto left0  = connectivityRow[1];
        const auto right0  = connectivityRow[3];
        if ((left0 == -1) && (right0 == -1)) {
            throw std::runtime_error("Should not have walls to left and right of same cell");
        }
        if ((left0 == -1) || (right0 == -1)) {
            for (int i = 0; i < numDofPerCell; ++i) {
                ghostValues[i] = 0.0;
            }
        }

        // TODO: extend to WENO5
        if (stencilSize1D > 1) {
            const auto left1  = connectivityRow[5];
            const auto right1  = connectivityRow[7];
            if ((left1 == -1) && (right1 == -1)) {
                throw std::runtime_error("Should not have walls to left and right of same cell");
            }

            // TODO: I don't think this is actually valid for cells that are more than 1 cell away from the boundary?
            if (left1 == -1) {
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = 0.0;
                }
            }
            if (right1 == -1) {
                const auto ind = left0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = 0.0;
                }
            }
        }
    }

    template<class ConnecRowType, class FactorsType>
    void HomogDirichletVertBC(
        ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        int numDofPerCell, FactorsType & factorsForBCJac) const
    {
        for (int i = 0; i < numDofPerCell; ++i) {
            factorsForBCJac[i] = 0.0;
        }
    }

    template<class ConnecRowType, class StateT, class T>
    void HomogDirichletHorizBC(
        const int /*unused*/, ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        const StateT & currentState, int numDofPerCell,
        const double cellWidth, T & ghostValues) const
    {
        // TODO: generalize to 3D

        // this operates under the assumption that this cell does not have ghost cells in two parallel walls
        int stencilSize1D = ghostValues.cols() / numDofPerCell;

        const auto front0  = connectivityRow[2];
        const auto back0  = connectivityRow[4];
        if ((front0 == -1) && (back0 == -1)) {
            throw std::runtime_error("Should not have walls to left and right of same cell");
        }

        if ((front0 == -1) || (back0 == -1)) {
            for (int i = 0; i < numDofPerCell; ++i) {
                ghostValues[i] = 0.0;
            }
        }

        // TODO: extend to WENO5
        if (stencilSize1D > 1) {
            const auto front1  = connectivityRow[6];
            const auto back1  = connectivityRow[8];
            if ((front1 == -1) && (back1 == -1)) {
                throw std::runtime_error("Should not have walls to left and right of same cell");
            }

            if (front1 == -1) {
                const auto ind = back0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = 0.0;
                }
            }
            if (back1 == -1) {
                const auto ind = front0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = 0.0;
                }
            }
        }
    }

    template<class ConnecRowType, class FactorsType>
    void HomogDirichletHorizBC(
        ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        int numDofPerCell, FactorsType & factorsForBCJac) const
    {
        for (int i = 0; i < numDofPerCell; ++i) {
            factorsForBCJac[i] = 0.0;
        }
    }

    template<class ConnecRowType, class StateT, class T>
    void SlipWallVertBC(
        const int /*unused*/, ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        const StateT & currentState, int numDofPerCell,
        const double cellWidth, T & ghostValues) const
    {
        // TODO: generalize to 1D/3D

        // this operates under the assumption that this cell does not have ghost cells in two parallel walls
        int stencilSize1D = ghostValues.cols() / numDofPerCell;
        const int cellGID = connectivityRow[0];
        const auto uIndex  = cellGID * numDofPerCell;

        const auto left0  = connectivityRow[1];
        const auto right0  = connectivityRow[3];
        if ((left0 == -1) && (right0 == -1)) {
            throw std::runtime_error("Should not have walls to left and right of same cell");
        }

        if ((left0 == -1) || (right0 == -1)) {
            for (int i = 0; i < numDofPerCell; ++i) {
                ghostValues[i] = currentState(uIndex+i);
            }
            ghostValues[1] *= -1.0; // reverse x-momentum
        }

        // TODO: extend to WENO5
        if (stencilSize1D > 1) {
            const auto left1  = connectivityRow[5];
            const auto right1  = connectivityRow[7];
            if ((left1 == -1) && (right1 == -1)) {
                throw std::runtime_error("Should not have walls to left and right of same cell");
            }

            // TODO: I don't think this is actually valid for cells that are more than 1 cell away from the boundary?
            if (left1 == -1) {
                const auto ind = right0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
                ghostValues[numDofPerCell + 1] *= -1.0; // reverse x-momentum
            }
            if (right1 == -1) {
                const auto ind = left0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
                ghostValues[numDofPerCell + 1] *= -1.0; // reverse x-momentum
            }
        }
    }

    template<class ConnecRowType, class FactorsType>
    void SlipWallVertBC(
        ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        int numDofPerCell, FactorsType & factorsForBCJac) const
    {
        for (int i = 0; i < numDofPerCell; ++i) {
            factorsForBCJac[i] = 1.0;
        }
        factorsForBCJac[1] = -1.0;
    }

    template<class ConnecRowType, class StateT, class T>
    void SlipWallHorizBC(
        const int /*unused*/, ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        const StateT & currentState, int numDofPerCell,
        const double cellWidth, T & ghostValues) const
    {
        // TODO: generalize to 3D

        // this operates under the assumption that this cell does not have ghost cells in two parallel walls
        int stencilSize1D = ghostValues.cols() / numDofPerCell;
        const int cellGID = connectivityRow[0];
        const auto uIndex  = cellGID * numDofPerCell;

        const auto front0  = connectivityRow[2];
        const auto back0  = connectivityRow[4];
        if ((front0 == -1) && (back0 == -1)) {
            throw std::runtime_error("Should not have walls to left and right of same cell");
        }

        if ((front0 == -1) || (back0 == -1)) {
            for (int i = 0; i < numDofPerCell; ++i) {
                ghostValues[i] = currentState(uIndex+i);
            }
            ghostValues[2] *= -1.0; // reverse y-momentum
        }

        // TODO: extend to WENO5
        if (stencilSize1D > 1) {
            const auto front1  = connectivityRow[6];
            const auto back1  = connectivityRow[8];
            if ((front1 == -1) && (back1 == -1)) {
                throw std::runtime_error("Should not have walls to left and right of same cell");
            }

            if (front1 == -1) {
                const auto ind = back0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
                ghostValues[numDofPerCell + 2] *= -1.0; // reverse y-momentum
            }
            if (back1 == -1) {
                const auto ind = front0*numDofPerCell;
                for (int i = 0; i < numDofPerCell; ++i) {
                    ghostValues[numDofPerCell + i] = currentState(ind+i);
                }
                ghostValues[numDofPerCell + 2] *= -1.0; // reverse y-momentum
            }
        }
    }

    template<class ConnecRowType, class FactorsType>
    void SlipWallHorizBC(
        ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        int numDofPerCell, FactorsType & factorsForBCJac) const
    {
        for (int i = 0; i < numDofPerCell; ++i) {
            factorsForBCJac[i] = 1.0;
        }
        factorsForBCJac[2] = -1.0;
    }

    /*=========================
        SCHWARZ BOUNDARIES
    =========================*/

    template<class ConnecRowType, class StateT, class T>
    void SchwarzDirichletBC(
        const int gRow, ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        const StateT & currentState, int numDofPerCell,
        const double cellWidth, T & ghostValues) const
    {
        if ((m_stateBcs == nullptr) || (m_graphBcs == nullptr)) {
            std::runtime_error("m_stateBcs or m_graphBcs not set");
        }

        // gRow: the index of current cell within graphRowsOfCellsNearBd()
        // connectivityRow: the stencil mesh graph associated with the current cell
        // ghostValues: the row of m_ghost(Left/Right/etc) associated with this cell

        int stencilSize1D = ghostValues.cols() / numDofPerCell;
        for (int stencilIdx = 0; stencilIdx < stencilSize1D; ++stencilIdx) {
            auto bcIdx = (*m_graphBcs)(gRow, stencilIdx);
            if (bcIdx != -1) {
                for (int dofIdx = 0; dofIdx < numDofPerCell; ++dofIdx) {
                    ghostValues[stencilIdx * numDofPerCell + dofIdx] = (*m_stateBcs)(bcIdx * numDofPerCell + dofIdx);
                }
            }
        }

    }

    template<class ConnecRowType, class FactorsType>
    void SchwarzDirichletBC(
        ConnecRowType const & connectivityRow,
        const double cellX, const double cellY,
        int numDofPerCell, FactorsType & factorsForBCJac) const
    {
        // assumes that FactorsType can be indexed by [], which is true for demoapps (std::array)
        for (int i = 0; i < numDofPerCell; ++i) {
            factorsForBCJac[i] = 0.0;
        }
    }

};

/*============================
    DEFAULT SPECIFICATIONS
=============================*/

auto getPhysBCs(pda::Euler2d probId, pda::impl::GhostRelativeLocation rloc)
{

    switch(probId)
    {
        case pda::Euler2d::Riemann:
            // All boundaries are homogeneous Neumann
            if ((rloc == pda::impl::GhostRelativeLocation::Left) || (rloc == pda::impl::GhostRelativeLocation::Right)) {
                return BCType::HomogNeumannVert;
            }
            else if ((rloc == pda::impl::GhostRelativeLocation::Front) || (rloc == pda::impl::GhostRelativeLocation::Back)) {
                return BCType::HomogNeumannHoriz;
            }
            else {
                throw std::runtime_error("Unexpected GhostRelativeLocation");
            }
            break;

        default:
            throw std::runtime_error("Invalid probId for getPhysBCs()");

    }
}

auto getPhysBCs(pda::Swe2d probId, pda::impl::GhostRelativeLocation rloc)
{

    switch(probId)
    {
        // SWE is really only SlipWall for now, but PDA doesn't permit Swe2d::SlipWall with custom BCs
        // Just use CustomBCs as a proxy for the normal slip wall case
        case pda::Swe2d::CustomBCs:
            if ((rloc == pda::impl::GhostRelativeLocation::Left) || (rloc == pda::impl::GhostRelativeLocation::Right)) {
                return BCType::SlipWallVert;
            }
            else if ((rloc == pda::impl::GhostRelativeLocation::Front) || (rloc == pda::impl::GhostRelativeLocation::Back)) {
                return BCType::SlipWallHoriz;
            }
            else {
                throw std::runtime_error("Unexpected GhostRelativeLocation");
            }
            break;

        default:
            throw std::runtime_error("Invalid probId for getPhysBCs()");

    }
}

auto getPhysBCs(pda::AdvectionDiffusion2d probId, pda::impl::GhostRelativeLocation rloc)
{

    switch(probId)
    {
        case pda::AdvectionDiffusion2d::BurgersOutflow:
            // left and bottom are homogeneous Dirichlet
            // top and right are homogeneous Neumann
            if (rloc == pda::impl::GhostRelativeLocation::Left) {
                return BCType::HomogDirichletVert;
            }
            else if (rloc == pda::impl::GhostRelativeLocation::Right) {
                return BCType::HomogNeumannVert;
            }
            else if (rloc == pda::impl::GhostRelativeLocation::Front) {
                return BCType::HomogNeumannHoriz;
            }
            else if (rloc == pda::impl::GhostRelativeLocation::Back) {
                return BCType::HomogDirichletHoriz;
            }
            else {
                throw std::runtime_error("Unexpected GhostRelativeLocation");
            }
            break;

        default:
            throw std::runtime_error("Invalid probId for getPhysBCs()");

    }
}

}

#endif
