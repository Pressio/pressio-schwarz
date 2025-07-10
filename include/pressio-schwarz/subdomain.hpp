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

#ifndef PRESSIODEMOAPPS_SCHWARZ_SUBDOMAIN_HPP_
#define PRESSIODEMOAPPS_SCHWARZ_SUBDOMAIN_HPP_

#include <string>
#include <vector>
#include <math.h>

#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"

#include "./tiling.hpp"
#include "./custom_bcs.hpp"
#include "./rom_utils.hpp"


namespace pschwarz {

namespace pda = pressiodemoapps;
namespace pode = pressio::ode;
namespace pls = pressio::linearsolvers;
namespace prom = pressio::rom;
namespace plspg = pressio::rom::lspg;

template<class mesh_t>
std::array<int, 3> calc_mesh_dims(const mesh_t & meshObj)
{
    // WARNING: this is only valid for full meshes
    const auto dx = meshObj.dx();
    const auto dy = meshObj.dy();
    const auto dz = meshObj.dz();
    const auto & xcoords = meshObj.viewX();
    const auto & ycoords = meshObj.viewY();
    const auto & zcoords = meshObj.viewZ();
    const auto ndim = meshObj.dimensionality();

    std::array<int, 3> meshdims = {0, 0, 0};

    // number cells in x
    auto xmin = xcoords.minCoeff();
    auto xmax = xcoords.maxCoeff();
    meshdims[0] = std::round((xmax - xmin) / dx) + 1;

    if (ndim > 1) {
        auto ymin = ycoords.minCoeff();
        auto ymax = ycoords.maxCoeff();
        meshdims[1] = std::round((ymax - ymin) / dy) + 1;
    }
    if (ndim == 3) {
        auto zmin = zcoords.minCoeff();
        auto zmax = zcoords.maxCoeff();
        meshdims[2] = std::round((zmax - zmin) / dz) + 1;
    }

    return meshdims;
}

template<class mesh_type, class state_type>
class SubdomainBase{
public:
    using state_t = state_type;
    using mesh_t = mesh_type;
    using graph_t = typename mesh_t::graph_t;
    using stencil_t  = decltype(create_cell_gids_vector_and_fill_from_ascii(std::declval<std::string>()));

    virtual void allocateStorageForHistory(const int) = 0;
    virtual void doStep(pode::StepStartAt<double>, pode::StepCount, pode::StepSize<double>) = 0;
    virtual void storeStateHistory(const int) = 0;
    virtual void resetStateFromHistory() = 0;
    virtual void updateFullState() = 0;
    virtual const mesh_t & getMeshStencil() const = 0;
    virtual const mesh_t & getMeshFull() const = 0;
    virtual const std::array<int, 3> getFullMeshDims() const = 0;
    virtual const stencil_t * getSampleGids() const = 0;
    virtual void setStencilGids(std::vector<int>) = 0;
    virtual void genHyperMesh(std::string &) = 0;
    virtual void setNeighborGraph(graph_t &) = 0;
    virtual const graph_t & getNeighborGraph() const = 0;
    virtual int getDofPerCell() const = 0;
    virtual void finalize_subdomain(std::string &) = 0;
    virtual state_t * getStateStencil() = 0;
    virtual state_t * getStateFull() = 0;
    virtual state_t * getStateReduced() = 0;
    virtual state_t * getStateBCs() = 0;
    virtual void setBCPointer(pda::impl::GhostRelativeLocation, state_t * ) = 0;
    virtual void setBCPointer(pda::impl::GhostRelativeLocation, graph_t *) = 0;
    virtual state_t & getLastStateInHistory() = 0;
};


template<class mesh_t, class app_type, class prob_t>
class SubdomainFOM: public SubdomainBase<mesh_t, typename app_type::state_type>
{
    using base_t = SubdomainBase<mesh_t, typename app_type::state_type>;

public:
    using app_t   = app_type;
    using scalar_t = typename app_t::scalar_type;
    using graph_t = typename mesh_t::graph_t;
    using state_t = typename app_t::state_type;
    using jacob_t = typename app_t::jacobian_type;
    using stencil_t = typename base_t::stencil_t;

    using stepper_t  =
        decltype(pressio::ode::create_implicit_stepper(pressio::ode::StepScheme(),
            std::declval<app_t&>())
        );

    using lin_solver_tag = pressio::linearsolvers::iterative::Bicgstab;
    using linsolver_t    = pressio::linearsolvers::Solver<lin_solver_tag, jacob_t>;
    using nonlinsolver_t =
        decltype( pressio::create_newton_solver( std::declval<stepper_t &>(),
                            std::declval<linsolver_t&>()) );

public:
    SubdomainFOM(
        const int domainIndex,
        const mesh_t & mesh,
        BCType bcLeft, BCType bcFront,
        BCType bcRight, BCType bcBack,
        prob_t probId,
        pressio::ode::StepScheme odeScheme,
        pda::InviscidFluxReconstruction fluxOrder,
        const int icflag,
        const std::string & icFileRoot,
        const std::unordered_map<std::string, scalar_t> & userParams)
    : m_domIdx(domainIndex)
    , m_mesh(&mesh)
    , m_app(std::make_shared<app_t>(pda::create_problem_eigen(
            mesh, probId, fluxOrder,
            BCFunctor<mesh_t>(bcLeft),  BCFunctor<mesh_t>(bcFront),
            BCFunctor<mesh_t>(bcRight), BCFunctor<mesh_t>(bcBack),
            icflag, userParams)))
    , m_state(m_app->initialCondition())
    , m_stepper(pressio::ode::create_implicit_stepper(odeScheme, *(m_app)))
    , m_linSolverObj(std::make_shared<linsolver_t>())
    , m_nonlinSolver(pressio::create_newton_solver(m_stepper, *m_linSolverObj))
    {
        m_fullMeshDims = calc_mesh_dims(*m_mesh);

        pda::resize(m_sampleGids, m_mesh->sampleMeshSize());
        for (int i = 0; i < m_mesh->sampleMeshSize(); ++i) {
            m_sampleGids(i) = i;
        }

        if (!(icFileRoot.empty())) {
            // load from file
            std::string icFile = icFileRoot + std::to_string(domainIndex) + ".bin";
            auto instate = read_vector_from_binary<scalar_t>(icFile);
            int nrows = instate.rows();
            if (nrows == m_state.rows()) {
                m_state = instate;
            }
            else {
                throw std::runtime_error("Invalid icFile dimensions: " + std::to_string(nrows));
            }
        }

        m_nonlinSolver.setStopCriterion(pressio::nonlinearsolvers::Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance);
        m_nonlinSolver.setStopTolerance(1e-5);
    }

    state_t & getLastStateInHistory() final { return m_stateHistVec.back(); }

    void setBCPointer(pda::impl::GhostRelativeLocation grl, state_t * v) final {
        m_app->setBCPointer(grl, v);
    }
    void setBCPointer(pda::impl::GhostRelativeLocation grl, graph_t * v) final {
      m_app->setBCPointer(grl, v);
    }

    state_t * getStateBCs() final { return &m_stateBCs; }
    state_t * getStateStencil() final { return &m_state; }
    state_t * getStateFull() final { return &m_state; }
    state_t * getStateReduced() final {
        throw std::runtime_error("Called getStateReduced() on a FOM subdomain, fix");
    }

    int getDofPerCell() const final { return m_app->numDofPerCell(); }
    const mesh_t & getMeshStencil() const final { return *m_mesh; }
    const mesh_t & getMeshFull() const final { return *m_mesh; }
    const std::array<int, 3> getFullMeshDims() const final { return m_fullMeshDims; }
    const stencil_t * getSampleGids() const final { return &m_sampleGids; }
    void setStencilGids(std::vector<int> gids_vec) final { /*noop*/ }
    void genHyperMesh(std::string & subdom_dir) final { /*noop*/ }
    const graph_t & getNeighborGraph() const final { return m_neighborGraph; }

    void setNeighborGraph(graph_t & graph_in) {
        pda::resize(m_neighborGraph, graph_in.rows(), graph_in.cols());
        for (int rowIdx = 0; rowIdx < graph_in.rows(); ++rowIdx) {
            for (int colIdx = 0; colIdx < graph_in.cols(); ++colIdx) {
                m_neighborGraph(rowIdx, colIdx) = graph_in(rowIdx, colIdx);
            }
        }
    }

    void init_bc_state()
    {
        // count number of neighbor ghost cells in neighborGraph
        int numGhostCells = 0;
        const auto & rowsBd = m_mesh->graphRowsOfCellsNearBd();
        for (int bdIdx = 0; bdIdx < rowsBd.size(); ++bdIdx) {
            auto rowIdx = rowsBd[bdIdx];
            // start at 1 to ignore own ID
            for (int colIdx = 1; colIdx < m_neighborGraph.cols(); ++colIdx) {
                if (m_neighborGraph(rowIdx, colIdx) != -1) {
                    numGhostCells++;
                }
            }
        }
        const int numDofStencilBc = m_app->numDofPerCell() * numGhostCells;
        pda::resize(m_stateBCs, numDofStencilBc);
        m_stateBCs.fill(0.0);
    }

    void finalize_subdomain(std::string &) final {
        init_bc_state();
    }

    void allocateStorageForHistory(const int count) final {
        for (int histIdx = 0; histIdx < count + 1; ++histIdx) {
            // createState creates a new state with all elements equal to zero
            m_stateHistVec.emplace_back(m_app->createState());
        }
    }

    void doStep(pode::StepStartAt<double> startTime,
        pode::StepCount step,
        pode::StepSize<double> dt) final
    {
        m_stepper(m_state, startTime, step, dt, m_nonlinSolver);
    }

    void storeStateHistory(const int step) final {
        m_stateHistVec[step] = m_state;
    }

    void resetStateFromHistory() final {
        m_state = m_stateHistVec[0];
    }

    void updateFullState() final {
        // noop
    }

public:
    int m_domIdx;
    mesh_t const * m_mesh;
    std::array<int, 3> m_fullMeshDims;
    stencil_t m_sampleGids;
    graph_t m_neighborGraph;
    std::shared_ptr<app_t> m_app;
    state_t m_state;
    state_t m_stateBCs;
    std::vector<state_t> m_stateHistVec;

    stepper_t m_stepper;
    std::shared_ptr<linsolver_t> m_linSolverObj;
    nonlinsolver_t m_nonlinSolver;

};



template<class mesh_t, class app_type, class prob_t>
class SubdomainROM: public SubdomainBase<mesh_t, typename app_type::state_type>
{
    using base_t = SubdomainBase<mesh_t, typename app_type::state_type>;

public:
    using app_t    = app_type;
    using graph_t = typename mesh_t::graph_t;
    using scalar_t = typename app_t::scalar_type;
    using state_t  = typename app_t::state_type;
    using stencil_t = typename base_t::stencil_t;

    using trans_t = decltype(read_vector_from_binary<scalar_t>(std::declval<std::string>()));
    using basis_t = decltype(read_matrix_from_binary<scalar_t>(std::declval<std::string>(), std::declval<int>()));
    using trial_t = decltype(prom::create_trial_column_subspace<
        state_t>(std::declval<basis_t&&>(), std::declval<trans_t&&>(), true));

public:
    SubdomainROM(
        const int domainIndex,
        const mesh_t & mesh,
        BCType bcLeft, BCType bcFront,
        BCType bcRight, BCType bcBack,
        prob_t probId,
        pressio::ode::StepScheme odeScheme,
        pda::InviscidFluxReconstruction fluxOrder,
        const int icflag,
        const std::string icFileRoot,
        const std::unordered_map<std::string, typename mesh_t::scalar_type> & userParams,
        const std::string & transRoot,
        const std::string & basisRoot,
        int nmodes)
    : m_domIdx(domainIndex)
    , m_mesh(&mesh)
    , m_app(std::make_shared<app_t>(pda::create_problem_eigen(
            mesh, probId, fluxOrder,
            BCFunctor<mesh_t>(bcLeft),  BCFunctor<mesh_t>(bcFront),
            BCFunctor<mesh_t>(bcRight), BCFunctor<mesh_t>(bcBack),
            icflag, userParams)))
    , m_state(m_app->initialCondition())
    , m_nmodes(nmodes)
    , m_trans(read_vector_from_binary<scalar_t>(transRoot + "_" + std::to_string(domainIndex) + ".bin"))
    , m_basis(read_matrix_from_binary<scalar_t>(basisRoot + "_" + std::to_string(domainIndex) + ".bin", nmodes))
    , m_trialSpace(prom::create_trial_column_subspace<state_t>(std::move(m_basis), std::move(m_trans), true))
    , m_stateReduced(m_trialSpace.createReducedState())
    {
        m_fullMeshDims = calc_mesh_dims(*m_mesh);

        pda::resize(m_sampleGids, m_mesh->sampleMeshSize());
        for (int i = 0; i < m_mesh->sampleMeshSize(); ++i) {
            m_sampleGids(i) = i;
        }

        // initial conditions
        if (icFileRoot.empty()) {
            // project full state initial conditions
            auto u = pressio::ops::clone(m_state);
            pressio::ops::update(u, 0., m_state, 1, m_trialSpace.translationVector(), -1);
            pressio::ops::product(::pressio::transpose(), 1., m_trialSpace.basis(), u, 0., m_stateReduced);
        }
        else {
            // load from file
            std::string icFile = icFileRoot + std::to_string(domainIndex) + ".bin";
            auto instate = read_vector_from_binary<scalar_t>(icFile);
            int nrows = instate.rows();
            if (nrows == m_stateReduced.rows()) {
                m_stateReduced = instate;
            }
            else if (nrows == m_state.rows()) {
                // project full state initial conditions
                auto u = pressio::ops::clone(instate);
                pressio::ops::update(u, 0., instate, 1, m_trialSpace.translationVector(), -1);
                pressio::ops::product(::pressio::transpose(), 1., m_trialSpace.basis(), u, 0., m_stateReduced);
            }
            else {
                throw std::runtime_error("Invalid icFile dimensions: " + std::to_string(nrows));
            }
        }
        m_trialSpace.mapFromReducedState(m_stateReduced, m_state);

    }

    state_t & getLastStateInHistory() final { return m_stateHistVec.back(); }
    void setBCPointer(pda::impl::GhostRelativeLocation grl, state_t * v) final{
        m_app->setBCPointer(grl, v);
    }
    void setBCPointer(pda::impl::GhostRelativeLocation grl, graph_t * v) final{
        m_app->setBCPointer(grl, v);
    }

    state_t * getStateBCs() final { return &m_stateBCs; }
    state_t * getStateStencil() final { return &m_state; }
    state_t * getStateFull() final { return &m_state; }
    state_t * getStateReduced() final { return &m_stateReduced; }

    int getDofPerCell() const final { return m_app->numDofPerCell(); }
    const mesh_t & getMeshStencil() const final { return *m_mesh; }
    const mesh_t & getMeshFull() const final { return *m_mesh; }
    const std::array<int, 3> getFullMeshDims() const final { return m_fullMeshDims; }
    const stencil_t * getSampleGids() const final { return &m_sampleGids; }
    void setStencilGids(std::vector<int> gids_vec) final { /*noop*/ }
    void genHyperMesh(std::string & subdom_dir) final { /*noop*/ }
    const graph_t & getNeighborGraph() const final { return m_neighborGraph; }

    void setNeighborGraph(graph_t & graph_in) {
        pda::resize(m_neighborGraph, graph_in.rows(), graph_in.cols());
        for (int rowIdx = 0; rowIdx < graph_in.rows(); ++rowIdx) {
            for (int colIdx = 0; colIdx < graph_in.cols(); ++colIdx) {
                m_neighborGraph(rowIdx, colIdx) = graph_in(rowIdx, colIdx);
            }
        }
    }

    void init_bc_state()
    {
        // count number of neighbor ghost cells in neighborGraph
        int numGhostCells = 0;
        const auto & rowsBd = m_mesh->graphRowsOfCellsNearBd();
        for (int bdIdx = 0; bdIdx < rowsBd.size(); ++bdIdx) {
            auto rowIdx = rowsBd[bdIdx];
            // start at 1 to ignore own ID
            for (int colIdx = 1; colIdx < m_neighborGraph.cols(); ++colIdx) {
                if (m_neighborGraph(rowIdx, colIdx) != -1) {
                    numGhostCells++;
                }
            }
        }
        const int numDofStencilBc = m_app->numDofPerCell() * numGhostCells;
        pda::resize(m_stateBCs, numDofStencilBc);
        m_stateBCs.fill(0.0);
    }

    void finalize_subdomain(std::string &) final {
        init_bc_state();
    }

    void allocateStorageForHistory(const int count){
        for (int histIdx = 0; histIdx < count + 1; ++histIdx) {
            m_stateHistVec.emplace_back(m_app->createState());
            m_stateReducedHistVec.emplace_back(m_trialSpace.createReducedState());
        }
    }

    void storeStateHistory(const int step) final {
        m_stateHistVec[step] = m_state;
        m_stateReducedHistVec[step] = m_stateReduced;
    }

    void resetStateFromHistory() final {
        m_state = m_stateHistVec[0];
        m_stateReduced = m_stateReducedHistVec[0];
    }

    void updateFullState() final {
        m_trialSpace.mapFromReducedState(m_stateReduced, m_state);
    }

protected:
    int m_domIdx;
    mesh_t const * m_mesh;
    std::array<int, 3> m_fullMeshDims;
    stencil_t m_sampleGids;
    graph_t m_neighborGraph;
    std::shared_ptr<app_t> m_app;
    state_t m_state;
    state_t m_stateBCs;
    std::vector<state_t> m_stateHistVec;

    int m_nmodes;
    trans_t m_trans;
    basis_t m_basis;
    trial_t m_trialSpace;
    state_t m_stateReduced;
    std::vector<state_t> m_stateReducedHistVec;

};

template<class mesh_t, class app_type, class prob_t>
class SubdomainLSPG: public SubdomainROM<mesh_t, app_type, prob_t>
{
    using base_t = SubdomainROM<mesh_t, app_type, prob_t>;

public:
    using app_t    = app_type;
    using graph_t  = typename mesh_t::graph_t;
    using scalar_t = typename app_t::scalar_type;
    using state_t  = typename app_t::state_type;

    using trans_t = typename base_t::trans_t;
    using basis_t = typename base_t::basis_t;
    using trial_t = typename base_t::trial_t;

    using hessian_t   = Eigen::Matrix<scalar_t, -1, -1>; // TODO: generalize?
    using solver_tag  = pressio::linearsolvers::direct::HouseholderQR;
    using linsolver_t = pressio::linearsolvers::Solver<solver_tag, hessian_t>;

    using problem_t       = decltype(plspg::create_unsteady_problem(pressio::ode::StepScheme(), std::declval<trial_t&>(), std::declval<app_t&>()));
    using nonlinsolver_t  = decltype(pressio::create_gauss_newton_solver(std::declval<problem_t&>(), std::declval<linsolver_t&>()));

public:

    SubdomainLSPG(
        const int domainIndex,
        const mesh_t & mesh,
        BCType bcLeft, BCType bcFront,
        BCType bcRight, BCType bcBack,
        prob_t probId,
        pressio::ode::StepScheme odeScheme,
        pda::InviscidFluxReconstruction fluxOrder,
        const int icflag,
        const std::string & icFileRoot,
        const std::unordered_map<std::string, typename mesh_t::scalar_type> & userParams,
        const std::string & transRoot,
        const std::string & basisRoot,
        const int nmodes)
    : base_t(domainIndex, mesh,
             bcLeft, bcFront, bcRight, bcBack,
             probId, odeScheme, fluxOrder, icflag, icFileRoot, userParams,
             transRoot, basisRoot, nmodes)
    , m_problem(plspg::create_unsteady_problem(odeScheme, this->m_trialSpace, *(this->m_app)))
    , m_linSolverObj(std::make_shared<linsolver_t>())
    , m_nonlinSolver(pressio::create_gauss_newton_solver(m_problem, *m_linSolverObj))
    {

        m_nonlinSolver.setStopCriterion(pressio::nonlinearsolvers::Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance);
        m_nonlinSolver.setStopTolerance(1e-5);

    }

    void doStep(pode::StepStartAt<double> startTime, pode::StepCount step, pode::StepSize<double> dt) final {
        m_problem(this->m_stateReduced, startTime, step, dt, m_nonlinSolver);
    }

private:
    problem_t m_problem;
    std::shared_ptr<linsolver_t> m_linSolverObj;
    nonlinsolver_t m_nonlinSolver;
};


template<class mesh_t, class app_type, class prob_t>
class SubdomainHyper: public SubdomainBase<mesh_t, typename app_type::state_type>
{
    // TODO: really need to add some error checking
    //      to protect against the circuitious initialization

    using base_t = SubdomainBase<mesh_t, typename app_type::state_type>;

public:
    using app_t    = app_type;
    using graph_t  = typename mesh_t::graph_t;
    using scalar_t = typename app_t::scalar_type;
    using state_t  = typename app_t::state_type;

    using trans_t = decltype(read_vector_from_binary<scalar_t>(std::declval<std::string>()));
    using basis_t = decltype(read_matrix_from_binary<scalar_t>(std::declval<std::string>(), std::declval<int>()));
    using trial_t = decltype(prom::create_trial_column_subspace<
        state_t>(std::declval<basis_t&&>(), std::declval<trans_t&&>(), true));

    using stencil_t  = typename base_t::stencil_t;
    using transHyp_t = decltype(reduce_vector_on_stencil_mesh(std::declval<trans_t&>(), std::declval<stencil_t&>(), std::declval<int>()));
    using basisHyp_t = decltype(reduce_matrix_on_stencil_mesh(std::declval<basis_t&>(), std::declval<stencil_t&>(), std::declval<int>()));
    using trialHyp_t = decltype(prom::create_trial_column_subspace<
        state_t>(std::declval<basisHyp_t&&>(), std::declval<transHyp_t&&>(), true));

public:

    SubdomainHyper(
        const int domainIndex,
        const mesh_t & meshFull,
        BCType bcLeft, BCType bcFront,
        BCType bcRight, BCType bcBack,
        prob_t probId,
        pressio::ode::StepScheme odeScheme,
        pda::InviscidFluxReconstruction fluxOrder,
        const int icflag,
        const std::string & icFileRoot,
        const std::unordered_map<std::string, typename mesh_t::scalar_type> & userParams,
        const std::string & transRoot,
        const std::string & basisRoot,
        const int nmodes,
        const std::string & sampleFile)
    : m_domIdx(domainIndex)
    , m_meshFull(&meshFull)
    , m_probId(probId)
    , m_fluxOrder(fluxOrder)
    , m_bcLeft(bcLeft), m_bcFront(bcFront)
    , m_bcRight(bcRight), m_bcBack(bcBack)
    , m_icflag(icflag)
    , m_userParams(userParams)
    , m_appFull(std::make_shared<app_t>(pda::create_problem_eigen(
            meshFull, probId, fluxOrder,
            BCFunctor<mesh_t>(bcLeft),  BCFunctor<mesh_t>(bcFront),
            BCFunctor<mesh_t>(bcRight), BCFunctor<mesh_t>(bcBack),
            icflag, userParams)))
    , m_sampleFile(sampleFile)
    , m_sampleGids(create_cell_gids_vector_and_fill_from_ascii(m_sampleFile))
    , m_nmodes(nmodes)
    , m_transFull(read_vector_from_binary<scalar_t>(transRoot + "_" + std::to_string(domainIndex) + ".bin"))
    , m_basisFull(read_matrix_from_binary<scalar_t>(basisRoot + "_" + std::to_string(domainIndex) + ".bin", nmodes))
    , m_transRead(read_vector_from_binary<scalar_t>(transRoot + "_" + std::to_string(domainIndex) + ".bin"))
    , m_basisRead(read_matrix_from_binary<scalar_t>(basisRoot + "_" + std::to_string(domainIndex) + ".bin", nmodes))
    , m_trialSpaceFull(prom::create_trial_column_subspace<
       state_t>(std::move(m_basisFull), std::move(m_transFull), true))
    {

        m_stateFull = m_appFull->initialCondition();
        m_stateReduced = m_trialSpaceFull.createReducedState();

        m_fullMeshDims = calc_mesh_dims(*m_meshFull);

        // initial conditions
        if (icFileRoot.empty()) {
            // project full state initial conditions
            auto u = pressio::ops::clone(m_stateFull);
            pressio::ops::update(u, 0., m_stateFull, 1, m_trialSpaceFull.translationVector(), -1);
            pressio::ops::product(::pressio::transpose(), 1., m_trialSpaceFull.basis(), u, 0., m_stateReduced);
        }
        else {
            // load from file
            std::string icFile = icFileRoot + std::to_string(domainIndex) + ".bin";
            auto instate = read_vector_from_binary<scalar_t>(icFile);
            int nrows = instate.rows();
            if (nrows == m_stateReduced.rows()) {
                m_stateReduced = instate;
            }
            else if (nrows == m_stateFull.rows()) {
                // project full state initial conditions
                auto u = pressio::ops::clone(instate);
                pressio::ops::update(u, 0., instate, 1, m_trialSpaceFull.translationVector(), -1);
                pressio::ops::product(::pressio::transpose(), 1., m_trialSpaceFull.basis(), u, 0., m_stateReduced);
            }
            else {
                throw std::runtime_error("Invalid icFile dimensions: " + std::to_string(nrows));
            }
        }
        m_trialSpaceFull.mapFromReducedState(m_stateReduced, m_stateFull);

    }

    state_t & getLastStateInHistory() final { return m_stateHistVec.back(); }
    void setBCPointer(pda::impl::GhostRelativeLocation grl, state_t * v) final{
        m_appHyper->setBCPointer(grl, v);
    }
    void setBCPointer(pda::impl::GhostRelativeLocation grl, graph_t * v) final{
        m_appHyper->setBCPointer(grl, v);
    }

    state_t * getStateBCs() final { return &m_stateBCs; }
    state_t * getStateStencil() final { return &m_stateStencil; }
    state_t * getStateFull() final {
        m_trialSpaceFull.mapFromReducedState(m_stateReduced, m_stateFull);
        return &m_stateFull;
    }
    state_t * getStateReduced() final { return &m_stateReduced; }

    int getDofPerCell() const final { return m_appHyper->numDofPerCell(); }
    const mesh_t & getMeshStencil() const final {
        if (!m_hyperMeshSet) {
            throw std::runtime_error("Must call genHyperMesh() before getMeshStencil()");
        }
        return m_meshHyper;
    }
    const mesh_t & getMeshFull() const final { return *m_meshFull; }
    const std::array<int, 3> getFullMeshDims() const final { return m_fullMeshDims; }
    const stencil_t * getSampleGids() const final { return &m_sampleGids; }
    const graph_t & getNeighborGraph() const final { return m_neighborGraph; }

    void setStencilGids(std::vector<int> gids_vec) final {
        m_stencilGidsSet = true;

        pda::resize(m_stencilGids, (int) gids_vec.size());
        for (int stencilIdx = 0; stencilIdx < gids_vec.size(); ++stencilIdx) {
            m_stencilGids(stencilIdx) = gids_vec[stencilIdx];
        }
    }

    void genHyperMesh(std::string & subdom_dir) final {
        m_hyperMeshSet = true;

        m_meshHyper = pda::load_cellcentered_uniform_mesh_eigen(subdom_dir);
    }

    void setNeighborGraph(graph_t & graph_in) {
        pda::resize(m_neighborGraph, graph_in.rows(), graph_in.cols());
        for (int rowIdx = 0; rowIdx < graph_in.rows(); ++rowIdx) {
            for (int colIdx = 0; colIdx < graph_in.cols(); ++colIdx) {
                m_neighborGraph(rowIdx, colIdx) = graph_in(rowIdx, colIdx);
            }
        }
    }

    void init_bc_state()
    {

        if (!m_hyperMeshSet) {
            throw std::runtime_error("Must call genHyperMesh() before init_bc_state()");
        }

        // count number of neighbor ghost cells in neighborGraph
        int numGhostCells = 0;
        const auto & rowsBd = m_meshHyper.graphRowsOfCellsNearBd();
        for (int bdIdx = 0; bdIdx < rowsBd.size(); ++bdIdx) {
            auto rowIdx = rowsBd[bdIdx];
            // start at 1 to ignore own ID
            for (int colIdx = 1; colIdx < m_neighborGraph.cols(); ++colIdx) {
                if (m_neighborGraph(rowIdx, colIdx) != -1) {
                    numGhostCells++;
                }
            }
        }
        const int numDofStencilBc = m_appHyper->numDofPerCell() * numGhostCells;
        pda::resize(m_stateBCs, numDofStencilBc);
        m_stateBCs.fill(0.0);
    }

    // All this junk has to happen AFTER construction,
    //      as the hyper-reduction mesh hasn't been generated at that point
    void finalize_subdomain(std::string &) {

        if (!m_hyperMeshSet) {
            throw std::runtime_error("Must call genHyperMesh() before finalize_subdomain()");
        }
        if (!m_stencilGidsSet) {
            throw std::runtime_error("Must call setStencilGids() before finalize_subdomain()");
        }

        m_appHyper = std::make_shared<app_t>(pda::create_problem_eigen(
            m_meshHyper, m_probId, m_fluxOrder,
            BCFunctor<mesh_t>(m_bcLeft),  BCFunctor<mesh_t>(m_bcFront),
            BCFunctor<mesh_t>(m_bcRight), BCFunctor<mesh_t>(m_bcBack),
            m_icflag, m_userParams));

        auto m_transHyper = reduce_vector_on_stencil_mesh(m_transRead, m_stencilGids, m_appFull->numDofPerCell());
        auto m_basisHyper = reduce_matrix_on_stencil_mesh(m_basisRead, m_stencilGids, m_appFull->numDofPerCell());
        m_trialSpaceHyper = std::make_shared<trialHyp_t>(prom::create_trial_column_subspace<
                                                         state_t>(std::move(m_basisHyper),
                                                         std::move(m_transHyper),
                                                         true));

        // initialize stencil mesh state
        m_stateStencil = reduce_vector_on_stencil_mesh(m_stateFull, m_stencilGids, m_appFull->numDofPerCell());

        updateFullState();
        init_bc_state();

    }

    void allocateStorageForHistory(const int count){
        for (int histIdx = 0; histIdx < count + 1; ++histIdx) {
            m_stateHistVec.emplace_back(m_appHyper->createState());
            m_stateReducedHistVec.emplace_back(m_trialSpaceHyper->createReducedState());
        }
    }

    void storeStateHistory(const int step) final {
        m_stateHistVec[step] = m_stateStencil;
        m_stateReducedHistVec[step] = m_stateReduced;
    }

    void resetStateFromHistory() final {
        m_stateStencil = m_stateHistVec[0];
        m_stateReduced = m_stateReducedHistVec[0];
    }

    void updateFullState() final {
        m_trialSpaceHyper->mapFromReducedState(m_stateReduced, m_stateStencil);
    }

public:
    int m_domIdx;
    mesh_t const * m_meshFull;
    mesh_t m_meshHyper;
    std::array<int, 3> m_fullMeshDims;
    graph_t m_neighborGraph;
    std::shared_ptr<app_t> m_appFull;
    std::shared_ptr<app_t> m_appHyper;

    prob_t m_probId;
    pda::InviscidFluxReconstruction m_fluxOrder;
    int m_icflag;
    const std::unordered_map<std::string, typename mesh_t::scalar_type> m_userParams;
    BCType m_bcLeft;
    BCType m_bcFront;
    BCType m_bcRight;
    BCType m_bcBack;

    state_t m_stateStencil;  // on stencil mesh
    state_t m_stateFull;     // on full, unsampled mesh (required for projection)
    state_t m_stateReduced;  // latent state
    state_t m_stateBCs;
    std::vector<state_t> m_stateHistVec;
    std::vector<state_t> m_stateReducedHistVec;

    // for error checking
    bool m_hyperMeshSet = false;
    bool m_stencilGidsSet = false;

    std::string m_sampleFile;
    stencil_t m_sampleGids;
    stencil_t m_stencilGids;

    int m_nmodes;
    trans_t m_transFull;
    basis_t m_basisFull;
    trans_t m_transRead;
    basis_t m_basisRead;
    trial_t m_trialSpaceFull;

    std::shared_ptr<trialHyp_t> m_trialSpaceHyper;

};

template<class mesh_t, class app_type, class prob_t>
class SubdomainLSPGHyper: public SubdomainHyper<mesh_t, app_type, prob_t>
{
    using base_t = SubdomainHyper<mesh_t, app_type, prob_t>;

public:
    using app_t    = typename base_t::app_t;
    using graph_t  = typename mesh_t::graph_t;
    using scalar_t = typename app_t::scalar_type;
    using state_t  = typename app_t::state_type;
    using weigh_t  = Weigher<scalar_t>;

    using hessian_t   = Eigen::Matrix<scalar_t, -1, -1>; // TODO: generalize?
    using solver_tag  = pressio::linearsolvers::direct::HouseholderQR;
    using linsolver_t = pressio::linearsolvers::Solver<solver_tag, hessian_t>;

    using trialHyp_t = typename base_t::trialHyp_t;

    using updaterHyp_t = HypRedUpdater<scalar_t>;
    using problemHyp_t =
      decltype(plspg::create_unsteady_problem(pressio::ode::StepScheme(),
                                              std::declval<trialHyp_t&>(),
                                              std::declval<app_t&>(),
                                              std::declval<updaterHyp_t&>()));

    using tag_t = pressio::nonlinearsolvers::impl::CompactWeightedGaussNewtonNormalEqTag;

    using nonlinsolverHyp_t =
        decltype(pressio::create_gauss_newton_solver(
            std::declval<problemHyp_t&>(),
            std::declval<linsolver_t&>(),
            std::declval<weigh_t&>(),
            std::declval<tag_t>()
        ));

public:
    SubdomainLSPGHyper(
        const int domainIndex,
        const mesh_t & meshFull,
        BCType bcLeft, BCType bcFront,
        BCType bcRight, BCType bcBack,
        prob_t probId,
        pressio::ode::StepScheme odeScheme,
        pda::InviscidFluxReconstruction fluxOrder,
        const int icflag,
        const std::string & icFileRoot,
        const std::unordered_map<std::string, typename mesh_t::scalar_type> & userParams,
        const std::string & transRoot,
        const std::string & basisRoot,
        const int nmodes,
        const std::string & sampleFile,
        const std::string & weigher_type,
        const std::string & basisRoot_gpod,
        const int nmodes_gpod)
    : base_t(domainIndex, meshFull,
             bcLeft, bcFront, bcRight, bcBack,
             probId, odeScheme, fluxOrder, icflag, icFileRoot, userParams,
             transRoot, basisRoot, nmodes,
             sampleFile)
    {
        m_odeScheme = odeScheme;
        m_weigher_type = weigher_type;
        m_basisRoot_gpod = basisRoot_gpod;
        m_nmodes_gpod = nmodes_gpod;
    }

    void doStep(pode::StepStartAt<double> startTime, pode::StepCount step, pode::StepSize<double> dt) final {
        (*m_problemHyper)(this->m_stateReduced, startTime, step, dt, *m_nonlinSolverHyper);
    }

    // Again, this has to be done because the hyper-reduced mesh
    //      has not been initialized on construction
    // tempdir is the temporary directory that stores the true stencil mesh
    void finalize_subdomain(std::string & tempdir) final
    {
        SubdomainHyper<mesh_t, app_t, prob_t>::finalize_subdomain(tempdir);

        std::string stencilFile = tempdir + "/domain_" + std::to_string(this->m_domIdx) + "/stencil_mesh_gids.dat";

        m_updaterHyper = std::make_shared<updaterHyp_t>
            (create_hyper_updater<mesh_t>(this->getDofPerCell(),
                                          stencilFile,
                                          this->m_sampleFile));

        m_problemHyper = std::make_shared<problemHyp_t>
            (plspg::create_unsteady_problem(m_odeScheme,
                                            *(this->m_trialSpaceHyper),
                                            *(this->m_appHyper),
                                            *m_updaterHyper));

        m_linSolverObjHyper = std::make_shared<linsolver_t>();
        m_tag = std::make_shared<tag_t>();

        // residual weighting
        std::string basisfile_gpod = m_basisRoot_gpod + "_" + std::to_string(this->m_domIdx) + ".bin";
        m_weigher = std::make_shared<weigh_t>(
            m_weigher_type,
            basisfile_gpod,
            this->m_sampleFile,
            m_nmodes_gpod,
            this->getDofPerCell()
        );

        m_nonlinSolverHyper = std::make_shared<nonlinsolverHyp_t>(
            pressio::create_gauss_newton_solver(
                *m_problemHyper, *m_linSolverObjHyper, *m_weigher, *m_tag
            )
        );

        m_nonlinSolverHyper->setStopCriterion(pressio::nonlinearsolvers::Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance);
        m_nonlinSolverHyper->setStopTolerance(1e-5);
    }

// TODO: to protected
public:
    pressio::ode::StepScheme m_odeScheme;
    std::string m_weigher_type;
    std::string m_basisRoot_gpod;
    int m_nmodes_gpod;
    std::shared_ptr<updaterHyp_t> m_updaterHyper;
    std::shared_ptr<problemHyp_t> m_problemHyper;
    std::shared_ptr<linsolver_t> m_linSolverObjHyper;
    std::shared_ptr<weigh_t> m_weigher;
    std::shared_ptr<tag_t> m_tag;
    std::shared_ptr<nonlinsolverHyp_t> m_nonlinSolverHyper;
};

//
// auxiliary function to create a vector of meshes given a count and meshRoot
//
auto create_meshes(std::string const & meshRoot, const int n)
{
    using mesh_t = pda::cellcentered_uniform_mesh_eigen_type;

    std::vector<mesh_t> meshes;
    std::vector<std::string> meshPaths;

    for (int domIdx = 0; domIdx < n; ++domIdx) {
        // read mesh
        meshPaths.emplace_back(meshRoot + "/domain_" + std::to_string(domIdx));
        meshes.emplace_back( pda::load_cellcentered_uniform_mesh_eigen(meshPaths.back()) );
    }

    return std::tuple(meshes, meshPaths);
}

//
// all domains are assumed to be FOM domains
//
template<class app_t, class mesh_t, class prob_t>
auto create_subdomains(
    const std::vector<mesh_t> & meshes,
    const Tiling & tiling,
    prob_t probId,
    std::vector<pode::StepScheme> & odeSchemes,
    std::vector<pda::InviscidFluxReconstruction> & fluxOrders,
    int icFlag = 0,
    const std::unordered_map<std::string, typename app_t::scalar_type> & userParams = {})
{
    auto ndomains = tiling.count();
    std::vector<std::string> domFlagVec(ndomains, "FOM");

    // dummy arguments
    std::vector<int> nmodesVec(ndomains, -1);
    std::vector<std::string> samplePaths(ndomains, "");
    std::vector<int> nmodesVec_gpod(ndomains, -1);

    return create_subdomains<app_t>(
        meshes, tiling,
        probId, odeSchemes, fluxOrders,
        domFlagVec, "", "", nmodesVec,
        icFlag, "", samplePaths,
        "identity", "", nmodesVec_gpod,
        userParams);

}

//
// Subdomain type specified by domFlagVec
//
template<class app_t, class mesh_t, class prob_t>
auto create_subdomains(
    const std::vector<mesh_t> & meshes,
    const Tiling & tiling,
    prob_t probId,
    std::vector<pode::StepScheme> & odeSchemes,
    std::vector<pda::InviscidFluxReconstruction> & fluxOrders,
    const std::vector<std::string> & domFlagVec,
    const std::string & transRoot,
    const std::string & basisRoot,
    const std::vector<int> & nmodesVec,
    int icFlag = 0,
    const std::string & icFileRoot = "",
    const std::vector<std::string> & samplePaths = {},
    const std::string & weigher_type = "identity",
    const std::string & basisRoot_gpod = "",
    const std::vector<int> & nmodesVec_gpod = {},
    const std::unordered_map<std::string, typename app_t::scalar_type> & userParams = {})
{

    using subdomain_t = SubdomainBase<mesh_t, typename app_t::state_type>;
    std::vector<std::shared_ptr<subdomain_t>> result;

    const int ndomX = tiling.countX();
    const int ndomY = tiling.countY();
    const int ndomZ = tiling.countZ();
    const int ndomains = tiling.count();

    // check sizes
    if (meshes.size() != ndomains) { throw std::runtime_error("Incorrect number of mesh objects"); }
    if (odeSchemes.size() != ndomains) { throw std::runtime_error("Incorrect number of ODE schemes"); }
    if (fluxOrders.size() != ndomains) { throw std::runtime_error("Incorrect number of flux order"); }
    if (domFlagVec.size() != ndomains) { throw std::runtime_error("Incorrect number of domain flags"); }
    if (nmodesVec.size() != ndomains) { throw std::runtime_error("Incorrect number of ROM mode counts"); }
    if (!samplePaths.empty()) {
        if (samplePaths.size() != ndomains) { throw std::runtime_error("Incorrect number of sample mesh paths"); }
    }

    // Gappy POD modes are a bit finicky
    // TODO: generalize to finding substring "Hyper" if Galerkin implemented
    std::vector<int> nmodesVec_gpod_in(ndomains, 0);
    if (std::find(domFlagVec.begin(), domFlagVec.end(), "LSPGHyper") != domFlagVec.end()) {
        if (nmodesVec_gpod.empty()) {
            if (weigher_type != "identity") {
                throw std::runtime_error("Got empty nmodesVec_gpod, must use identity weigher");
            }
        }
        else {
            if (nmodesVec_gpod.size() != ndomains) { throw std::runtime_error("Incorrect number of Gappy POD mode counts"); }
            nmodesVec_gpod_in = nmodesVec_gpod;
        }
    }

    // determine boundary conditions for each subdomain, specify app type
    for (int domIdx = 0; domIdx < ndomains; ++domIdx)
    {

        // the actual BC used are defaulted to Dirichlet, and modified below
        // when they need to be physical BCs
        BCType bcLeft  = BCType::SchwarzDirichlet;
        BCType bcRight = BCType::SchwarzDirichlet;
        BCType bcFront = BCType::SchwarzDirichlet;
        BCType bcBack  = BCType::SchwarzDirichlet;

        const int i = domIdx % ndomX;
        const int j = domIdx / ndomX;

        // left physical boundary
        if (i == 0) {
            bcLeft = getPhysBCs(probId, pda::impl::GhostRelativeLocation::Left);
        }

        // right physical boundary
        if (i == (ndomX - 1)) {
            bcRight = getPhysBCs(probId, pda::impl::GhostRelativeLocation::Right);
        }

        // back physical boundary
        if (j == 0) {
            bcBack = getPhysBCs(probId, pda::impl::GhostRelativeLocation::Back);
        }

        // front physical boundary
        if (j == (ndomY - 1)) {
            bcFront = getPhysBCs(probId, pda::impl::GhostRelativeLocation::Front);
        }

        if (domFlagVec[domIdx] == "FOM") {
            result.emplace_back(std::make_shared<SubdomainFOM<mesh_t, app_t, prob_t>>(
                domIdx, meshes[domIdx],
                bcLeft, bcFront, bcRight, bcBack,
                probId, odeSchemes[domIdx], fluxOrders[domIdx], icFlag, icFileRoot, userParams));
        }
        else if (domFlagVec[domIdx] == "LSPG") {
            result.emplace_back(std::make_shared<SubdomainLSPG<mesh_t, app_t, prob_t>>(
                domIdx, meshes[domIdx],
                bcLeft, bcFront, bcRight, bcBack,
                probId, odeSchemes[domIdx], fluxOrders[domIdx], icFlag, icFileRoot, userParams,
                transRoot, basisRoot, nmodesVec[domIdx]));
        }
        else if (domFlagVec[domIdx] == "LSPGHyper") {
            result.emplace_back(std::make_shared<SubdomainLSPGHyper<mesh_t, app_t, prob_t>>(
                domIdx, meshes[domIdx],
                bcLeft, bcFront, bcRight, bcBack,
                probId, odeSchemes[domIdx], fluxOrders[domIdx], icFlag, icFileRoot, userParams,
                transRoot, basisRoot, nmodesVec[domIdx],
                samplePaths[domIdx],
                weigher_type, basisRoot_gpod, nmodesVec_gpod_in[domIdx]));
        }
        else {
            std::runtime_error("Invalid subdomain flag value: " + domFlagVec[domIdx]);
        }
    }

    return result;
}

}

#endif
