
#include "pressio/ode_steppers.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/advection_diffusion2d.hpp"
#include "../observer.hpp"
#include <chrono>

int main()
{
    PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

    namespace pda = pressiodemoapps;
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
    const auto schemeInv   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
    const auto schemeInv   = pda::InviscidFluxReconstruction::Weno3;
#else
    const auto schemeInv   = pda::InviscidFluxReconstruction::FirstOrder;
#endif

    const auto schemeVisc = pda::ViscousFluxReconstruction::FirstOrder;

    const auto probId  = pda::AdvectionDiffusion2d::BurgersOutflow;

    auto appObj = pda::create_problem_eigen(
        meshObj,
        probId,
        schemeInv,
        schemeVisc
    );
    using app_t = decltype(appObj);
    using state_t = typename app_t::state_type;
    using jacob_t = typename app_t::jacobian_type;

    state_t state = appObj.initialCondition();

    auto stepperObj = pressio::ode::create_implicit_stepper(
        pressio::ode::StepScheme::BDF1, appObj);

    using lin_solver_t = pressio::linsol::Solver<
        pressio::linsol::iterative::Bicgstab, jacob_t>;
    lin_solver_t linSolverObj;
    auto NonLinSolver = pressio::nlsol::create_newton_solver(stepperObj, linSolverObj);
    NonLinSolver.setStopTolerance(1e-5);

    FomObserver<state_t> Obs("burgers_outflow2d_solution.bin", 1);
    RuntimeObserver Obs_run("runtime.bin");

    const double tf = 1.0;
    const double dt = 0.02;
    const auto Nsteps = pressio::ode::StepCount(tf/dt);

    auto runtimeStart = std::chrono::high_resolution_clock::now();
    auto policy = pressio::ode::steps_fixed_dt(0., Nsteps, dt);
    pressio::ode::advance(stepperObj, state, policy, NonLinSolver, Obs);
    auto runtimeEnd = std::chrono::high_resolution_clock::now();
    auto nsElapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(runtimeEnd - runtimeStart).count();
    double secElapsed = static_cast<double>(nsElapsed) * 1e-9;

    Obs_run(secElapsed);

    PRESSIOLOG_FINALIZE();
    return 0;
}
