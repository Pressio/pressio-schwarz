
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "../../observer.hpp"
#include <chrono>

int main()
{
    PRESSIOLOG_FINALIZE();

    namespace pda = pressiodemoapps;
    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");
#ifdef USE_WENO5
    const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
    const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
    const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif
    const auto probId  = pda::Swe2d::SlipWall;
    int icFlag = 1;

    auto appObj = pda::create_problem_eigen(meshObj, probId, order, icFlag);
    using app_t = decltype(appObj);
    using state_t = typename app_t::state_type;
    using jacob_t = typename app_t::jacobian_type;

    state_t state = appObj.initialCondition();

    auto stepperObj = pressio::ode::create_implicit_stepper(
        pressio::ode::StepScheme::BDF1, appObj);

    using lin_solver_t = pressio::linearsolvers::Solver<
        pressio::linearsolvers::iterative::Bicgstab, jacob_t>;
    lin_solver_t linSolverObj;
    auto NonLinSolver = pressio::create_newton_solver(stepperObj, linSolverObj);
    NonLinSolver.setStopTolerance(1e-5);

    FomObserver<state_t> Obs("swe_slipWall2d_solution.bin", 1);
    RuntimeObserver Obs_run("runtime.bin");

    const double tf = 1.0;
    const double dt = 0.02;
    const auto Nsteps = pressio::ode::StepCount(tf/dt);

    auto runtimeStart = std::chrono::high_resolution_clock::now();
    pressio::ode::advance_n_steps(stepperObj, state, 0., dt, Nsteps, Obs, NonLinSolver);
    auto runtimeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = runtimeEnd - runtimeStart;
    Obs_run(duration.count() * 1e-3);

    PRESSIOLOG_FINALIZE();
    return 0;
}
