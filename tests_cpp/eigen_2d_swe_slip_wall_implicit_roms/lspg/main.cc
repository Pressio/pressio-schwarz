
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressiodemoapps/swe2d.hpp"
#include "../../observer.hpp"
#include <chrono>
#include "pressio/rom_subspaces.hpp"
#include "pressio/rom_lspg_unsteady.hpp"
#include "pressio-schwarz/rom_utils.hpp"

int main()
{
    PRESSIOLOG_INITIALIZE(pressiolog::LogLevel::debug);

    namespace pda = pressiodemoapps;
    namespace pode = pressio::ode;
    namespace prom = pressio::rom;
    namespace plspg = pressio::rom::lspg;
    namespace pnlins = pressio::nonlinearsolvers;

    const auto meshObj = pda::load_cellcentered_uniform_mesh_eigen(".");

    std::string transfile = "./trial_space/center.bin";
    std::string basisfile = "./trial_space/basis.bin";
    const int nmodes = 25;

#ifdef USE_WENO5
    const auto order   = pda::InviscidFluxReconstruction::Weno5;
#elif defined USE_WENO3
    const auto order   = pda::InviscidFluxReconstruction::Weno3;
#else
    const auto order   = pda::InviscidFluxReconstruction::FirstOrder;
#endif
    auto scheme = pode::StepScheme::BDF1;
    const int icFlag = 1;

    // create fomSystem
    const auto probId  = pda::Swe2d::SlipWall;
    auto fomSystem = pda::create_problem_eigen(meshObj, probId, order, icFlag);

    using FomSystemType = decltype(fomSystem);
    using scalar_type = FomSystemType::scalar_type;
    using reduced_state_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;

    // read and define trial space
    auto trans = pschwarz::read_vector_from_binary<scalar_type>(transfile);
    auto basis = pschwarz::read_matrix_from_binary<scalar_type>(basisfile, nmodes);
    const auto trialSpace = prom::create_trial_column_subspace<
        reduced_state_type>(std::move(basis), std::move(trans), true);

    // project initial condition
    auto state = fomSystem.initialCondition();
    auto u = pressio::ops::clone(state);
    pressio::ops::update(u, 0., state, 1, trialSpace.translationVector(), -1);
    auto reducedState = trialSpace.createReducedState();
    pressio::ops::product(::pressio::transpose(), 1., trialSpace.basisOfTranslatedSpace(), u, 0., reducedState);

    // define ROM problem
    auto problem = plspg::create_unsteady_problem(scheme, trialSpace, fomSystem);
    auto stepper = problem.lspgStepper();

    // define solver
    using hessian_t       = Eigen::Matrix<scalar_type, -1, -1>;
    using solver_tag      = pressio::linearsolvers::direct::HouseholderQR;
    using linear_solver_t = pressio::linearsolvers::Solver<solver_tag, hessian_t>;
    linear_solver_t linearSolver;

    auto solver = pressio::create_gauss_newton_solver(stepper, linearSolver);
    solver.setStopCriterion(pnlins::Stop::WhenAbsolutel2NormOfGradientBelowTolerance);
    solver.setStopTolerance(1e-5);

    // observer
    StateObserver Obs("swe_slipWall2d_solution.bin", 1);
    RuntimeObserver Obs_run("runtime.bin");

    const double tf = 1.0;
    const double dt = 0.02;
    const auto Nsteps = pressio::ode::StepCount(tf/dt);

    auto runtimeStart = std::chrono::high_resolution_clock::now();
    pode::advance_n_steps(stepper, reducedState, 0.0, dt, Nsteps, Obs, solver);
    auto runtimeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = runtimeEnd - runtimeStart;
    Obs_run(duration.count() * 1e-3);

    PRESSIOLOG_FINALIZE();
    return 0;
}
