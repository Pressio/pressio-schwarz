#include <chrono>
#include "pressiodemoapps/swe2d.hpp"
#include "pressio-schwarz/schwarz.hpp"
#include "../../../observer.hpp"
#include "../../../help_cmdline.hpp"

int main(int argc, char *argv[])
{
    namespace pda  = pressiodemoapps;
    namespace pode = pressio::ode;

#if defined SCHWARZ_ENABLE_THREADPOOL
    const int numthreads = parse_num_threads(argc, argv);
    std::string dir_suffix = "_tp";
#elif defined SCHWARZ_ENABLE_OMP
    std::string dir_suffix = "_omp";
#else
    std::string dir_suffix = "";
#endif

    // +++++ USER INPUTS +++++

    const int obsFreq = 1;

    // problem definition
    const auto probId = pda::Swe2d::CustomBCs;
#ifdef USE_WENO5
    static_assert(false);
    std::vector<pda::InviscidFluxReconstruction> orderVec(12, pda::InviscidFluxReconstruction::Weno5);
#elif defined USE_WENO3
    std::string outRoot = "./weno3" + dir_suffix;
    std::vector<pda::InviscidFluxReconstruction> orderVec(12, pda::InviscidFluxReconstruction::Weno3);
#else
    std::string outRoot = "./firstorder" + dir_suffix;
    std::vector<pda::InviscidFluxReconstruction> orderVec(12, pda::InviscidFluxReconstruction::FirstOrder);
#endif
    std::string obsRoot = outRoot + "/swe_slipWall2d_solution";
    std::string meshRootFull = outRoot + "/full_mesh_decomp";
    std::string meshRootHyper = outRoot + "/sample_mesh_decomp";

    std::vector<pode::StepScheme> schemeVec(12, pode::StepScheme::BDF1);
    const int icFlag = 1;
    using app_t = pschwarz::swe2d_app_type;

    // ROM definition
    std::vector<std::string> domFlagVec(12, "LSPGHyper");
    std::string transRoot = outRoot + "/trial_space/center";
    std::string basisRoot = outRoot + "/trial_space/basis";
    std::vector<int> nmodesVec(12, 25);

    // time stepping
    const double tf = 1.0;
    std::vector<double> dt(1, 0.02);
    const int convergeStepMax = 10;
    const double abs_err_tol = 1e-11;
    const double rel_err_tol = 1e-11;

    // +++++ END USER INPUTS +++++

    auto tiling = std::make_shared<pschwarz::Tiling>(meshRootFull);
    auto [meshObjsFull, meshPathsFull] = pschwarz::create_meshes(meshRootFull, tiling->count());
    std::vector<std::string> samplePaths;
    for (int domIdx = 0; domIdx < meshPathsFull.size(); ++ domIdx) {
        samplePaths.emplace_back(meshRootHyper + "/domain_" + std::to_string(domIdx) + "/sample_mesh_gids.dat");
    }
    auto subdomains = pschwarz::create_subdomains<app_t>(
        meshObjsFull, *tiling, probId, schemeVec, orderVec,
        domFlagVec, transRoot, basisRoot, nmodesVec, icFlag, "", samplePaths);
    pschwarz::SchwarzDecomp decomp(subdomains, tiling, dt);

    // observers
    using state_t = decltype(decomp)::state_t;
    using obs_t = FomObserver<state_t>;
    std::vector<obs_t> obsVec((*decomp.m_tiling).count());
    for (int domIdx = 0; domIdx < (*decomp.m_tiling).count(); ++domIdx) {
        obsVec[domIdx] = obs_t(obsRoot + "_" + std::to_string(domIdx) + ".bin", obsFreq);
        obsVec[domIdx](::pressio::ode::StepCount(0), 0.0, *decomp.m_subdomainVec[domIdx]->getStateFull());
    }
    RuntimeObserver obs_time(outRoot + "/runtime.bin");

// -----------------------------------------
// OMP
// -----------------------------------------
#if defined SCHWARZ_ENABLE_OMP

    const int numSteps = tf / decomp.m_dtMax;
    std::chrono::time_point<std::chrono::high_resolution_clock> runtimeStart;
    double secsElapsed;

#pragma omp parallel firstprivate(numSteps, rel_err_tol, abs_err_tol, convergeStepMax)
{

    double simultime = 0.0;
    for (int outerStep = 1; outerStep <= numSteps; ++outerStep)
    {
#pragma omp single
        {
            std::cout << "Step " << outerStep << std::endl;
        }

        auto runtimeStart = std::chrono::high_resolution_clock::now();
        auto numSubiters = decomp.additive_step(outerStep, simultime, rel_err_tol, abs_err_tol, convergeStepMax);
        simultime += decomp.m_dtMax;
#pragma omp barrier
#pragma omp master
        {
            const auto runtimeEnd = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> duration = runtimeEnd - runtimeStart;
            obs_time(duration.count() * 1e-3, numSubiters);

            // output observer
            if ((outerStep % obsFreq) == 0) {
                const auto stepWrap = pode::StepCount(outerStep);
                for (int domIdx = 0; domIdx < (*decomp.m_tiling).count(); ++domIdx) {
                    obsVec[domIdx](stepWrap, time, *decomp.m_subdomainVec[domIdx]->getStateFull());
                }
            }
        }

    }
}

// -----------------------------------------
// Thread pool
// -----------------------------------------
#elif defined SCHWARZ_ENABLE_THREADPOOL

    // solve
    BS::thread_pool pool(numthreads);
    const int numSteps = tf / decomp.m_dtMax;
    double time = 0.0;
    for (int outerStep = 1; outerStep <= numSteps; ++outerStep)
    {
        std::cout << "Step " << outerStep << std::endl;

        auto runtimeStart = std::chrono::high_resolution_clock::now();
        auto numSubiters = decomp.additive_step(outerStep, time, rel_err_tol,
                             abs_err_tol, convergeStepMax, pool);
        const auto runtimeEnd = std::chrono::high_resolution_clock::now();
        const auto nsDuration = std::chrono::duration_cast<std::chrono::nanoseconds>(runtimeEnd - runtimeStart);
        const double secsElapsed = static_cast<double>(nsDuration.count()) * 1e-9;

        time += decomp.m_dtMax;

        if ((outerStep % obsFreq) == 0) {
            const auto stepWrap = pode::StepCount(outerStep);
            for (int domIdx = 0; domIdx < (*decomp.m_tiling).count(); ++domIdx) {
                obsVec[domIdx](stepWrap, time, *decomp.m_subdomainVec[domIdx]->getStateFull());
            }
        }
        obs_time(secsElapsed, numSubiters);
    }

#endif

    return 0;
}
