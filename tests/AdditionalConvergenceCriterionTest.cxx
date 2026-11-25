/*!
 * \file   tests/AdditionalConvergenceCriterionTest.cxx
 * \brief  Tests of the `AbstractAdditionalConvergenceCriterion` interface
 *         across `NonLinearEvolutionProblem`, `AbstractNonLinearSolver`,
 *         `NewtonSolver`, and `Simulation`.
 * \date   25/11/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>
#include "mfem/general/optparser.hpp"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearSolvers/AbstractAdditionalConvergenceCriterion.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/Simulation.hxx"

namespace {

  struct TestOptions {
    const char* mesh_file = nullptr;
    const char* library = nullptr;
    int order = 1;
    int parallel = 0;
  };

  TestOptions test_options;

  //! \brief helper criterion that records lifecycle calls and enforces a
  //! minimum number of check calls per solve
  struct RecordingVetoCriterion final
      : mfem_mgis::AbstractAdditionalConvergenceCriterion {
    explicit RecordingVetoCriterion(
        const int min_checks, std::vector<std::string>& event_log) noexcept
        : required_checks(min_checks), events(event_log) {}

    void helper() noexcept override { this->events.emplace_back("helper"); }

    void reset() noexcept override {
      this->checks_in_step = 0;
      this->events.emplace_back("reset");
    }

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context&, const CheckArguments& args) noexcept override {
      this->events.push_back("check:" + std::to_string(args.iter) + ":" +
                             (args.converged ? "1" : "0"));
      this->last_ref_norm = args.reference_residual_norm;
      this->last_max_iter = args.max_iter;
      this->last_u_size = args.u.Size();
      ++(this->checks_in_step);
      if (!args.converged) {
        return false;
      }
      return this->checks_in_step >= this->required_checks;
    }

    ~RecordingVetoCriterion() noexcept override = default;

    int required_checks;
    int checks_in_step = 0;
    mfem_mgis::real last_ref_norm = mfem_mgis::real{0};
    int last_max_iter = 0;
    int last_u_size = 0;
    std::vector<std::string>& events;
  };

  //! \brief continuation criterion adjusting a spatially-indexed ShearModulus
  //! field between prediction (nominal), iteration 0 (stiffened), and converged
  //! intermediate state (relaxed back to nominal)
  struct ShearModulusContinuationCriterion final
      : mfem_mgis::AbstractAdditionalConvergenceCriterion {
    ShearModulusContinuationCriterion(std::vector<mfem_mgis::real>& s0_field,
                                      std::vector<mfem_mgis::real>& s1_field,
                                      const mfem_mgis::real nominal_mu,
                                      const mfem_mgis::real stiff_mu) noexcept
        : mu0(s0_field),
          mu1(s1_field),
          mu_nominal(nominal_mu),
          mu_stiff(stiff_mu) {}

    void setShearModulus(const mfem_mgis::real val) noexcept {
      std::fill(this->mu0.begin(), this->mu0.end(), val);
      std::fill(this->mu1.begin(), this->mu1.end(), val);
    }

    void helper() noexcept override {
      this->helper_called = true;
      this->setShearModulus(this->mu_nominal);
    }

    void reset() noexcept override {
      this->reset_called = true;
      this->relaxed = false;
      this->setShearModulus(this->mu_stiff);
    }

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context&, const CheckArguments& args) noexcept override {
      if (args.iter == 0) {
        this->iter0_converged = args.converged;
      }
      if (!args.converged) {
        return false;
      }
      if (!this->relaxed) {
        this->setShearModulus(this->mu_nominal);
        this->relaxed = true;
        return false;
      }
      return true;
    }

    ~ShearModulusContinuationCriterion() noexcept override = default;

    std::vector<mfem_mgis::real>& mu0;
    std::vector<mfem_mgis::real>& mu1;
    mfem_mgis::real mu_nominal;
    mfem_mgis::real mu_stiff;
    bool helper_called = false;
    bool reset_called = false;
    bool relaxed = false;
    bool iter0_converged = true;
  };

  //! \brief criterion updating the linear solver tolerance during iterations
  struct LinearToleranceContinuationCriterion final
      : mfem_mgis::AbstractAdditionalConvergenceCriterion {
    LinearToleranceContinuationCriterion(
        mfem_mgis::NonLinearEvolutionProblem& prob,
        mfem_mgis::Context& c,
        const mfem_mgis::real target) noexcept
        : problem(prob), ctx(c), target_tol(target) {}

    void applyLinearSolver(const mfem_mgis::real tol) noexcept {
      const auto ok = this->problem.setLinearSolver(
          this->ctx, "CGSolver",
          {{"VerbosityLevel", 0},
           {"AbsoluteTolerance", mfem_mgis::real{0}},
           {"RelativeTolerance", tol},
           {"MaximumNumberOfIterations", 1000}});
      if (!ok) {
        this->linear_solver_update_failed = true;
      }
    }

    void helper() noexcept override {
      this->applyLinearSolver(this->target_tol);
    }

    void reset() noexcept override {
      this->tightened = false;
      this->applyLinearSolver(mfem_mgis::real{1e-1});
    }

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context&, const CheckArguments& args) noexcept override {
      if (!args.converged) {
        return false;
      }
      if (!this->tightened) {
        this->applyLinearSolver(this->target_tol);
        this->tightened = true;
        return false;
      }
      return true;
    }

    ~LinearToleranceContinuationCriterion() noexcept override = default;

    mfem_mgis::NonLinearEvolutionProblem& problem;
    mfem_mgis::Context& ctx;
    mfem_mgis::real target_tol;
    bool tightened = false;
    bool linear_solver_update_failed = false;
  };

  //! \brief criterion checking call count and iteration threshold
  struct IterationThresholdCriterion final
      : mfem_mgis::AbstractAdditionalConvergenceCriterion {
    explicit IterationThresholdCriterion(const int min_iter) noexcept
        : target_iter(min_iter) {}

    void helper() noexcept override {}

    void reset() noexcept override { this->calls = 0; }

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context&, const CheckArguments& args) noexcept override {
      ++(this->calls);
      return args.iter >= this->target_iter;
    }

    ~IterationThresholdCriterion() noexcept override = default;

    int target_iter;
    int calls = 0;
  };

  //! \brief criterion that unconditionally returns true
  struct AlwaysTrueCriterion final
      : mfem_mgis::AbstractAdditionalConvergenceCriterion {
    void helper() noexcept override {}

    void reset() noexcept override {}

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context&, const CheckArguments& args) noexcept override {
      if (args.iter == 0) {
        this->saw_unconverged_iter0 = !args.converged;
      }
      ++(this->calls);
      return true;
    }

    ~AlwaysTrueCriterion() noexcept override = default;

    bool saw_unconverged_iter0 = false;
    int calls = 0;
  };

  //! \brief criterion that aborts by returning an invalid optional
  struct FailingCriterion final
      : mfem_mgis::AbstractAdditionalConvergenceCriterion {
    void helper() noexcept override {}

    void reset() noexcept override {}

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context& c, const CheckArguments&) noexcept override {
      return c.registerErrorMessage("custom additional convergence failure");
    }

    ~FailingCriterion() noexcept override = default;
  };

  //! \brief criterion that rejects only the first solve attempt to trigger
  //! substepping in Simulation
  struct RejectFirstAttemptCriterion final
      : mfem_mgis::AbstractAdditionalConvergenceCriterion {
    void helper() noexcept override { ++(this->helper_calls); }

    void reset() noexcept override { ++(this->reset_calls); }

    [[nodiscard]] std::optional<bool> check(
        mfem_mgis::Context& c, const CheckArguments& args) noexcept override {
      if (this->reset_calls == 1) {
        return c.registerErrorMessage("initial step rejected for substepping");
      }
      return args.converged;
    }

    ~RejectFirstAttemptCriterion() noexcept override = default;

    int helper_calls = 0;
    int reset_calls = 0;
  };

}  // namespace

struct AdditionalConvergenceCriterionTest final : public tfel::tests::TestCase {
  AdditionalConvergenceCriterionTest()
      : tfel::tests::TestCase("MFEMMGIS",
                              "AdditionalConvergenceCriterionTest") {}

  tfel::tests::TestResult execute() override {
    this->testLifecycleAndVeto();
    this->testParameterContinuation();
    this->testLinearSolverUpdate();
    this->testMultipleCriteriaAndNoShortCircuit();
    this->testCannotBypassResidualConvergence();
    this->testErrorReporting();
    this->testSubSteppingWithSimulation();
    return this->result;
  }

 private:
  //! \brief build an isotropic linear elastic problem on the unit cube
  std::unique_ptr<mfem_mgis::NonLinearEvolutionProblem> buildElasticProblem(
      mfem_mgis::Context& ctx,
      const mfem_mgis::PredictionStrategy strategy,
      const bool set_uniform_shear = true) {
    using namespace mfem_mgis;
    auto oproblem = construct<NonLinearEvolutionProblem>(
        ctx, Parameters{{"MeshFileName", test_options.mesh_file},
                        {"FiniteElementFamily", "H1"},
                        {"FiniteElementOrder", test_options.order},
                        {"UnknownsSize", 3},
                        {"NumberOfUniformRefinements",
                         test_options.parallel ? 1 : 0},
                        {"Hypothesis", "Tridimensional"},
                        {"Parallel", test_options.parallel != 0}});
    TFEL_TESTS_ASSERT(isValid(oproblem));
    auto& problem = *oproblem;
    TFEL_TESTS_ASSERT(
        isValid(problem.addBehaviourIntegrator(ctx, "Mechanics", 1,
                                               test_options.library,
                                               "IsotropicLinearElasticity")));
    auto om1 = problem.getMaterial(ctx, 1, 0);
    TFEL_TESTS_ASSERT(isValid(om1));
    auto& mat = *om1;
    for (auto* ps : {&mat.s0, &mat.s1}) {
      TFEL_TESTS_ASSERT(mgis::behaviour::setMaterialProperty(
          ctx, *ps, "FirstLameCoefficient", real{100e9}));
      if (set_uniform_shear) {
        TFEL_TESTS_ASSERT(mgis::behaviour::setMaterialProperty(
            ctx, *ps, "ShearModulus", real{75e9}));
      }
      TFEL_TESTS_ASSERT(mgis::behaviour::setExternalStateVariable(
          ctx, *ps, "Temperature", real{293.15}));
    }
    TFEL_TESTS_ASSERT(problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", 1}, {"Component", 1}}));
    TFEL_TESTS_ASSERT(problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", 2}, {"Component", 2}}));
    TFEL_TESTS_ASSERT(problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", 5}, {"Component", 0}}));
    TFEL_TESTS_ASSERT(problem.addUniformDirichletBoundaryCondition(
        ctx, {{"Boundary", 3},
              {"Component", 0},
              {"LoadingEvolution",
               [](const auto t) noexcept { return real{3e-2} * t; }}}));
    problem.setPredictionPolicy({.strategy = strategy});
    TFEL_TESTS_ASSERT(problem.setLinearSolver(
        ctx, "CGSolver",
        {{"VerbosityLevel", 0},
         {"AbsoluteTolerance", real{1e-16}},
         {"RelativeTolerance", real{1e-16}},
         {"MaximumNumberOfIterations", 1000}}));
    TFEL_TESTS_ASSERT(problem.setSolverParameters(
        ctx, {{"VerbosityLevel", 0},
              {"RelativeTolerance", real{1e-12}},
              {"AbsoluteTolerance", real{0}},
              {"MaximumNumberOfIterations", 20}}));
    return std::make_unique<NonLinearEvolutionProblem>(std::move(problem));
  }

  //! \brief verify helper -> prediction -> reset -> check call order and veto
  void testLifecycleAndVeto() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto problem = this->buildElasticProblem(
        ctx, PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION);
    auto events = std::vector<std::string>{};
    auto crit = std::make_shared<RecordingVetoCriterion>(2, events);
    problem->getSolver().addAdditionalConvergenceCriterion(crit);

    const auto r1 = problem->solve(ctx, 0, 1);
    TFEL_TESTS_CHECK(r1.status);
    TFEL_TESTS_CHECK_EQUAL(r1.iterations, 1);
    TFEL_TESTS_CHECK(crit->last_ref_norm > real{0});
    TFEL_TESTS_CHECK_EQUAL(crit->last_max_iter, 20);
    TFEL_TESTS_CHECK(crit->last_u_size > 0);
    const auto expected1 =
        std::vector<std::string>{"helper", "reset", "check:0:1", "check:1:1"};
    TFEL_TESTS_CHECK(events == expected1);

    TFEL_TESTS_ASSERT(problem->update(ctx));
    const auto r2 = problem->solve(ctx, 1, 1);
    TFEL_TESTS_CHECK(r2.status);
    TFEL_TESTS_CHECK_EQUAL(r2.iterations, 1);
    const auto expected2 = std::vector<std::string>{
        "helper", "reset", "check:0:1", "check:1:1",
        "helper", "reset", "check:0:1", "check:1:1"};
    TFEL_TESTS_CHECK(events == expected2);
  }

  //! \brief verify parameter continuation across helper, reset, and check
  void testParameterContinuation() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto ref_problem = this->buildElasticProblem(
        ctx, PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION);
    const auto ref_out = ref_problem->solve(ctx, 0, 1);
    TFEL_TESTS_ASSERT(ref_out.status);
    const auto& u_ref = ref_problem->getUnknowns(ets);
    auto oref_mat = ref_problem->getMaterial(ctx, 1, 0);
    TFEL_TESTS_ASSERT(isValid(oref_mat));

    auto problem = this->buildElasticProblem(
        ctx, PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION, false);
    auto omat = problem->getMaterial(ctx, 1, 0);
    TFEL_TESTS_ASSERT(isValid(omat));
    auto& mat = *omat;
    TFEL_TESTS_ASSERT(mat.n > 1);
    auto mu0 = std::vector<real>(static_cast<std::size_t>(mat.n), real{75e9});
    auto mu1 = std::vector<real>(static_cast<std::size_t>(mat.n), real{75e9});
    TFEL_TESTS_ASSERT(mgis::behaviour::setMaterialProperty(
        ctx, mat.s0, "ShearModulus", std::span<real>(mu0),
        mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE));
    TFEL_TESTS_ASSERT(mgis::behaviour::setMaterialProperty(
        ctx, mat.s1, "ShearModulus", std::span<real>(mu1),
        mgis::behaviour::MaterialStateManager::EXTERNAL_STORAGE));

    auto crit = std::make_shared<ShearModulusContinuationCriterion>(
        mu0, mu1, real{75e9}, real{150e9});
    problem->getSolver().addAdditionalConvergenceCriterion(crit);

    const auto out = problem->solve(ctx, 0, 1);
    TFEL_TESTS_CHECK(out.status);
    TFEL_TESTS_CHECK(crit->helper_called);
    TFEL_TESTS_CHECK(crit->reset_called);
    TFEL_TESTS_CHECK(!crit->iter0_converged);
    TFEL_TESTS_CHECK(out.iterations >= 2);

    const auto& u = problem->getUnknowns(ets);
    TFEL_TESTS_ASSERT(u.Size() == u_ref.Size());
    for (size_type i = 0; i != u.Size(); ++i) {
      TFEL_TESTS_CHECK(std::abs(u[i] - u_ref[i]) < real{1e-11});
    }
    if ((mat.n > 0) && (oref_mat->n > 0)) {
      const auto s = mat.s1.thermodynamic_forces[0];
      const auto s_ref = oref_mat->s1.thermodynamic_forces[0];
      TFEL_TESTS_CHECK(std::abs(s - s_ref) < real{1e-10} * std::abs(s_ref));
    }
  }

  //! \brief verify updating the linear solver from within check()
  void testLinearSolverUpdate() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto ref_problem = this->buildElasticProblem(
        ctx, PredictionStrategy::DEFAULT_PREDICTION);
    const auto ref_out = ref_problem->solve(ctx, 0, 1);
    TFEL_TESTS_ASSERT(ref_out.status);
    const auto& u_ref = ref_problem->getUnknowns(ets);

    auto problem = this->buildElasticProblem(
        ctx, PredictionStrategy::DEFAULT_PREDICTION);
    TFEL_TESTS_ASSERT(problem->setSolverParameters(
        ctx, {{"VerbosityLevel", 0},
              {"RelativeTolerance", real{1e-12}},
              {"AbsoluteTolerance", real{0}},
              {"MaximumNumberOfIterations", 100}}));
    auto crit = std::make_shared<LinearToleranceContinuationCriterion>(
        *problem, ctx, real{1e-16});
    problem->getSolver().addAdditionalConvergenceCriterion(crit);

    const auto out = problem->solve(ctx, 0, 1);
    TFEL_TESTS_CHECK(out.status);
    TFEL_TESTS_CHECK(!crit->linear_solver_update_failed);
    TFEL_TESTS_CHECK(crit->tightened);
    TFEL_TESTS_CHECK(out.iterations >= 2);

    const auto& u = problem->getUnknowns(ets);
    TFEL_TESTS_ASSERT(u.Size() == u_ref.Size());
    for (size_type i = 0; i != u.Size(); ++i) {
      TFEL_TESTS_CHECK(std::abs(u[i] - u_ref[i]) < real{1e-11});
    }
  }

  //! \brief verify multiple criteria are all evaluated without short-circuiting
  void testMultipleCriteriaAndNoShortCircuit() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto problem = this->buildElasticProblem(
        ctx, PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION);
    auto c1 = std::make_shared<IterationThresholdCriterion>(1);
    auto c2 = std::make_shared<IterationThresholdCriterion>(2);
    problem->getSolver().addAdditionalConvergenceCriterion(c1);
    problem->getSolver().addAdditionalConvergenceCriterion(c2);

    const auto out = problem->solve(ctx, 0, 1);
    TFEL_TESTS_CHECK(out.status);
    TFEL_TESTS_CHECK_EQUAL(out.iterations, 2);
    TFEL_TESTS_CHECK_EQUAL(c1->calls, 3);
    TFEL_TESTS_CHECK_EQUAL(c2->calls, 3);
  }

  //! \brief verify a returning true criterion cannot force convergence
  void testCannotBypassResidualConvergence() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto problem = this->buildElasticProblem(
        ctx, PredictionStrategy::DEFAULT_PREDICTION);
    auto crit = std::make_shared<AlwaysTrueCriterion>();
    problem->getSolver().addAdditionalConvergenceCriterion(crit);

    const auto out = problem->solve(ctx, 0, 1);
    TFEL_TESTS_CHECK(out.status);
    TFEL_TESTS_CHECK(crit->saw_unconverged_iter0);
    TFEL_TESTS_CHECK(out.iterations >= 1);
    TFEL_TESTS_CHECK(crit->calls >= 2);
  }

  //! \brief verify returning std::nullopt aborts resolution and preserves error
  void testErrorReporting() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto problem = this->buildElasticProblem(
        ctx, PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION);
    auto failing = std::make_shared<FailingCriterion>();
    auto trailing = std::make_shared<IterationThresholdCriterion>(0);
    problem->getSolver().addAdditionalConvergenceCriterion(failing);
    problem->getSolver().addAdditionalConvergenceCriterion(trailing);

    const auto out = problem->solve(ctx, 0, 1);
    TFEL_TESTS_CHECK(!out.status);
    TFEL_TESTS_CHECK_EQUAL(trailing->calls, 0);
    const auto msg = ctx.getErrorMessage();
    TFEL_TESTS_CHECK(msg.find("custom additional convergence failure") !=
                     std::string::npos);
  }

  //! \brief verify Simulation substepping re-invokes helper() and reset()
  void testSubSteppingWithSimulation() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto problem = this->buildElasticProblem(
        ctx, PredictionStrategy::BEGINNING_OF_TIME_STEP_PREDICTION);
    auto crit = std::make_shared<RejectFirstAttemptCriterion>();
    problem->getSolver().addAdditionalConvergenceCriterion(crit);

    auto sim = construct<Simulation>(
        ctx, *problem,
        Parameters{{"Times", std::vector<Parameter>{real{0}, real{1}}},
                   {"KeepOutputs", true}});
    TFEL_TESTS_ASSERT(isValid(sim));

    const auto [status, output] = sim->run(ctx);
    TFEL_TESTS_CHECK(status == ExitStatus::success);
    TFEL_TESTS_ASSERT(isValid(output));
    TFEL_TESTS_CHECK_EQUAL(crit->helper_calls, 3);
    TFEL_TESTS_CHECK_EQUAL(crit->reset_calls, 3);
    const auto nsteps = get<size_type>(ctx, *output, "NumberOfTimeSteps");
    const auto nsubsteps = get<size_type>(ctx, *output, "NumberOfSubSteps");
    TFEL_TESTS_ASSERT(isValid(nsteps));
    TFEL_TESTS_ASSERT(isValid(nsubsteps));
    TFEL_TESTS_CHECK_EQUAL(*nsteps, 2);
    TFEL_TESTS_CHECK_EQUAL(*nsubsteps, 1);
  }
};

TFEL_TESTS_GENERATE_PROXY(AdditionalConvergenceCriterionTest,
                          "AdditionalConvergenceCriterionTest");

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&test_options.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&test_options.library, "-l", "--library",
                 "Material library.");
  args.AddOption(&test_options.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&test_options.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (test_options.mesh_file == nullptr) ||
      (test_options.library == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  if (test_options.parallel != 0) {
    m.addXMLTestOutput(
        "ParallelAdditionalConvergenceCriterionTest-" +
        std::to_string(mfem_mgis::getMPISizeForDefaultCommunicator()) + "-" +
        std::to_string(mfem_mgis::getMPIRankForDefaultCommunicator()) + ".xml");
  } else {
    m.addXMLTestOutput("AdditionalConvergenceCriterionTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
