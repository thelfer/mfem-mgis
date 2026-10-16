/*!
 * \file   tests/SubSteppingMessageTest.cxx
 * \brief  This test checks that, when a time step fails and is cut, the
 * message of the simulation describes the failed time step, not the
 * sub-step which replaces it
 * \date   16/10/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <sstream>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/ModelBase.hxx"
#include "MFEMMGIS/LoopCouplingScheme.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/Simulation.hxx"

const char* mesh_file = nullptr;

//! \brief a model whose first time step fails
struct FailingOnceModel final : public mfem_mgis::ModelBase {
  FailingOnceModel(mfem_mgis::Context& ctx,
                   const mfem_mgis::MeshDiscretization& m) noexcept
      : ModelBase(ctx, m) {}
  std::string getName() const noexcept override { return "FailingOnceModel"; }
  bool addPostProcessing(mfem_mgis::Context& ctx,
                         std::string_view,
                         const mfem_mgis::Parameters&) noexcept override {
    return ctx.registerErrorMessage("no post-processing available");
  }
  std::vector<std::string> getAvailablePostProcessings()
      const noexcept override {
    return {};
  }
  std::pair<mfem_mgis::ExitStatus,
            std::optional<mfem_mgis::ComputeNextStateOutput>>
  computeNextState(mfem_mgis::Context& ctx,
                   const mfem_mgis::TimeStep& ts) noexcept override {
    auto r = ModelBase::computeNextState(ctx, ts);
    if (this->first) {
      this->first = false;
      r.first.update(mfem_mgis::ExitStatus::recoverableError);
    }
    return r;
  }

 protected:
  using mfem_mgis::ModelBase::addPostProcessing;

 private:
  //! \brief if true, the next time step is the first one
  bool first = true;
};

struct SubSteppingMessageTest final : public tfel::tests::TestCase {
  SubSteppingMessageTest()
      : tfel::tests::TestCase("MFEMMGIS", "SubSteppingMessageTest") {
  }  // end of SubSteppingMessageTest
  tfel::tests::TestResult execute() override {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto log = std::make_shared<std::ostringstream>();
    ctx.setLogStream(log);
    ctx.setVerbosityLevel(mfem_mgis::verboseLevel1);
    auto mesh = construct<MeshDiscretization>(
        ctx, Parameters{{"MeshFileName", mesh_file}, {"Parallel", false}});
    TFEL_TESTS_ASSERT(isValid(mesh));
    auto ps = construct<PhysicalSystem>(ctx, *mesh);
    TFEL_TESTS_ASSERT(isValid(ps));
    auto c = make_shared<LoopCouplingScheme>(
        ctx, *mesh, Parameters{{"NumberOfIterations", 1}});
    TFEL_TESTS_ASSERT(isValid(c));
    TFEL_TESTS_ASSERT(
        c->addModel(ctx, std::make_shared<FailingOnceModel>(ctx, *mesh)));
    TFEL_TESTS_ASSERT(ps->setCouplingScheme(ctx, c));
    auto s = construct<Simulation>(
        ctx, *ps, Parameters{{"Times", std::vector<Parameter>{0., 1.}}});
    TFEL_TESTS_ASSERT(isValid(s));
    TFEL_TESTS_CHECK(s->run(ctx).first.shallContinue());
    // the failed time step is the whole temporal sequence, from 0 to 1,
    // then cut in two
    const auto m = log->str();
    TFEL_TESTS_CHECK(m.find("did not converge for the time step from "
                            "0.000000 to 1.000000") != std::string::npos);
    TFEL_TESTS_CHECK(m.find("did not converge for the time step from "
                            "0.000000 to 0.500000") == std::string::npos);
    return this->result;
  }
};

TFEL_TESTS_GENERATE_PROXY(SubSteppingMessageTest, "SubSteppingMessageTest");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  //
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.Parse();
  if (!args.Good()) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  //
  if (mesh_file == nullptr) {
    mfem_mgis::getOutputStream() << "no mesh file specified\n";
    return EXIT_FAILURE;
  }
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("SubSteppingMessageTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
