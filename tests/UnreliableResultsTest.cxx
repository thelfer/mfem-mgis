/*!
 * \file   tests/UnreliableResultsTest.cxx
 * \brief  This test checks that a simulation whose time steps return
 * unreliable results returns the `unreliableResults` status
 * \date   01/10/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

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

//! \brief a model whose time steps always return unreliable results
struct UnreliableModel final : public mfem_mgis::ModelBase {
  UnreliableModel(mfem_mgis::Context& ctx,
                  const mfem_mgis::MeshDiscretization& m) noexcept
      : ModelBase(ctx, m) {}
  std::string getName() const noexcept override { return "UnreliableModel"; }
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
    r.first.update(mfem_mgis::ExitStatus::unreliableResults);
    return r;
  }

 protected:
  using mfem_mgis::ModelBase::addPostProcessing;
};

struct UnreliableResultsTest final : public tfel::tests::TestCase {
  UnreliableResultsTest()
      : tfel::tests::TestCase("MFEMMGIS", "UnreliableResultsTest") {
  }  // end of UnreliableResultsTest
  tfel::tests::TestResult execute() override {
    this->check(false);
    this->check(true);
    return this->result;
  }

 private:
  /*!
   * \brief run a simulation and check its status
   * \param[in] keep_outputs: value of the `KeepOutputs` parameter
   */
  void check(const bool keep_outputs) {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto mesh = construct<MeshDiscretization>(
        ctx, Parameters{{"MeshFileName", mesh_file}, {"Parallel", false}});
    TFEL_TESTS_ASSERT(isValid(mesh));
    auto ps = construct<PhysicalSystem>(ctx, *mesh);
    TFEL_TESTS_ASSERT(isValid(ps));
    auto c = make_shared<LoopCouplingScheme>(
        ctx, *mesh, Parameters{{"NumberOfIterations", 1}});
    TFEL_TESTS_ASSERT(isValid(c));
    TFEL_TESTS_ASSERT(
        c->addModel(ctx, std::make_shared<UnreliableModel>(ctx, *mesh)));
    TFEL_TESTS_ASSERT(ps->setCouplingScheme(ctx, c));
    auto s = construct<Simulation>(
        ctx, *ps,
        Parameters{{"Times", std::vector<Parameter>{0., 0.5, 1.}},
                   {"KeepOutputs", keep_outputs}});
    TFEL_TESTS_ASSERT(isValid(s));
    TFEL_TESTS_CHECK(s->run(ctx).first == ExitStatus::unreliableResults);
  }  // end of check
};

TFEL_TESTS_GENERATE_PROXY(UnreliableResultsTest, "UnreliableResultsTest");

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
  m.addXMLTestOutput("UnreliableResultsTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
