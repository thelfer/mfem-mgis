/*!
 * \file   StateManagerTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   01/04/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "mfem/general/optparser.hpp"
#include "mfem/fem/intrules.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/StateManager.hxx"

const char* mesh_file = nullptr;
int parallel = 0;
int order = 1;

struct StateManagerTest final : public tfel::tests::TestCase {
  StateManagerTest()
      : tfel::tests::TestCase("MFEMMGIS", "StateManagerTest") {
  }  // end of StateManagerTest
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    auto ctx = mfem_mgis::Context{};
    auto or_die = ctx.getFatalFailureHandler();
    auto mesh = mfem_mgis::construct<mfem_mgis::MeshDiscretization>(
                    ctx, mfem_mgis::Parameters{{"MeshFileName", mesh_file},
                                               {"Parallel", bool(parallel)}}) |
                or_die;
    auto sm = mfem_mgis::StateManager{mesh};
    // creating a partial quadrature function for a given finite element
    // discretization and registring it to the state manager
    auto fed1 = mfem_mgis::construct<mfem_mgis::FiniteElementDiscretization>(
                    ctx, mesh,
                    mfem_mgis::Parameters{{{"FiniteElementFamily", "H1"},
                                           {"FiniteElementOrder", order},
                                           {"UnknownsSize", 3}}}) |
                or_die;
    auto qspace1 =
        mfem_mgis::make_shared<mfem_mgis::PartialQuadratureSpace>(
            ctx, fed1, 1,
            [](const mfem::FiniteElement& e, const mfem::ElementTransformation&)
                -> const mfem::IntegrationRule& {
              return mfem::IntRules.Get(e.GetGeomType(), 2);
            }) |
        or_die;
    auto f = mfem_mgis::construct<mfem_mgis::PartialQuadratureFunction>(
                 ctx, qspace1) |
             or_die;
    sm.add(ctx, "function1", f, mfem_mgis::ets) | or_die;
    TFEL_TESTS_ASSERT(sm.contains(ctx, qspace1, "function1", mfem_mgis::ets) |
                      or_die);
    // creating a second finite element discretization and a second quadrature
    // space, which is equivalent to the first one
    auto fed2 = mfem_mgis::construct<mfem_mgis::FiniteElementDiscretization>(
                    ctx, mesh,
                    mfem_mgis::Parameters{{{"FiniteElementFamily", "H1"},
                                           {"FiniteElementOrder", order},
                                           {"UnknownsSize", 2}}}) |
                or_die;
    auto qspace2 =
        mfem_mgis::make_shared<mfem_mgis::PartialQuadratureSpace>(
            ctx, fed2, 1,
            [](const mfem::FiniteElement& e, const mfem::ElementTransformation&)
                -> const mfem::IntegrationRule& {
              return mfem::IntRules.Get(e.GetGeomType(), 2);
            }) |
        or_die;
    // retrieving a view on the function
    const auto oview = sm.get(ctx, qspace2, "function1", mfem_mgis::ets);
    TFEL_TESTS_ASSERT(mfem_mgis::isValid(oview));
    if (!isValid(oview)) {
      std::cerr << ctx.getErrorMessage() << std::endl;
      return;
    }
    // we now check that f and the view points to the values
    TFEL_TESTS_ASSERT(f.getValues().data() == oview->getValues().data());
  }
};

TFEL_TESTS_GENERATE_PROXY(StateManagerTest, "StateManagerTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main(int argc, char** argv) {
  //
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  //
  auto &m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  m.addXMLTestOutput("StateManagerTest.xml");
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
