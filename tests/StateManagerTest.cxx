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
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/StateManager.hxx"

struct {
  const char* mesh_file = nullptr;
  int parallel = 0;
  int order = 1;
} parameters;

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
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto omesh = construct<MeshDiscretization>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(omesh));
    auto sm = StateManager{*omesh};
    // creating a partial quadrature function for a given finite element
    // discretization and registring it to the state manager
    auto ofed1 = construct<FiniteElementDiscretization>(
        ctx, *omesh,
        dict{{{"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", parameters.order},
              {"UnknownsSize", 3}}});
    TFEL_TESTS_ASSERT(isValid(ofed1));
    auto qspace1 = make_shared<PartialQuadratureSpace>(
        ctx, *ofed1, 1,
        [](const mfem::FiniteElement& e,
           const mfem::ElementTransformation&) -> const mfem::IntegrationRule& {
          return mfem::IntRules.Get(e.GetGeomType(), 2);
        });
    TFEL_TESTS_ASSERT(isValid(qspace1));
    auto of = construct<PartialQuadratureFunction>(ctx, qspace1);
    TFEL_TESTS_ASSERT(isValid(of));
    TFEL_TESTS_CHECK(sm.add(ctx, "function1", *of, ets));
    const auto ook = sm.contains(ctx, qspace1, "function1", ets);
    TFEL_TESTS_ASSERT(isValid(ook));
    TFEL_TESTS_CHECK(*ook);
    // creating a second finite element discretization and a second quadrature
    // space, which is equivalent to the first one
    auto ofed2 = construct<FiniteElementDiscretization>(
                     ctx, *omesh,
                     dict{{{"FiniteElementFamily", "H1"},
                           {"FiniteElementOrder", parameters.order},
                           {"UnknownsSize", 2}}});
    TFEL_TESTS_ASSERT(isValid(ofed2));
    auto qspace2 =
        make_shared<PartialQuadratureSpace>(
            ctx, *ofed2, 1,
            [](const mfem::FiniteElement& e, const mfem::ElementTransformation&)
                -> const mfem::IntegrationRule& {
              return mfem::IntRules.Get(e.GetGeomType(), 2);
            });
    TFEL_TESTS_ASSERT(isValid(qspace2));
    TFEL_TESTS_CHECK(areEquivalent(*qspace1, *qspace2));
    // retrieving a view on the function
    const auto oview = sm.get(ctx, qspace2, "function1", ets);
    TFEL_TESTS_ASSERT(isValid(oview));
    // we now check that f and the view points to the values
    TFEL_TESTS_CHECK(of->getValues().data() == oview->getValues().data());
  }
};

TFEL_TESTS_GENERATE_PROXY(StateManagerTest, "StateManagerTest");

/* coverity [UNCAUGHT_EXCEPT]*/
int main(int argc, char** argv) {
  //
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&parameters.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&parameters.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (parameters.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  if (parameters.parallel) {
    m.addXMLTestOutput("ParallelStateManagerTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("StateManagerTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}
