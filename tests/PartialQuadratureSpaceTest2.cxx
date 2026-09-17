/*!
 * \file   tests/PartialQuadratureSpaceTest2.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/03/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <memory>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "mfem/fem/intrules.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureSpaceIdentifiersManager.hxx"

struct {
  const char* mesh_file = nullptr;
  int parallel = 0;
  int order = 1;
} parameters;

struct PartialQuadratureSpaceTest2 final : public tfel::tests::TestCase {
  PartialQuadratureSpaceTest2()
      : tfel::tests::TestCase("MFEMMGIS", "PartialQuadratureSpaceTest2") {
  }  // end of PartialQuadratureSpaceTest2
  tfel::tests::TestResult execute() override {
    this->test1();
    return this->result;
  }

 private:
  void test1() {
    using namespace mfem_mgis;
    auto ctx = Context{};
    auto om = construct<MeshDiscretization>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(om));
    auto ofed1 = construct<FiniteElementDiscretization>(
        ctx, *om,
        dict{{{"FiniteElementFamily", "H1"},
              {"FiniteElementOrder", parameters.order},
              {"UnknownsSize", 3}}});
    TFEL_TESTS_ASSERT(isValid(ofed1));
    auto ofed2 = construct<FiniteElementDiscretization>(
        ctx, ofed1->getFiniteElementSpacesManager(),
        dict{{{"UnknownsSize", 2}}});
    TFEL_TESTS_ASSERT(isValid(ofed2));
    auto qspace1 = make_shared<PartialQuadratureSpace>(
        ctx, *ofed1, 1,
        [](const mfem::FiniteElement& e,
           const mfem::ElementTransformation&) -> const mfem::IntegrationRule& {
          return mfem::IntRules.Get(e.GetGeomType(), 2);
        });
    TFEL_TESTS_ASSERT(isValid(qspace1));
    auto qspace2 = make_shared<PartialQuadratureSpace>(
        ctx, *ofed2, 1,
        [](const mfem::FiniteElement& e,
           const mfem::ElementTransformation&) -> const mfem::IntegrationRule& {
          return mfem::IntRules.Get(e.GetGeomType(), 2);
        });
    TFEL_TESTS_ASSERT(isValid(qspace2));
    TFEL_TESTS_CHECK(areEquivalent(*qspace1, *qspace2));
    auto qspace3 = make_shared<PartialQuadratureSpace>(
        ctx, *ofed1, 1,
        [](const mfem::FiniteElement& e,
           const mfem::ElementTransformation&) -> const mfem::IntegrationRule& {
          return mfem::IntRules.Get(e.GetGeomType(), 4);
        });
    TFEL_TESTS_ASSERT(isValid(qspace3));
    TFEL_TESTS_CHECK(!areEquivalent(*qspace1, *qspace3));
    TFEL_TESTS_CHECK(!areEquivalent(*qspace2, *qspace3));
    //
    auto oids = construct<PartialQuadratureSpaceIdentifiersManager>(ctx, *om);
    TFEL_TESTS_ASSERT(isValid(oids));
    const auto oi1 = oids->getIdentifier(ctx, qspace1);
    TFEL_TESTS_ASSERT(isValid(oi1));
    const auto oi1b = oids->getIdentifier(ctx, qspace1);
    TFEL_TESTS_ASSERT(isValid(oi1b));
    const auto oi2 = oids->getIdentifier(ctx, qspace2);
    TFEL_TESTS_ASSERT(isValid(oi2));
    const auto oi3 = oids->getIdentifier(ctx, qspace3);
    TFEL_TESTS_ASSERT(isValid(oi3));
    TFEL_TESTS_CHECK(*oi1 == 0);
    TFEL_TESTS_CHECK(*oi1 == *oi1b);
    TFEL_TESTS_CHECK(*oi1 == *oi2);
    TFEL_TESTS_CHECK(*oi1 != *oi3);
    TFEL_TESTS_CHECK(*oi3 == 1);
  }
};

TFEL_TESTS_GENERATE_PROXY(PartialQuadratureSpaceTest2,
                          "PartialQuadratureSpaceTest2");

int main(int argc, char** argv) {
  //
  mfem_mgis::initialize(argc, argv);
  // options treatment
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
  m.addTestOutput(std::cout);  // mfem_mgis::getDefaultLogStream());
  if (parameters.parallel) {
    m.addXMLTestOutput("ParallelPartialQuadratureSpaceTest2-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("PartialQuadratureSpaceTest2.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
