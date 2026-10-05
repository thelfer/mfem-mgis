/*!
 * \file   tests/FiniteElementDiscretizationTest.cxx
 * \brief  Tests of the `FiniteElementDiscretization` and
 * `FiniteElementSpacesManager` classes
 * \author Thomas Helfer
 * \date   14/09/2026
 */

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

struct TestParameters {
  const char* mesh_file = nullptr;
  int parallel = 0;
};  // end of struct TestParameters

auto parameters = TestParameters{};
auto ctx = mgis::Context{};

static void parseCommandLineOptions(int argc, char** argv) {
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&parameters.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&parameters.parallel, "-p", "--parallel",
                 "choose between serial (-p 0) and parallel (-p 1)");
  args.Parse();
  if ((!args.Good()) || (parameters.mesh_file == nullptr)) {
    args.PrintUsage(mfem_mgis::getOutputStream());
    mfem_mgis::abort(EXIT_FAILURE);
  }
}  // end of parseCommandLineOptions

struct FiniteElementDiscretizationTest final : public tfel::tests::TestCase {
  FiniteElementDiscretizationTest()
      : tfel::tests::TestCase("MFEMMGIS", "FiniteElementDiscretizationTest") {
  }  // end of FiniteElementDiscretizationTest
  tfel::tests::TestResult execute() override {
#ifdef MFEM_USE_MPI
    if (parameters.parallel == 1) {
      this->test1<true>();
      this->test2<true>();
      this->test3<true>();
      this->test4<true>();
    }
#endif /* MFEM_USE_MPI */
    if (parameters.parallel == 0) {
      this->test1<false>();
      this->test2<false>();
      this->test3<false>();
      this->test4<false>();
    }
    return this->result;
  }

 private:
  template <bool parallel>
  void test1() {
    using namespace mfem_mgis;
    auto ofed = construct<FiniteElementDiscretization>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"FiniteElementFamily", "H1"},
                  {"FiniteElementOrder", 2},
                  {"UnknownsSize", 2},
                  {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                  {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                  {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(ofed));
    auto& m = ofed->template getMesh<parallel>();
    auto& fes = ofed->template getFiniteElementSpace<parallel>();
    auto fespaces = ofed->getFiniteElementSpacesManager();
    //
    TFEL_TESTS_CHECK(fespaces.manages(fes));
    //
    this->template common_tests_on_materials<parallel>(fespaces);
    this->template common_tests_on_boundaries<parallel>(fespaces);
    // check that if a finite element space is created on all materials, the
    // initial one is returned.
    auto fes1 = fespaces.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = list{"Attr1", "Attr2"},
              .number_of_components = 2});
    TFEL_TESTS_ASSERT(isValid(fes1));
    if constexpr (parallel) {
      TFEL_TESTS_CHECK_EQUAL(fes1->GetParMesh(), &m);
    } else {
      TFEL_TESTS_CHECK_EQUAL(fes1->GetMesh(), &m);
    }
    TFEL_TESTS_CHECK_EQUAL(fes1.get(), &fes);
    TFEL_TESTS_CHECK_EQUAL(fes1->FEColl(),
                           ofed->getFiniteElementCollectionPointer().get());
    // check that if a finite element space is created on all materials with a
    // different number of components, a new finite element space is created
    auto fes2 = fespaces.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = list{"Attr1", "Attr2"},
              .number_of_components = 3});
    TFEL_TESTS_ASSERT(isValid(fes2));
    TFEL_TESTS_CHECK(fes2.get() != &fes);
    TFEL_TESTS_CHECK_EQUAL(fes2->FEColl(),
                           ofed->getFiniteElementCollectionPointer().get());
    // check reuse of fespaces defined on all materials
    auto fes3 = fespaces.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = list{2, 1},
              .number_of_components = 3});
    TFEL_TESTS_ASSERT(isValid(fes3));
    TFEL_TESTS_CHECK(fes2.get() == fes3.get());
    TFEL_TESTS_CHECK_EQUAL(fes3->FEColl(),
                           ofed->getFiniteElementCollectionPointer().get());
  }  // end of test1
  template <bool parallel>
  void test2() {
    using namespace mfem_mgis;
    auto om = construct<FiniteElementSpacesManager>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"FiniteElementFamily", "H1"},
                  {"FiniteElementOrder", 2},
                  {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                  {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                  {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(om));
    this->template common_tests_on_materials<parallel>(*om);
    this->template common_tests_on_boundaries<parallel>(*om);
  }
  template <bool parallel>
  void test3() {
    using namespace mfem_mgis;
    auto omesh = construct<MeshDiscretization>(
        ctx, dict{{"MeshFileName", parameters.mesh_file},
                  {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                  {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                  {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                  {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(omesh));
    auto om = construct<FiniteElementSpacesManager>(
        ctx, *omesh,
        dict{{"FiniteElementFamily", "H1"}, {"FiniteElementOrder", 2}});
    this->template common_tests_on_materials<parallel>(*om);
    this->template common_tests_on_boundaries<parallel>(*om);
  }
  //! \brief tests of the `FiniteElementSpaceOrdering` parameter
  template <bool parallel>
  void test4() {
    using namespace mfem_mgis;
    const auto mesh_parameters =
        dict{{"MeshFileName", parameters.mesh_file},
             {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
             {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
             {"Boundaries", dict{{"left", 1}, {"right", 2}}},
             {"Parallel", bool(parameters.parallel)}};
    // mesh parameters completed by the given parameters
    auto merge = [&mesh_parameters](const Parameters& p) {
      auto r = mesh_parameters;
      r.insert(throwing, p);
      return r;
    };
    // check the ordering of all the finite element spaces of a manager
    auto check = [this](const FiniteElementSpacesManager& m,
                        const mfem::Ordering::Type o) {
      const auto ofes1 = m.template getFiniteElementSpace<parallel>(ctx, 3);
      TFEL_TESTS_ASSERT(isValid(ofes1));
      TFEL_TESTS_CHECK_EQUAL(ofes1->GetOrdering(), o);
      const auto ofes2 = m.template getFiniteElementSpace<parallel>(
          ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
                .identifiers = list{"Attr1"},
                .number_of_components = 3});
      TFEL_TESTS_ASSERT(isValid(ofes2));
      TFEL_TESTS_CHECK_EQUAL(ofes2->GetOrdering(), o);
      const auto ofes3 = m.template getFiniteElementSpace<parallel>(
          ctx, {.location = MeshDiscretization::Location::ON_BOUNDARIES,
                .identifiers = list{"left"},
                .number_of_components = 2});
      TFEL_TESTS_ASSERT(isValid(ofes3));
      TFEL_TESTS_CHECK_EQUAL(ofes3->GetOrdering(), o);
    };
    // default ordering
    auto om1 = construct<FiniteElementSpacesManager>(
        ctx, merge(dict{{"FiniteElementOrder", 1}}));
    TFEL_TESTS_ASSERT(isValid(om1));
    check(*om1, mfem::Ordering::byNODES);
    // explicit orderings
    for (const auto& [n, o] : {std::pair{"byNODES", mfem::Ordering::byNODES},
                               std::pair{"byVDIM", mfem::Ordering::byVDIM}}) {
      auto om2 = construct<FiniteElementSpacesManager>(
          ctx, merge(dict{{"FiniteElementOrder", 1},
                          {"FiniteElementSpaceOrdering", n}}));
      TFEL_TESTS_ASSERT(isValid(om2));
      check(*om2, o);
    }
    // finite element space manager built on an existing mesh
    auto omesh = construct<MeshDiscretization>(ctx, mesh_parameters);
    TFEL_TESTS_ASSERT(isValid(omesh));
    auto om3 = construct<FiniteElementSpacesManager>(
        ctx, *omesh,
        dict{{"FiniteElementOrder", 1},
             {"FiniteElementSpaceOrdering", "byVDIM"}});
    TFEL_TESTS_ASSERT(isValid(om3));
    check(*om3, mfem::Ordering::byVDIM);
    // finite element discretization
    auto ofed = construct<FiniteElementDiscretization>(
        ctx, merge(dict{{"FiniteElementOrder", 1},
                        {"FiniteElementSpaceOrdering", "byVDIM"},
                        {"UnknownsSize", 3}}));
    TFEL_TESTS_ASSERT(isValid(ofed));
    TFEL_TESTS_CHECK_EQUAL(
        ofed->template getFiniteElementSpace<parallel>().GetOrdering(),
        mfem::Ordering::byVDIM);
    check(ofed->getFiniteElementSpacesManager(), mfem::Ordering::byVDIM);
    // invalid ordering
    auto om4 = construct<FiniteElementSpacesManager>(
        ctx, *omesh,
        dict{{"FiniteElementOrder", 1},
             {"FiniteElementSpaceOrdering", "byComponents"}});
    TFEL_TESTS_CHECK(isInvalid(om4));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "FiniteElementSpacesManager::FiniteElementSpacesManager: invalid "
        "finite element space ordering 'byComponents' (expected 'byNODES' or "
        "'byVDIM')");
  }  // end of test4
  //
  template <bool parallel>
  void common_tests_on_materials(
      const mfem_mgis::FiniteElementSpacesManager& m) {
    using namespace mfem_mgis;
    const auto mesh = m.getMeshDiscretization();
    const auto osm1 = mesh.template getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(osm1));
    auto fes1 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = list{"Attr1"},
              .number_of_components = 3});
    TFEL_TESTS_ASSERT(isValid(fes1));
    TFEL_TESTS_CHECK(m.manages(*fes1));
    if constexpr (parallel) {
      TFEL_TESTS_CHECK(mesh.manages(*(fes1->GetParMesh())));
      TFEL_TESTS_CHECK_EQUAL(fes1->GetParMesh(), &(*osm1));
    } else {
      TFEL_TESTS_CHECK(mesh.manages(*(fes1->GetMesh())));
      TFEL_TESTS_CHECK_EQUAL(fes1->GetMesh(), &(*osm1));
    }
    TFEL_TESTS_CHECK_EQUAL(fes1->GetVDim(), 3);
    TFEL_TESTS_CHECK_EQUAL(fes1->FEColl(),
                           m.getFiniteElementCollectionPointer().get());
    // check that we return the same finite element space as before
    auto fes2 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = 1,
              .number_of_components = 3});
    TFEL_TESTS_ASSERT(isValid(fes2));
    TFEL_TESTS_CHECK(m.manages(*fes2));
    TFEL_TESTS_CHECK_EQUAL(fes1.get(), fes2.get());
    // invalid number of components
    auto fes3 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = 1,
              .number_of_components = -3});
    TFEL_TESTS_CHECK(isInvalid(fes3));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "invalid number of components ('-3')");
    // invalid identifer
    auto fes4 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = -1,
              .number_of_components = 3});
    TFEL_TESTS_CHECK(isInvalid(fes4));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "getMaterialsIdentifiers: no material associated with attribute '-1'");
    // empty list of materials
    auto fes5 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_MATERIALS,
              .identifiers = list{},
              .number_of_components = 3});
    TFEL_TESTS_CHECK(isInvalid(fes5));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "getMaterialsIdentifiers: empty list of identifiers");
  }  // end of common_tests_on_materials
  template <bool parallel>
  void common_tests_on_boundaries(
      const mfem_mgis::FiniteElementSpacesManager& m) {
    using namespace mfem_mgis;
    const auto mesh = m.getMeshDiscretization();
    const auto osm1 = mesh.template getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_ASSERT(isValid(osm1));
    auto fes1 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_BOUNDARIES,
              .identifiers = list{"left"},
              .number_of_components = 3});
    TFEL_TESTS_ASSERT(isValid(fes1));
    TFEL_TESTS_CHECK(m.manages(*fes1));
    if constexpr (parallel) {
      TFEL_TESTS_CHECK(mesh.manages(*(fes1->GetParMesh())));
      const auto ook = mesh.isDefinedOnBoundaries(ctx, *(fes1->GetParMesh()));
      TFEL_TESTS_ASSERT(isValid(ook));
      TFEL_TESTS_CHECK(*ook);
      TFEL_TESTS_CHECK_EQUAL(fes1->GetParMesh(), &(*osm1));
    } else {
      TFEL_TESTS_CHECK(mesh.manages(*(fes1->GetMesh())));
      const auto ook = mesh.isDefinedOnBoundaries(ctx, *(fes1->GetMesh()));
      TFEL_TESTS_ASSERT(isValid(ook));
      TFEL_TESTS_CHECK(*ook);
      TFEL_TESTS_CHECK_EQUAL(fes1->GetMesh(), &(*osm1));
    }
    TFEL_TESTS_CHECK_EQUAL(fes1->GetVDim(), 3);
    TFEL_TESTS_CHECK_EQUAL(fes1->FEColl(),
                           m.getFiniteElementCollectionPointer().get());
    // check that we return the same finite element space as before
    auto fes2 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_BOUNDARIES,
              .identifiers = 1,
              .number_of_components = 3});
    TFEL_TESTS_ASSERT(isValid(fes2));
    TFEL_TESTS_CHECK(m.manages(*fes2));
    TFEL_TESTS_CHECK_EQUAL(fes1.get(), fes2.get());
    // invalid number of components
    auto fes3 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_BOUNDARIES,
              .identifiers = 1,
              .number_of_components = -3});
    TFEL_TESTS_CHECK(isInvalid(fes3));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "invalid number of components ('-3')");
    // invalid identifer
    auto fes4 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_BOUNDARIES,
              .identifiers = -1,
              .number_of_components = 3});
    TFEL_TESTS_CHECK(isInvalid(fes4));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "getBoundariesIdentifiers: no boundary associated with attribute '-1'");
    // empty list of materials
    auto fes5 = m.template getFiniteElementSpace<parallel>(
        ctx, {.location = MeshDiscretization::Location::ON_BOUNDARIES,
              .identifiers = list{},
              .number_of_components = 3});
    TFEL_TESTS_CHECK(isInvalid(fes5));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "getBoundariesIdentifiers: empty list of identifiers");
  }  // end of common_tests_on_boundaries
};

TFEL_TESTS_GENERATE_PROXY(FiniteElementDiscretizationTest,
                          "FiniteElementDiscretizationTest");

int main(int argc, char** argv) {
  // options treatment
  mfem_mgis::initialize(argc, argv);
  parseCommandLineOptions(argc, argv);
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  if (parameters.parallel == 1) {
    m.addXMLTestOutput("ParallelFiniteElementDiscretizationTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("FiniteElementDiscretizationTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
