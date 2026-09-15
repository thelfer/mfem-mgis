/*!
 * \file   tests/MeshDiscretizationTest.cxx
 * \brief
 * \author Thomas Helfer
 * \date   14/09/2026
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "TFEL/Tests/TestCase.hxx"
#include "TFEL/Tests/TestProxy.hxx"
#include "TFEL/Tests/TestManager.hxx"
#include "mfem/general/optparser.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"

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

struct MeshDiscretizationTest final : public tfel::tests::TestCase {
  MeshDiscretizationTest()
      : tfel::tests::TestCase("MFEMMGIS", "MeshDiscretizationTest") {
  }  // end of MeshDiscretizationTest
  tfel::tests::TestResult execute() override {
#ifdef MFEM_USE_MPI
    if (parameters.parallel == 1) {
      this->test1<true>();
    }
#endif /* MFEM_USE_MPI */
    if (parameters.parallel == 0) {
      this->test1<false>();
    }
    return this->result;
  }

 private:
  template <bool parallel>
  void test1() {
    using namespace mfem_mgis;
    auto omesh = construct<MeshDiscretization>(
        ctx,
        Parameters{{"MeshFileName", parameters.mesh_file},
                   {"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
                   {"Materials", dict{{"Attr1", 1}, {"Attr2", 2}}},
                   {"Boundaries", dict{{"left", 1}, {"right", 2}}},
                   {"Parallel", bool(parameters.parallel)}});
    TFEL_TESTS_ASSERT(isValid(omesh));
    TFEL_TESTS_CHECK(omesh->manages(omesh->template getMesh<parallel>()));
    const auto osm1a = omesh->template getSubMesh<parallel>(
        ctx, "Attr1", MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(osm1a));
    TFEL_TESTS_CHECK(omesh->manages(*osm1a));
    const auto ook1 = omesh->isDefinedOnMaterials(ctx, *osm1a);
    TFEL_TESTS_ASSERT(isValid(ook1));
    TFEL_TESTS_CHECK(*ook1);
    const auto ook2 = omesh->isDefinedOnBoundaries(ctx, *osm1a);
    TFEL_TESTS_ASSERT(isValid(ook2));
    TFEL_TESTS_CHECK(!(*ook2));
    const auto& attr1 = osm1a->attributes;
    TFEL_TESTS_ASSERT(attr1.Size() == 1);
    TFEL_TESTS_CHECK_EQUAL(attr1[0], 1);
    //
    const auto osm1b = omesh->template getSubMesh<parallel>(
        ctx, 1, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isValid(osm1b));
    TFEL_TESTS_CHECK(omesh->manages(*osm1b));
    TFEL_TESTS_CHECK(&(*osm1a) == &(*osm1b));
    //
    const auto osm2 = omesh->template getSubMesh<parallel>(
        ctx, "Attr2", MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_CHECK(omesh->manages(*osm2));
    TFEL_TESTS_CHECK(&(*osm1a) != &(*osm2));
    const auto& attr2 = osm2->attributes;
    TFEL_TESTS_ASSERT(attr2.Size() == 1);
    TFEL_TESTS_CHECK_EQUAL(attr2[0], 2);
    const auto osm3 = omesh->template getSubMesh<parallel>(
        ctx, list{"Attr1", "Attr2"},
        MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_CHECK(isInvalid(osm3));
    if (parameters.parallel) {
      TFEL_TESTS_CHECK_EQUAL(
          ctx.getRawErrorMessage(),
          "can't create a parallel sub mesh on all materials");
    } else {
      TFEL_TESTS_CHECK_EQUAL(
          ctx.getRawErrorMessage(),
          "can't create a sequential sub mesh on all materials");
    }
    //
    const auto osm4 = omesh->template getSubMesh<parallel>(
        ctx, "Att1", MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isInvalid(osm4));
    TFEL_TESTS_CHECK_EQUAL(ctx.getRawErrorMessage(),
                           "getMaterialsIdentifiers: no material matching "
                           "regular expression 'Att1'");
    //
    const auto osm5 = omesh->template getSubMesh<parallel>(
        ctx, list{-1, -2, -3}, MeshDiscretization::Location::ON_MATERIALS);
    TFEL_TESTS_ASSERT(isInvalid(osm5));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "getMaterialsIdentifiers: no material associated with attribute '-1'");
    //
    const auto obsm1 = omesh->template getSubMesh<parallel>(
        ctx, list{"left", "right"},
        MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_ASSERT(isValid(obsm1));
    TFEL_TESTS_CHECK(omesh->manages(*obsm1));
    const auto ook3 = omesh->isDefinedOnMaterials(ctx, *obsm1);
    TFEL_TESTS_ASSERT(isValid(ook3));
    TFEL_TESTS_CHECK(!(*ook3));
    const auto ook4 = omesh->isDefinedOnBoundaries(ctx, *obsm1);
    TFEL_TESTS_ASSERT(isValid(ook4));
    TFEL_TESTS_CHECK(*ook4);
    const auto& attr3 = obsm1->attributes;
    TFEL_TESTS_ASSERT(attr3.Size() == 2);
    TFEL_TESTS_CHECK_EQUAL(attr3[0], 1);
    TFEL_TESTS_CHECK_EQUAL(attr3[1], 2);
    // check that attributes are indeed sorted
    const auto obsm2 = omesh->template getSubMesh<parallel>(
        ctx, list{"right", "left"},
        MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_ASSERT(isValid(obsm2));
    TFEL_TESTS_CHECK(omesh->manages(*obsm2));
    TFEL_TESTS_CHECK_EQUAL(&(*obsm1), &(*obsm2));
    // an attribute can't appear twice or more
    const auto obsm3 = omesh->template getSubMesh<parallel>(
        ctx, list{"right", "right"},
        MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_CHECK(isInvalid(obsm3));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "getBoundariesIdentifiers: identifier '2' multiply selected");
    // bad boundary attribute
    const auto obsm4 = omesh->template getSubMesh<parallel>(
        ctx, -1, MeshDiscretization::Location::ON_BOUNDARIES);
    TFEL_TESTS_CHECK(isInvalid(obsm4));
    TFEL_TESTS_CHECK_EQUAL(
        ctx.getRawErrorMessage(),
        "getBoundariesIdentifiers: no boundary associated with attribute '-1'");
    // creating a sub mesh not managed by the mesh description
    auto mids = mfem::Array<size_type>{1};
    mids[0] = 1;
    auto sm = SubMesh<parallel>::CreateFromDomain(
        *(omesh->template getMutableMeshPointer<parallel>()), mids);
    TFEL_TESTS_CHECK(!omesh->manages(sm));
  }  // end of test1
};

TFEL_TESTS_GENERATE_PROXY(MeshDiscretizationTest, "MeshDiscretizationTest");

int main(int argc, char** argv) {
  // options treatment
  mfem_mgis::initialize(argc, argv);
  parseCommandLineOptions(argc, argv);
  //
  auto& m = tfel::tests::TestManager::getTestManager();
  m.addTestOutput(std::cout);
  if (parameters.parallel == 1) {
    m.addXMLTestOutput("ParallelMeshDiscretizationTest-" +
                       std::to_string(mfem_mgis::getMPIsize()) + "-" +
                       std::to_string(mfem_mgis::getMPIrank()) + ".xml");
  } else {
    m.addXMLTestOutput("MeshDiscretizationTest.xml");
  }
  return m.execute().success() ? EXIT_SUCCESS : EXIT_FAILURE;
}  // end of main
