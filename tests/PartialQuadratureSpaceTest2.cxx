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
#ifdef DO_USE_MPI
#include <mpi.h>
#endif
#include "mfem/general/optparser.hpp"
#include "mfem/fem/intrules.hpp"
#include "mfem/fem/fe/fe_base.hpp"
#include "mfem/fem/eltrans.hpp"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureSpaceIdentifiersManager.hxx"

int main(int argc, char** argv) {
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  //
  const char* mesh_file = nullptr;
  int parallel = 0;
  int order = 1;
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
  auto m = mfem_mgis::construct<mfem_mgis::MeshDiscretization>(
               ctx, mfem_mgis::Parameters{{"MeshFileName", mesh_file},
                                          {"Parallel", bool(parallel)}}) |
           or_die;
  auto fed1 = mfem_mgis::construct<mfem_mgis::FiniteElementDiscretization>(
                  ctx, m,
                  mfem_mgis::Parameters{{{"FiniteElementFamily", "H1"},
                                         {"FiniteElementOrder", order},
                                         {"UnknownsSize", 3}}}) |
              or_die;
  auto fed2 = mfem_mgis::construct<mfem_mgis::FiniteElementDiscretization>(
                  ctx, m,
                  mfem_mgis::Parameters{{{"FiniteElementFamily", "H1"},
                                         {"FiniteElementOrder", order},
                                         {"UnknownsSize", 2}}}) |
              or_die;
  auto qspace1 =
      mfem_mgis::make_shared<mfem_mgis::PartialQuadratureSpace>(
          ctx, fed1, 1,
          [](const mfem::FiniteElement& e, const mfem::ElementTransformation&)
              -> const mfem::IntegrationRule& {
            return mfem::IntRules.Get(e.GetGeomType(), 2);
          }) |
      or_die;
  auto qspace2 =
      mfem_mgis::make_shared<mfem_mgis::PartialQuadratureSpace>(
          ctx, fed2, 1,
          [](const mfem::FiniteElement& e, const mfem::ElementTransformation&)
              -> const mfem::IntegrationRule& {
            return mfem::IntRules.Get(e.GetGeomType(), 2);
          }) |
      or_die;
  assert(mfem_mgis::areEquivalent(*qspace1, *qspace2));
  auto qspace3 =
      mfem_mgis::make_shared<mfem_mgis::PartialQuadratureSpace>(
          ctx, fed1, 1,
          [](const mfem::FiniteElement& e, const mfem::ElementTransformation&)
              -> const mfem::IntegrationRule& {
            return mfem::IntRules.Get(e.GetGeomType(), 4);
          }) |
      or_die;
  assert(!mfem_mgis::areEquivalent(*qspace1, *qspace3));
  assert(!mfem_mgis::areEquivalent(*qspace2, *qspace3));
  //
  auto ids = mfem_mgis::PartialQuadratureSpaceIdentifiersManager{m};
  const auto i1 = ids.getIdentifier(ctx, qspace1) | or_die;
  const auto i1b = ids.getIdentifier(ctx, qspace1) | or_die;
  const auto i2 = ids.getIdentifier(ctx, qspace2) | or_die;
  const auto i3 = ids.getIdentifier(ctx, qspace3) | or_die;
  assert(i1 == 0);
  assert(i1 == i1b);
  assert(i1 == i2);
  assert(i1 != i3);
  assert(i3 == 1);
  return EXIT_SUCCESS;
}
