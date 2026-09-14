/*!
 * \file   src/GridFunctionUtilities.cxx
 * \brief  Implementation of the functions defined in
 * `MFEMMGIS/GridFunctionUtilities.hxx` \author Thomas Helfer \date   03/09/2026
 */

#include "mfem/fem/gridfunc.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pgridfunc.hpp"
#endif /* MFEM_USE_MPI */

#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/GridFunctionUtilities.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  size_type getNumberOfComponents(const GridFunction<true>& f) noexcept {
    return f.ParFESpace()->GetVDim();
  }  // end of getNumberOfComponents

#endif /* MFEM_USE_MPI */

  size_type getNumberOfComponents(const GridFunction<false>& f) noexcept {
    return f.FESpace()->GetVDim();
  }  // end of getNumberOfComponents

  template <>
  std::unique_ptr<GridFunction<true>> makeGridFunction<true>(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type nc) noexcept {
#ifdef MFEM_USE_MPI
    if (!fed.describesAParallelComputation()) {
      return ctx.registerErrorMessage(
          "can't create a parallel grid function on a finite element "
          "discretization describing a sequential computation");
    }
    auto m = fed.getFiniteElementSpacesManager();
    auto fes = m.getFiniteElementSpace<true>(ctx, nc);
    if (isInvalid(fes)) {
      return {};
    }
    return make_unique<GridFunction<true>>(ctx, fes.get());
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of makeGridFunction<true>

  template <>
  std::unique_ptr<GridFunction<false>> makeGridFunction<false>(
      Context& ctx,
      const FiniteElementDiscretization& fed,
      const size_type nc) noexcept {
    if (fed.describesAParallelComputation()) {
      return ctx.registerErrorMessage(
          "can't create a sequential grid function on a finite element "
          "discretization describing a parallel computation");
    }
    auto m = fed.getFiniteElementSpacesManager();
    auto fes = m.getFiniteElementSpace<false>(ctx, nc);
    if (isInvalid(fes)) {
      return {};
    }
    return make_unique<GridFunction<false>>(ctx, fes.get());
  }  // end of makeGridFunction<false>

}  // end of namespace mfem_mgis