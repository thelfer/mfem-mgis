/*!
 * \file   src/GridFunctionUtilities.cxx
 * \brief  Implementation of the functions defined in `MFEMMGIS/GridFunctionUtilities.hxx`
 * \author Thomas Helfer
 * \date   03/09/2026
 */

#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/GridFunctionUtilities.hxx"

namespace mfem_mgis {

  template <>
  std::optional<MakeGridFunctionResult<true>>
  makeGridFunction<true>(Context& ctx,
                         const FiniteElementDiscretization& fed,
                         const size_type nc) noexcept {
#ifdef MFEM_USE_MPI
    if (!fed.describesAParallelComputation()) {
      return ctx.registerErrorMessage(
          "can't create a parallel grid function on a finite element "
          "discretization describing a sequential computation");
    }
    return makeGridFunction<true>(ctx, fed.getFiniteElementSpace<true>(), nc);
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of makeGridFunction<true>

  template <>
  std::optional<MakeGridFunctionResult<false>>
  makeGridFunction<false>(Context& ctx,
                         const FiniteElementDiscretization& fed,
                         const size_type nc) noexcept {
    if (fed.describesAParallelComputation()) {
      return ctx.registerErrorMessage(
          "can't create a sequential grid function on a finite element "
          "discretization describing a parallel computation");
    }
    return makeGridFunction<false>(ctx, fed.getFiniteElementSpace<false>(), nc);
  } // end of makeGridFunction<false>

  template <bool parallel>
  std::optional<MakeGridFunctionResult<parallel>> makeGridFunction_impl(
      Context& ctx,
      const FiniteElementSpace<parallel>& fe_space,
      const size_type nc) noexcept {
    if (nc <= 0) {
      return ctx.registerErrorMessage("invalid number of components");
    }
    auto fes = make_unique<FiniteElementSpace<parallel>>(
        ctx, static_cast<Mesh<parallel>*>(fe_space.GetMesh()),
        fe_space.FEColl(), nc, fe_space.GetOrdering());
    if (isInvalid(fes)) {
      return {};
    }
    auto f = make_unique<GridFunction<parallel>>(ctx, fes.get());
    if (isInvalid(f)) {
      return {};
    }
    return MakeGridFunctionResult<parallel>{.fe_space = std::move(fes),
                                            .f = std::move(f)};
  }  // end of makeGridFunction<parallel>

#ifdef MFEM_USE_MPI
  template <>
  std::optional<MakeGridFunctionResult<true>> makeGridFunction<true>(
      Context& ctx,
      const FiniteElementSpace<true>& fe_space,
      const size_type nc) noexcept {
    return makeGridFunction_impl<true>(ctx, fe_space, nc);
  }  // end of makeGridFunction<true>
#else
  template <>
  std::optional<MakeGridFunctionResult<true>> makeGridFunction<true>(
      Context&,
      const FiniteElementSpace<true>&,
      const size_type ) noexcept {
    reportUnsupportedParallelComputations();
  }  // end of makeGridFunction<true>
#endif

  template <>
  std::optional<MakeGridFunctionResult<false>> makeGridFunction<false>(
      Context& ctx,
      const FiniteElementSpace<false>& fe_space,
      const size_type nc) noexcept {
    return makeGridFunction_impl<false>(ctx, fe_space, nc);
  }  // end of makeGridFunction<false>

}  // end of namespace mfem_mgis