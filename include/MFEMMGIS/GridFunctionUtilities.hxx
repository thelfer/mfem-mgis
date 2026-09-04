/*!
 * \file   MFEMMGIS/GridFunctionUtilities.hxx
 * \brief  Various utility functions around MFEM's `GridFunction`
 * \author Thomas Helfer
 * \date   03/09/2026
 */

#ifndef LIB_MFEMMGIS_GRIDFUNCTIONUTILITIES_HXX
#define LIB_MFEMMGIS_GRIDFUNCTIONUTILITIES_HXX

#include <memory>
#include <optional>

#include "mfem/fem/gridfunc.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pgridfunc.hpp"
#endif /* MFEM_USE_MPI */

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  //! \brief result `makeGridFunction`
  template <bool parallel>
  struct MakeGridFunctionResult {
    //! \brief finite element space on which the grid function is defined
    std::unique_ptr<FiniteElementSpace<parallel>> fe_space;
    //! \brief created grid function
    std::unique_ptr<GridFunction<parallel>> f;
  };  // end of MakeGridFunctionResult

  /*!
   * \return a GridFunction with the given number of components, creating a new
   * finite element space if required.
   *
   * In MFEM, a GridFunction has the number of components (VDIM) of the
   * underlying finite element space, which is quite limiting in practice.
   *
   * \param[in] ctx: execution context
   * \param[in] fed: finite element discretization
   * \param[in] nc: number of components
   *
   */
  template <bool parallel>
  [[nodiscard]] std::optional<MakeGridFunctionResult<parallel>>
  makeGridFunction(Context&,
                   const FiniteElementDiscretization&,
                   const size_type) noexcept;

  // partial specialisations
  template <>
  [[nodiscard]] std::optional<MakeGridFunctionResult<true>>
  makeGridFunction<true>(Context&,
                         const FiniteElementDiscretization&,
                         const size_type) noexcept;

  template <>
  [[nodiscard]] std::optional<MakeGridFunctionResult<false>>
  makeGridFunction<false>(Context&,
                          const FiniteElementDiscretization&,
                          const size_type) noexcept;

  /*!
   * \return a GridFunction with the given number of components, creating a new
   * finite element space if required.
   *
   * In MFEM, a GridFunction has the number of components (VDIM) of the
   * underlying finite element space, which is quite limiting in practice.
   *
   * \param[in] ctx: execution context
   * \param[in] fed: finite element discretization
   * \param[in] nc: number of components
   *
   */
  template <bool parallel>
  [[nodiscard]] std::optional<MakeGridFunctionResult<parallel>>
  makeGridFunction(Context&,
                   const FiniteElementSpace<parallel>&,
                   const size_type) noexcept;

  // partial specialisations
  template <>
  [[nodiscard]] std::optional<MakeGridFunctionResult<true>>
  makeGridFunction<true>(Context&,
                         const FiniteElementSpace<true>&,
                         const size_type) noexcept;
  template <>
  [[nodiscard]] std::optional<MakeGridFunctionResult<false>>
  makeGridFunction<false>(Context&,
                          const FiniteElementSpace<false>&,
                          const size_type) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_GRIDFUNCTIONUTILITIES_HXX */
