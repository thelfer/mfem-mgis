/*!
 * \file   MFEMMGIS/FiniteElementSpacesManager.ixx
 * \brief  This file implements the inline methods of the
 * `FiniteElementSpacesManager` class
 * \author Thomas Helfer
 * \date 09/09/2026
 */

#ifndef LIB_MFEM_MGIS_FINITEELEMENTSPACESMANAGER_IXX
#define LIB_MFEM_MGIS_FINITEELEMENTSPACESMANAGER_IXX

namespace mfem_mgis {

  template <bool parallel>
  std::shared_ptr<FiniteElementSpace<parallel>>
  FiniteElementSpacesManager::getFiniteElementSpace(
      Context& ctx, const size_type nc) const noexcept {
    if constexpr (parallel) {
      return this->getParallelFiniteElementSpace(ctx, nc);
    } else {
      return this->getSequentialFiniteElementSpace(ctx, nc);
    }
  }  // end of getFiniteElementSpace

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FINITEELEMENTSPACESMANAGER_IXX */
