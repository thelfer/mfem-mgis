/*!
 * \file   include/MFEMMGIS/BoundaryUtilities.ixx
 * \brief  This file implements the inline functions declared in
 * `MFEMMGIS/BoundaryUtilities.hxx`
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_BOUNDARYUTILITIES_IXX
#define LIB_MFEM_MGIS_BOUNDARYUTILITIES_IXX

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  inline std::vector<std::pair<size_type, size_type>> buildFacesDescription(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const size_type bid) {
    std::vector<std::pair<size_type, size_type>> faces;
    auto& fed = p.getFiniteElementDiscretization();
    auto& m = fed.template getMesh<parallel>();
    auto& fes = fed.template getFiniteElementSpace<parallel>();
    for (size_type i = 0; i < fes.GetNBE(); i++) {
      const auto bdr_attr = m.GetBdrAttribute(i);
      if (bdr_attr != bid) {
        continue;
      }
      const auto* const tr = m.GetBdrFaceTransformations(i);
      if (tr != nullptr) {
        faces.push_back({i, tr->Elem1No});
      }
    }
    return faces;
  }  // end of buildFacesDescription

  template <bool parallel>
  inline std::vector<std::pair<size_type, std::vector<std::vector<size_type>>>>
  getElementsDegreesOfFreedomOnBoundary(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const size_type bid) {
    auto& fed = p.getFiniteElementDiscretization();
    auto& m = fed.template getMesh<parallel>();
    auto& fes = fed.template getFiniteElementSpace<parallel>();
    const auto nc = fes.GetVDim();
    auto r = std::vector<
        std::pair<size_type, std::vector<std::vector<size_type>>>>{};
    // the boundary attributes are known by all processes
    if ((m.bdr_attributes.Size() == 0) || (bid < 1) ||
        (bid > m.bdr_attributes.Max())) {
      return r;
    }
    // in parallel, a boundary element is only given to the process owning the
    // adjacent element: GetEssentialVDofs synchronizes the marked degrees of
    // freedom between processes, so that an element touching the boundary by
    // an edge or a vertex is also found
    mfem::Array<int> boundary_markers(m.bdr_attributes.Max());
    boundary_markers = 0;
    boundary_markers[bid - 1] = 1;
    mfem::Array<int> boundary_dofs;
    fes.GetEssentialVDofs(boundary_markers, boundary_dofs);
    auto elts_dofs = std::vector<std::vector<size_type>>{};
    elts_dofs.resize(nc);
    for (size_type i = 0; i < fes.GetNE(); i++) {
      for (size_type c = 0; c != nc; ++c) {
        elts_dofs[c].clear();
      }
      const auto& e = *(fes.GetFE(i));
      const auto nnodes = e.GetDof();
      mfem::Array<int> element_dofs;
      fes.GetElementVDofs(i, element_dofs);
      for (size_type c = 0; c != nc; ++c) {
        auto& edofs = elts_dofs[c];
        for (auto j = 0; j != nnodes; ++j) {
          const auto d = element_dofs[c * nnodes + j];
          if (boundary_dofs[d >= 0 ? d : -1 - d] != 0) {
            edofs.push_back(j);
          }
        }
      }
      for (size_type c = 0; c != nc; ++c) {
        if (!elts_dofs[c].empty()) {
          r.push_back({i, elts_dofs});
          break;
        }
      }
    }
    return r;
  }  // end of getElementsDegreesOfFreedomOnBoundary

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_BOUNDARYUTILITIES_IXX */
