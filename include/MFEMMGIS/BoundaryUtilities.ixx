/*!
 * \file   include/MFEMMGIS/BoundaryUtilities.ixx
 * \brief    
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_BOUNDARYUTILITIES_IXX
#define LIB_MFEM_MGIS_BOUNDARYUTILITIES_IXX

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static std::vector<std::pair<size_type, size_type>> buildFacesDescription(
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
  }  // end of buildFaces

}  // end of namespace mfem_mgis


#endif /* LIB_MFEM_MGIS_BOUNDARYUTILITIES_IXX */
