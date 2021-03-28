/*!
 * \file   src/DirichletBoundaryConditionBase.cxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#include "mfem/general/array.hpp"
#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif MFEM_USE_MPI
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/DirichletBoundaryConditionBase.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static std::vector<size_type> buildDegreesOfFredomList(
      FiniteElementDiscretization& fed,
      const size_type bid,
      const size_type c) {
    auto& m = fed.getMesh<parallel>();
    auto boundaries_markers =
        mfem::Array<mfem_mgis::size_type>(m.bdr_attributes.Max());
    boundaries_markers = 0;
    boundaries_markers[bid - 1] = 1;
    auto tmp = mfem::Array<mfem_mgis::size_type>{};
    fed.getFiniteElementSpace<parallel>().GetEssentialTrueDofs(
        boundaries_markers, tmp, c);
    return std::vector<size_type>(tmp.GetData(), tmp.GetData() + tmp.Size());
  }  // end of buildDegreesOfFredomList

  static std::vector<size_type> buildDegreesOfFredomListDispatch(
      FiniteElementDiscretization& fed,
      const size_type bid,
      const size_type c) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return buildDegreesOfFredomList<true>(fed, bid, c);
#else
      mgis::raise(
          "DirichletBoundaryConditionBase::DirichletBoundaryConditionBase: "
          "unsupported parallel computations");
#endif
    }
    return buildDegreesOfFredomList<false>(fed, bid, c);
  }  // end of buildDegreesOfFredomListDispatch

  DirichletBoundaryConditionBase::DirichletBoundaryConditionBase(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const size_type bid,
      const size_type c)
      : dofs(buildDegreesOfFredomListDispatch(*fed, bid, c)) {
  }  // end of DirichletBoundaryConditionBase::DirichletBoundaryConditionBase

  std::vector<size_type>
  DirichletBoundaryConditionBase::getHandledDegreesOfFreedom() const {
    return this->dofs;
  } // end of getHandledDegreesOfFreedom

  DirichletBoundaryConditionBase::~DirichletBoundaryConditionBase() = default;

}  // end of namespace mfem_mgis
