/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX

#include "MFEMMGIS/BehaviourIntegrator.hxx"

namespace mfem_mgis {

  template <bool parallel>
  void computeResultantForceOnBoundary(
      mfem::Vector& F,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const std::vector<std::pair<size_type, size_type>>& faces) {
    auto& fed = p.getFiniteElementDiscretization();
    auto& fes = fed.template getFiniteElementSpace<parallel>();
    const auto nc = fes.GetVDim();
    mfem::Array<int> cell_dofs;
    mfem::Array<int> face_dofs;
    mfem::Vector cell_forces;
    mfem::Vector face_forces;
    F.SetSize(nc);
    F = real{0};
    for (const auto& f : faces) {
      const auto& fe = *(fes.GetFE(std::get<1>(f)));
      auto& tr = *(fes.GetElementTransformation(std::get<1>(f)));
      fes.GetElementVDofs(std::get<1>(f), cell_dofs);
      // compute the inner forces
      auto& bi = p.getBehaviourIntegrator(tr.Attribute);
      bi.computeInnerForces(cell_forces, fe, tr);
      // computation of the resultant by summing the contributions on the
      // considered boundary
      fes.GetBdrElementVDofs(std::get<0>(f), face_dofs);
      face_forces.SetSize(face_dofs.Size());
      const auto n_cell_nodes = fe.GetDof();
      const auto n_face_nodes = face_dofs.Size() / fe.GetDim();
      const auto cell_dofs_begin = cell_dofs.GetData();
      // here we take into account the fact that that the components are stored
      // contiguously
      for (size_type c = 0; c != nc; ++c) {
        const auto cell_dofs_c_begin = cell_dofs_begin + c * n_cell_nodes;
        const auto cell_dofs_c_end = cell_dofs_c_begin + n_cell_nodes;
        for (size_type i = 0; i != n_face_nodes; ++i) {
          const auto fid = i + c * n_face_nodes;
          const auto face_dof = face_dofs[fid];
          const auto pc =
              std::find(cell_dofs_c_begin, cell_dofs_c_end, face_dof);
          F[c] += cell_forces[pc - cell_dofs_begin];
        }
      }
    }
  }  // end of computeResultantForceOnBoundary

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX */
