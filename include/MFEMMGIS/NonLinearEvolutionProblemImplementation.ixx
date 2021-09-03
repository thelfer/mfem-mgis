/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

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

  template <bool parallel>
  std::pair<std::vector<std::vector<real>>, std::vector<real>>
  computeMeanThermodynamicForcesValues(
      NonLinearEvolutionProblemImplementation<parallel>& p) {
    auto nmax = p.getFiniteElementSpace().GetMesh()->attributes.Max() + 1;
    std::vector<std::vector<mfem_mgis::real>> stress_integrals(nmax);
    const auto& mis = p.getAssignedMaterialsIdentifiers();
    for (const auto& mi : mis) {
      const auto& bi = p.getBehaviourIntegrator(mi);
      const auto& s1 = bi.getMaterial().s1;
      const auto thsize = s1.thermodynamic_forces_stride;
      stress_integrals[mi].resize(thsize, mfem_mgis::real(0));
    }
    std::vector<mfem_mgis::real> volumes(stress_integrals.size(),
                                         mfem_mgis::real(0));
    const auto& fes = p.getFiniteElementSpace();
    for (mfem_mgis::size_type i = 0; i < fes.GetNE(); i++) {
      auto& e = *(fes.GetFE(i));
      auto& tr = *(fes.GetElementTransformation(i));
      const auto& bi = p.getBehaviourIntegrator(tr.Attribute);
      const auto& ir = bi.getIntegrationRule(e, tr);
      const auto& m = bi.getMaterial();
      const auto& s1 = bi.getMaterial().s1;
      const auto& qspace = m.getPartialQuadratureSpace();
      const auto thsize =
          static_cast<mfem_mgis::size_type>(s1.thermodynamic_forces_stride);
      auto& s = stress_integrals[tr.Attribute];
      auto& v = volumes[tr.Attribute];
      const auto eoffset = qspace.getOffset(tr.ElementNo);
      for (mfem_mgis::size_type j = 0; j < ir.GetNPoints(); j++) {
        const auto o = eoffset + j;
        const auto& ip = ir.IntPoint(j);
        tr.SetIntPoint(&ip);
        const auto thf = s1.thermodynamic_forces.subspan(o * thsize, thsize);
        const auto w = bi.getIntegrationPointWeight(tr, ip);
        for (mfem_mgis::size_type k = 0; k != thsize; ++k) {
          s[k] += w * thf[k];
        }
        v += w;
      }
    }
    return {stress_integrals, volumes};
  }  // end of computeMeanThermodynamicForcesValues

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX */
