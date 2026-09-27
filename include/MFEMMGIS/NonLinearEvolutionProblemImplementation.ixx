/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  template <bool parallel>
  bool computeResultantForceOnBoundary(
      Context& ctx,
      mfem::Vector& F,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const std::vector<
          std::pair<size_type, std::vector<std::vector<size_type>>>>&
          elements) noexcept {
    auto& fed = p.getFiniteElementDiscretization();
    auto& fes = fed.template getFiniteElementSpace<parallel>();
    const auto nc = fes.GetVDim();
    mfem::Vector elt_forces;
    F.SetSize(nc);
    F = real{0};
    for (const auto& e : elements) {
      const auto& fe = *(fes.GetFE(e.first));
      auto& tr = *(fes.GetElementTransformation(e.first));
      const auto nnodes = fe.GetDof();
      // compute the inner forces
#pragma message("FIXME: invalid if multiple behaviour integrators is defined")
      const auto obi = p.getBehaviourIntegrator(ctx, tr.Attribute, 0);
      if (isInvalid(obi)) {
        return ctx.registerErrorMessage("invalid behaviour integrator");
      }
      obi->computeInnerForces(elt_forces, fe, tr);
      for (size_type c = 0; c != nc; ++c) {
        const auto* const Fe = elt_forces.GetData() + c * nnodes;
        for (const auto& i : e.second[c]) {
          F[c] += Fe[i];
        }
      }
    }
    return true;
  }  // end of computeResultantForceOnBoundary

  template <bool parallel>
  std::optional<std::pair<std::vector<std::vector<real>>, std::vector<real>>>
  computeMeanThermodynamicForcesValues(
      Context& ctx,
      NonLinearEvolutionProblemImplementation<parallel>& p) noexcept {
    auto nmax = p.getFiniteElementSpace().GetMesh()->attributes.Max() + 1;
    std::vector<std::vector<mfem_mgis::real>> stress_integrals(nmax);
    const auto& mis = p.getAssignedMaterialsIdentifiers();
    for (const auto& mi : mis) {
#pragma message("FIXME: invalid if multiple behaviour integrators is defined")
      const auto obi = p.getBehaviourIntegrator(ctx, mi, 0);
      if (isInvalid(obi)) {
        return {};
      }
      const auto om = obi->getMaterial(ctx);
      if (isInvalid(om)) {
        return {};
      }
      const auto& s1 = om->s1;
      const auto thsize = s1.thermodynamic_forces_stride;
      stress_integrals[mi].resize(thsize, mfem_mgis::real(0));
    }
    std::vector<mfem_mgis::real> volumes(stress_integrals.size(),
                                         mfem_mgis::real(0));
    const auto& fes = p.getFiniteElementSpace();
    for (mfem_mgis::size_type i = 0; i < fes.GetNE(); i++) {
      auto& e = *(fes.GetFE(i));
      auto& tr = *(fes.GetElementTransformation(i));
#pragma message("FIXME: invalid if multiple behaviour integrators is defined")
      const auto obi = p.getBehaviourIntegrator(ctx, tr.Attribute, 0);
      if (isInvalid(obi)) {
        return {};
      }
      const auto& ir = obi->getIntegrationRule(e, tr);
      const auto& om = obi->getMaterial(ctx);
      if (isInvalid(om)) {
        return {};
      }
      const auto& s1 = om->s1;
      const auto& qspace = om->getPartialQuadratureSpace();
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
        const auto w = obi->getIntegrationPointWeight(tr, ip);
        if (om->b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
          const auto r = om->getRotationMatrixAtIntegrationPoint(o);
          std::vector<real> rthf(thf.begin(), thf.end());
          om->b.rotate_thermodynamic_forces_ptr(rthf.data(), rthf.data(),
                                                r.data());
          for (mfem_mgis::size_type k = 0; k != thsize; ++k) {
            s[k] += w * rthf[k];
          }
        } else {
          for (mfem_mgis::size_type k = 0; k != thsize; ++k) {
            s[k] += w * thf[k];
          }
        }
        v += w;
      }
    }
    return std::pair<std::vector<std::vector<real>>, std::vector<real>>{
        std::move(stress_integrals), std::move(volumes)};
  }  // end of computeMeanThermodynamicForcesValues

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX */
