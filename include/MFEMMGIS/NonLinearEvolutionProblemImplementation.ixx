/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx
 * \brief  This file implements the `computeResultantForceOnBoundary` and
 * `computeMeanThermodynamicForcesValues` functions
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis::internals {

  /*!
   * \return if two materials have the same thermodynamic forces
   * \param[in] m1: first material
   * \param[in] m2: second material
   */
  inline bool haveSameThermodynamicForces(const Material& m1,
                                          const Material& m2) noexcept {
    const auto& tfs1 = m1.b.thermodynamic_forces;
    const auto& tfs2 = m2.b.thermodynamic_forces;
    if (tfs1.size() != tfs2.size()) {
      return false;
    }
    for (std::size_t i = 0; i != tfs1.size(); ++i) {
      if ((tfs1[i].name != tfs2[i].name) || (tfs1[i].type != tfs2[i].type)) {
        return false;
      }
    }
    return m1.s1.thermodynamic_forces_stride ==
           m2.s1.thermodynamic_forces_stride;
  }  // end of haveSameThermodynamicForces

}  // end of namespace mfem_mgis::internals

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
      // the inner forces of all the behaviour integrators of the material are
      // summed, as in the residual
      const auto onbis = p.getNumberOfBehaviourIntegrators(ctx, tr.Attribute);
      if ((isInvalid(onbis)) || (*onbis == 0)) {
        return ctx.registerErrorMessage("invalid behaviour integrator");
      }
      for (size_type b = 0; b != *onbis; ++b) {
        const auto obi = p.getBehaviourIntegrator(ctx, tr.Attribute, b);
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
      const auto obi = p.getBehaviourIntegrator(ctx, mi, 0);
      if (isInvalid(obi)) {
        return {};
      }
      const auto om = obi->getMaterial(ctx);
      if (isInvalid(om)) {
        return {};
      }
      // the thermodynamic forces of all the behaviour integrators of the
      // material are summed, as in the residual: they must be the same
      const auto onbis = p.getNumberOfBehaviourIntegrators(ctx, mi);
      if (isInvalid(onbis)) {
        return {};
      }
      for (size_type b = 1; b < *onbis; ++b) {
        const auto obi2 = p.getBehaviourIntegrator(ctx, mi, b);
        if (isInvalid(obi2)) {
          return {};
        }
        const auto om2 = obi2->getMaterial(ctx);
        if (isInvalid(om2)) {
          return {};
        }
        if (!internals::haveSameThermodynamicForces(*om, *om2)) {
          return ctx.registerErrorMessage(
              "computeMeanThermodynamicForcesValues: the behaviour "
              "integrators of material '" +
              std::to_string(mi) +
              "' do not have the same thermodynamic forces");
        }
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
      const auto onbis = p.getNumberOfBehaviourIntegrators(ctx, tr.Attribute);
      if ((isInvalid(onbis)) || (*onbis == 0)) {
        return ctx.registerErrorMessage("invalid behaviour integrator");
      }
      auto& s = stress_integrals[tr.Attribute];
      auto& v = volumes[tr.Attribute];
      for (size_type b = 0; b != *onbis; ++b) {
        const auto obi = p.getBehaviourIntegrator(ctx, tr.Attribute, b);
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
          // the volume is only accumulated with the first behaviour integrator
          if (b == 0) {
            v += w;
          }
        }
      }
    }
    return std::pair<std::vector<std::vector<real>>, std::vector<real>>{
        std::move(stress_integrals), std::move(volumes)};
  }  // end of computeMeanThermodynamicForcesValues

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_IXX */
