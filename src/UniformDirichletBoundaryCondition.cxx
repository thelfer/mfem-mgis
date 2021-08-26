/*!
 * \file   UniformDirichletBoundaryCondition.cxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Parameter.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"

namespace mfem_mgis {

  UniformDirichletBoundaryCondition::UniformDirichletBoundaryCondition(
      AbstractNonLinearEvolutionProblem& p, const Parameters& params)
      : DirichletBoundaryConditionBase(
            p.getFiniteElementDiscretization(),
            getBoundariesIdentifiers(p, params, false),
            get<size_type>(params, "Component")) {
    this->ufct = get_if<std::function<real(const real)>>(
        params, ("LoadingEvolution"), [](const real) { return real(0); });
  }  // end of UniformDirichletBoundaryCondition

  UniformDirichletBoundaryCondition::UniformDirichletBoundaryCondition(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const size_type bid,
      const size_type c)
      : DirichletBoundaryConditionBase(*fed, {bid}, c),
        ufct([](const real) { return real(0); }) {
  }  // end of UniformDirichletBoundaryCondition

  UniformDirichletBoundaryCondition::UniformDirichletBoundaryCondition(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const size_type bid,
      const size_type c,
      std::function<real(const real)> uvalues)
      : DirichletBoundaryConditionBase(*fed, {bid}, c),
        ufct(uvalues) {}  // end of UniformDirichletBoundaryCondition

  void UniformDirichletBoundaryCondition::updateImposedValues(
      mfem::Vector& u, const real t) const {
    const auto uv = this->ufct(t);
    for (const auto& dof : this->dofs) {
      u[dof] = uv;
    }
  }  // end of updateImposedValues

  void UniformDirichletBoundaryCondition::setImposedValuesIncrements(
      mfem::Vector& du, const real ti, const real te) const {
    const auto duv = this->ufct(te) - this->ufct(ti);
    for (const auto& dof : this->dofs) {
      du[dof] = duv;
    }
  }  // end of setImposedValuesIncrements

  UniformDirichletBoundaryCondition::~UniformDirichletBoundaryCondition() =
      default;

}  // end of namespace mfem_mgis
