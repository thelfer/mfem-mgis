/*!
 * \file   UniformDirichletBoundaryCondition.cxx
 * \brief
 * \author Thomas Helfer
 * \date   18/03/2021
 */

#include <iostream>

#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"

namespace mfem_mgis {

  UniformDirichletBoundaryCondition::UniformDirichletBoundaryCondition(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const size_type bid,
      const size_type c)
      : DirichletBoundaryConditionBase(fed, bid, c),
        ufct([](const real) { return real(0); }) {
  }  // end of UniformDirichletBoundaryCondition

  UniformDirichletBoundaryCondition::UniformDirichletBoundaryCondition(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const size_type bid,
      const size_type c,
      std::function<real(const real)> uvalues)
      : DirichletBoundaryConditionBase(fed, bid, c),
        ufct(uvalues) {}  // end of UniformDirichletBoundaryCondition

  void UniformDirichletBoundaryCondition::updateImposedValues(
      mfem::Vector& u, const real t) const {
    const auto uv = this->ufct(t);
    for (const auto& dof : this->dofs) {
      u[dof] = uv;
    }
  }  // end of updateImposedValues

  void UniformDirichletBoundaryCondition::setImposedValuesIncrements(
      mfem::Vector& u, const real ti, const real te) const {
    std::cerr << "ti: " << ti << '\n';
    std::cerr << "te: " << te << '\n';
    std::cerr << "u(ti): " << this->ufct(ti) << '\n';
    std::cerr << "u(te): " << this->ufct(te) << '\n';
    const auto duv = this->ufct(te) - this->ufct(ti);
    std::cerr << "duv: " << duv << std::endl;
    for (const auto& dof : this->dofs) {
      u[dof] = duv;
    }
  }  // end of setImposedValuesIncrements

  UniformDirichletBoundaryCondition::~UniformDirichletBoundaryCondition() =
      default;

}  // end of namespace mfem_mgis
