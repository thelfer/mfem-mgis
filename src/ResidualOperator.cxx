/*!
 * \file   src/ResidualOperator.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/ResidualOperator.hxx"

namespace mfem_mgis {

  ResidualOperator::ResidualOperator(NonLinearEvolutionProblem &p)
      : Operator(p.getFiniteElementSpace().GetTrueVSize()),
        problem(p) {}  // end of ResidualOperator

  void ResidualOperator::Mult(const mfem::Vector &k, mfem::Vector &y) const {
    this->problem.Mult(k, y);
  }  // end of Mult

  mfem::Operator &ResidualOperator::GetGradient(const mfem::Vector &xp) const {
    return this->problem.GetGradient(xp);
  }  // end of GetGradient

  ResidualOperator::~ResidualOperator() = default;

}  // end of namespace mfem_mgis
