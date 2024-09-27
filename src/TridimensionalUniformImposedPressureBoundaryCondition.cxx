/*!
 * \file   src/TridimensionalUniformImposedPressureBoundaryCondition.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/09/2024
 */

#include "MFEMMGIS/TridimensionalUniformImposedPressureBoundaryCondition.hxx"

namespace mfem_mgis {

  TridimensionalUniformImposedPressureBoundaryCondition::
      TridimensionalUniformImposedPressureBoundaryCondition(
          std::function<real(const real)> prvalues)
      : prfct(prvalues) {
  }  // end of TridimensionalUniformImposedPressureBoundaryCondition

  void
  TridimensionalUniformImposedPressureBoundaryCondition::AssembleElementVector(
      const mfem::FiniteElement &e,
      mfem::ElementTransformation &tr,
      const mfem::Vector &,
      mfem::Vector &R) {
    const int dim = e.GetDim() + 1;
    const int dof = e.GetDof();
    mfem::Vector n(dim);
    if (dim == 1) {
      n[0] = 1.0;
    }
#ifdef MFEM_THREAD_SAFE
    mfem::Vector shape;
#endif
    //
    shape.SetSize(dof);
    // initialize the residual
    R.SetSize(dof * dim);
    R = 0.0;
    // selection of the integration rule
    const auto *ir = this->IntRule;
    if (ir == nullptr) {
      const int o = e.GetOrder();
      ir = &mfem::IntRules.Get(e.GetGeomType(), o);
    }
    // loop over the integration point
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const auto &ip = ir->IntPoint(i);
      tr.SetIntPoint(&ip);
      // computation of the normal
      if (dim > 1) {
        CalcOrtho(tr.Jacobian(), n);
      }
      // evaluation of the pressure
      const auto pr = this->prfct(0);
      // computation of the external forces
      e.CalcShape(ip, shape);
      //      R.Add(ip.weight * (pr * n), shape);
    }
  }  // end of AssembleElementVector

  void
  TridimensionalUniformImposedPressureBoundaryCondition::AssembleElementGrad(
      const mfem::FiniteElement &,
      mfem::ElementTransformation &,
      const mfem::Vector &,
      mfem::DenseMatrix &) {}

  TridimensionalUniformImposedPressureBoundaryCondition::
      ~TridimensionalUniformImposedPressureBoundaryCondition() = default;

}  // end of namespace mfem_mgis