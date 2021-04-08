/*!
 * \file   src/ComputeStressStrain.cxx
 * \brief
 * \author Guillaume Latu
 * \date   07/04/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/BoundaryUtilities.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ComputeStressStrain.hxx"

namespace mfem_mgis {

  template <bool parallel>
  class DiagCoefficient : public mfem::Coefficient
  {
  protected:
    const size_type _icomp; // component to evaluate, 0 <= _icom < max
    mfem::DenseMatrix _grad; // auxiliary matrix, used in Eval
    const Material _m; // problem
    
  public:
    DiagCoefficient(NonLinearEvolutionProblemImplementation<parallel> &p, size_type icomp = 0)
      : _icomp(icomp), _m(p.getMaterial()) { }
    
    void SetComponent(size_type icomp) { _icomp = icomp; }
    
    virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) = 0;

    //    DiagCoefficient::~DiagCoefficient() = default;
  };
  
  ComputeStressStrainCommon::ComputeStressStrainCommon(
     size_type picomp, std::shared_ptr<FiniteElementDiscretization> pfed)
    : icomp(picomp), fed(pfed) {}  // end of ComputeStressStrainCommon

#ifdef MFEM_USE_MPI
  ComputeStressStrain<true>::ComputeStressStrain(
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& params)
    : ComputeStressStrainCommon(get<int>(params, "StressComponent"),p.getFiniteElementDiscretizationPointer()) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      auto& fes = fed->getFiniteElementSpace<true>();
    }
    std::cout << "WIP" << std::endl;
    mfem_mgis::abort(3);
  }  // end of ComputeStressStrain

  void ComputeStressStrain<true>::execute(
      NonLinearEvolutionProblemImplementation<true>& p,
      const real t,
      const real dt) {
    mfem::Vector F;
    // Do the computation   computeXXX(F, p);
    //
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    mfem::Vector gF;
    gF.SetSize(F.Size());
    MPI_Reduce(F.GetData(), gF.GetData(), F.Size(), MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
    std::cout << "WIP" << std::endl;
    mfem_mgis::abort(3);
  }  // end of execute

  ComputeStressStrain<true>::~ComputeStressStrain() =
      default;

#endif /* MFEM_USE_MPI */

  ComputeStressStrain<false>::ComputeStressStrain(
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& params)
    : ComputeStressStrainCommon(get<int>(params, "StressComponent"),p.getFiniteElementDiscretizationPointer()) {
    auto& fes = fed->getFiniteElementSpace<false>();
    //    fieldspace = new FiniteElementSpace(mesh, fec, 1, Ordering::byVDIM);
  }  // end of ComputeStressStrain

  void ComputeStressStrain<false>::execute(
      NonLinearEvolutionProblemImplementation<false>& p,
      const real t,
      const real dt) {
    mfem::Vector F;
    // Do the computation   computeXXX(F, p);
    //
  }  // end of execute

  ComputeStressStrain<false>::~ComputeStressStrain() =
      default;

}  // end of namespace mfem_mgis
