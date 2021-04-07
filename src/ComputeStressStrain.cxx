/*!
 * \file   src/ComputeStressStrain.cxx
 * \brief
 * \author Guillaume Latu
 * \date   07/04/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BoundaryUtilities.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/ComputeStressStrain.hxx"

namespace mfem_mgis {


#ifdef MFEM_USE_MPI

  ComputeStressStrain<true>::ComputeStressStrain(
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& params) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      const auto& f = get<std::string>(params, "Stress");
      auto& fed = p.getFiniteElementDiscretization();
      auto& fes = fed.template getFiniteElementSpace<true>();
    }
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
  }  // end of execute

  ComputeStressStrain<true>::~ComputeStressStrain() =
      default;

#endif /* MFEM_USE_MPI */

  ComputeStressStrain<false>::ComputeStressStrain(
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& params) {
    auto& fed = p.getFiniteElementDiscretization();
    auto& fes = fed.template getFiniteElementSpace<false>();
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
