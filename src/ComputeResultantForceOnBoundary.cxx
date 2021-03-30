/*!
 * \file   src/ComputeResultantForceOnBoundary.cxx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/BoundaryUtilities.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/ComputeResultantForceOnBoundary.hxx"

namespace mfem_mgis {

  static void writeOuputFileHeader(std::ofstream& out,
                                   const size_type bid,
                                   const size_type nc) {
    out << "# first column: time\n";
    for (size_type i = 0; i != nc; ++i) {
      out << "# " << i + 1 << "th column: " << i + 1
          << " component of the resultant of the inner forces on boundary '"
          << bid << "'\n";
    }
  }  // end of writeOuputFileHeader

  static void writeResultantForce(std::ofstream& out,
                                  const mfem::Vector& F,
                                  const real t) {
    out << t;
    for (size_type i = 0; i != F.Size(); ++i) {
      out << " " << F[i];
    }
    out << "\n";
  }  // end of writeResultantForce

  ComputeResultantForceOnBoundaryCommon::ComputeResultantForceOnBoundaryCommon(
      std::vector<std::pair<size_type, size_type>> fd, const size_type i)
      : faces(std::move(fd)),
        bid(i) {}  // end of ComputeResultantForceOnBoundaryCommon

#ifdef MFEM_USE_MPI

  ComputeResultantForceOnBoundary<true>::ComputeResultantForceOnBoundary(
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& params)
      : ComputeResultantForceOnBoundaryCommon(
            buildFacesDescription<true>(p, get<int>(params, "Boundary")),
            get<int>(params, "Boundary")) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      const auto& f = get<std::string>(params, "OutputFileName");
      this->out.open(f);
      if (!this->out) {
        mgis::raise("can't open file '" + f + "'");
      }
      auto& fed = p.getFiniteElementDiscretization();
      auto& fes = fed.template getFiniteElementSpace<true>();
      writeOuputFileHeader(this->out, this->bid, fes.GetVDim());
    }
  }  // end of ComputeResultantForceOnBoundary

  void ComputeResultantForceOnBoundary<true>::execute(
      NonLinearEvolutionProblemImplementation<true>& p,
      const real t,
      const real dt) {
    mfem::Vector F;
    computeResultantForceOnBoundary(F, p, this->faces);
    //
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      mfem::Vector gF;
      gf.SetSize(F.GetSize());
      MPI_Reduce(F.data(), gF.data(), F.GetSize(), MPI_DOUBLE, MPI_SUM, 0,
                 MPI_COMM_WORLD);
      writeResultantForce(this->out, gF, t + dt);
    }
  }  // end of execute

  ComputeResultantForceOnBoundary<true>::~ComputeResultantForceOnBoundary() =
      default;

#endif /* MFEM_USE_MPI */

  ComputeResultantForceOnBoundary<false>::ComputeResultantForceOnBoundary(
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& params)
      : ComputeResultantForceOnBoundaryCommon(
            buildFacesDescription<false>(p, get<int>(params, "Boundary")),
            get<int>(params, "Boundary")) {
    const auto& f = get<std::string>(params, "OutputFileName");
    this->out.open(f);
    if (!this->out) {
      mgis::raise("can't open file '" + f + "'");
    }
    auto& fed = p.getFiniteElementDiscretization();
    auto& fes = fed.template getFiniteElementSpace<false>();
    writeOuputFileHeader(this->out, this->bid, fes.GetVDim());
  }  // end of ComputeResultantForceOnBoundary

  void ComputeResultantForceOnBoundary<false>::execute(
      NonLinearEvolutionProblemImplementation<false>& p,
      const real t,
      const real dt) {
    mfem::Vector F;
    computeResultantForceOnBoundary(F, p, this->faces);
    writeResultantForce(this->out, F, t + dt);
  }  // end of execute

  ComputeResultantForceOnBoundary<false>::~ComputeResultantForceOnBoundary() =
      default;

}  // end of namespace mfem_mgis