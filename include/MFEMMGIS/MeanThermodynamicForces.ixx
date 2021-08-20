/*!
 * \file   include/MFEMMGIS/MeanThermodynamicForces.ixx
 * \brief
 * \author Thomas Helfer, Hugo Copin
 * \date   08/04/2021
 */

#ifndef LIB_MFEM_MGIS_MEANTHERMODYNAMICFORCES_IXX
#define LIB_MFEM_MGIS_MEANTHERMODYNAMICFORCES_IXX

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  MeanThermodynamicForces<parallel>::MeanThermodynamicForces(
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const Parameters &params) {
    checkParameters(params, {"OutputFileName"});
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0) {
        this->openFile(p, get<std::string>(params, "OutputFileName"));
      }
#else /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->openFile(p, get<std::string>(params, "OutputFileName"));
    }
  }  // end of MeanThermodynamicForces

  template <bool parallel>
  void MeanThermodynamicForces<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const real t,
      const real dt) {
    const auto [tf_integrals, volumes] =
        computeMeanThermodynamicForcesValues(p);
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0) {
        this->out << t + dt;
      }
      for (const auto &mi : p.getAssignedMaterialsIdentifiers()) {
        real v = 0;
        std::vector<double> tf_integral(tf_integrals[mi].size(), 0);
        MPI_Reduce(tf_integrals[mi].data(), tf_integral.data(),
                   tf_integrals[mi].size(), MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        MPI_Reduce(&volumes[mi], &v, 1, MPI_DOUBLE, MPI_SUM, 0,
                   MPI_COMM_WORLD);
        if (rank == 0) {
          this->writeResults(tf_integral, v);
        }
      }
      if (rank == 0) {
        this->out << '\n';
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->out << t + dt;
      for (const auto &mi : p.getAssignedMaterialsIdentifiers()) {
        this->writeResults(tf_integrals[mi], volumes[mi]);
      }
      this->out << '\n';
    }
  }  // end of MeanThermodynamicForces

  template <bool parallel>
  void MeanThermodynamicForces<parallel>::openFile(
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const std::string &f) {
    this->out.open(f);
    if (!this->out) {
      raise("MeanThermodynamicForces::openFile: unable to open file '" + f +
            "'");
    }
    out << "# first column: time\n";
    auto c = mfem_mgis::size_type{2};
    for (const auto &mi : p.getAssignedMaterialsIdentifiers()) {
      const auto &bi = p.getBehaviourIntegrator(mi);
      const auto &s1 = bi.getMaterial().s1;
      const auto thsize =
          static_cast<mfem_mgis::size_type>(s1.thermodynamic_forces_stride);
      for (mfem_mgis::size_type k = 0; k != thsize; ++k, ++c) {
        out << "# " << c << " column: " << k
            << " component of the mean thermodynamic force for material " << mi
            << '\n';
      }
    }
  }

  template <bool parallel>
  void MeanThermodynamicForces<parallel>::writeResults(
      const std::vector<real> &tf_integral, const real v) {
    for (mfem_mgis::size_type k = 0; k != tf_integral.size(); ++k) {
      this->out << " " << tf_integral[k] / v;
    }
  }  // end of writeResults

  template <bool parallel>
  MeanThermodynamicForces<parallel>::~MeanThermodynamicForces() = default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MEANTHERMODYNAMICFORCES_IXX */
