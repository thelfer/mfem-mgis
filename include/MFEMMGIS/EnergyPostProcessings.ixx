/*!
 * \file   include/MFEMMGIS/EnergyPostProcessings.ixx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2021
 */

#ifndef LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_IXX
#define LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_IXX

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  EnergyPostProcessingBase<parallel>::EnergyPostProcessingBase(
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const Parameters &params,
      const std::string_view etype)
      : materials_identifiers(getMaterialsIdentifiers(p, params)) {
    checkParameters(params, {"OutputFileName", "Material", "Materials"});
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0) {
        this->openFile(get<std::string>(params, "OutputFileName"), etype);
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->openFile(get<std::string>(params, "OutputFileName"), etype);
    }
  }  // end of EnergyPostProcessingBase

  template <bool parallel>
  void EnergyPostProcessingBase<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const real t,
      const real dt) {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0) {
        this->out << t + dt;
      }
      const auto energies = this->computeEnergies(p);
      if (rank == 0) {
        this->writeResults(energies);
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->out << t + dt;
      const auto energies = this->computeEnergies(p);
      this->writeResults(energies);
    }
  }  // end of EnergyPostProcessingBase

  template <bool parallel>
  void EnergyPostProcessingBase<parallel>::openFile(
      const std::string &f, const std::string_view etype) {
    this->out.open(f);
    if (!this->out) {
      raise("EnergyPostProcessingBase::openFile: unable to open file '" + f +
            "'");
    }
    this->out << "# first column: time\n";
    auto c = size_type{2};
    for (const auto &m : this->materials_identifiers) {
      this->out << "# " << c << "column: " << etype  //
                << "energy of material (" << m << ")\n";
      ++c;
    }
  }  // end of openFile

  template <bool parallel>
  void EnergyPostProcessingBase<parallel>::writeResults(
      const std::vector<real> &energies) {
    for (const auto &e : energies) {
      this->out << " " << e;
    }
    this->out << std::endl;
  }  // end of writeResults

  template <bool parallel>
  EnergyPostProcessingBase<parallel>::~EnergyPostProcessingBase() = default;

  template <bool parallel>
  StoredEnergyPostProcessing<parallel>::StoredEnergyPostProcessing(
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const Parameters &params)
      : EnergyPostProcessingBase<parallel>(p, params, "stored") {
  }  // end of StoredEnergyPostProcessing

  template <bool parallel>
  std::vector<real> StoredEnergyPostProcessing<parallel>::computeEnergies(
      const AbstractNonLinearEvolutionProblem &p) const {
    auto energies = std::vector<real>{};
    energies.reserve(this->materials_identifiers.size());
    for (const auto &m : this->materials_identifiers) {
      energies.push_back(computeStoredEnergy(p.getBehaviourIntegrator(m)));
    }
    return energies;
  }  // end of computeEnergies

  template <bool parallel>
  StoredEnergyPostProcessing<parallel>::~StoredEnergyPostProcessing() = default;

  template <bool parallel>
  DissipatedEnergyPostProcessing<parallel>::DissipatedEnergyPostProcessing(
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const Parameters &params)
      : EnergyPostProcessingBase<parallel>(p, params, "dissipated") {
  }  // end of DissipatedEnergyPostProcessing

  template <bool parallel>
  std::vector<real> DissipatedEnergyPostProcessing<parallel>::computeEnergies(
      const AbstractNonLinearEvolutionProblem &p) const {
    auto energies = std::vector<real>{};
    energies.reserve(this->materials_identifiers.size());
    for (const auto &m : this->materials_identifiers) {
      energies.push_back(computeDissipatedEnergy(p.getBehaviourIntegrator(m)));
    }
    return energies;
  }  // end of computeEnergies

  template <bool parallel>
  DissipatedEnergyPostProcessing<parallel>::~DissipatedEnergyPostProcessing() =
      default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_IXX */
