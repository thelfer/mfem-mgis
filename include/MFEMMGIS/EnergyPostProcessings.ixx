/*!
 * \file   include/MFEMMGIS/EnergyPostProcessings.ixx
 * \brief  This file implements the inline functions declared in
 * `MFEMMGIS/EnergyPostProcessings.hxx`
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
      : materials_identifiers(getMaterialsIdentifiers(throwing, p, params)),
        behaviour_integrators(
            getBehaviourIntegratorsSelection(throwing, params)) {
    checkParameters(
        throwing, params,
        {"OutputFileName", "Material", "Materials", "BehaviourIntegrator"});
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      int rank;
      MPI_Comm_rank(getMPICommunicator(p), &rank);
      if (rank == 0) {
        this->openFile(get<std::string>(throwing, params, "OutputFileName"),
                       etype);
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->openFile(get<std::string>(throwing, params, "OutputFileName"),
                     etype);
    }
  }  // end of EnergyPostProcessingBase

  template <bool parallel>
  bool EnergyPostProcessingBase<parallel>::execute(
      Context &ctx,
      NonLinearEvolutionProblemImplementation<parallel> &p,
      const real t,
      const real dt) noexcept {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      int rank;
      MPI_Comm_rank(getMPICommunicator(p), &rank);
      const auto oenergies = this->computeEnergies(ctx, p);
      if (isInvalid(oenergies)) {
        return false;
      }
      if (rank == 0) {
        this->out << t + dt;
        this->writeResults(*oenergies);
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      const auto oenergies = this->computeEnergies(ctx, p);
      if (isInvalid(oenergies)) {
        return false;
      }
      this->out << t + dt;
      this->writeResults(*oenergies);
    }
    return true;
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
      this->out << "# column " << c << ": " << etype  //
                << " energy of material (" << m << ")\n";
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
  std::optional<std::vector<real>>
  StoredEnergyPostProcessing<parallel>::computeEnergies(
      Context &ctx, const AbstractNonLinearEvolutionProblem &p) const noexcept {
    auto energies = std::vector<real>{};
    energies.reserve(this->materials_identifiers.size());
    for (const auto &m : this->materials_identifiers) {
      // sum of the energies of the selected behaviour integrators
      const auto obis = getSelectedBehaviourIntegrators(
          ctx, p, m, this->behaviour_integrators);
      if (isInvalid(obis)) {
        return {};
      }
      auto energy = real{};
      for (auto b = obis->first; b != obis->second; ++b) {
        const auto obi = p.getBehaviourIntegrator(ctx, m, b);
        if (isInvalid(obi)) {
          return {};
        }
        const auto oe = computeStoredEnergy(ctx, *obi);
        if (isInvalid(oe)) {
          return {};
        }
        energy += *oe;
      }
      energies.push_back(energy);
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
  std::optional<std::vector<real>>
  DissipatedEnergyPostProcessing<parallel>::computeEnergies(
      Context &ctx, const AbstractNonLinearEvolutionProblem &p) const noexcept {
    auto energies = std::vector<real>{};
    energies.reserve(this->materials_identifiers.size());
    for (const auto &m : this->materials_identifiers) {
      // sum of the energies of the selected behaviour integrators
      const auto obis = getSelectedBehaviourIntegrators(
          ctx, p, m, this->behaviour_integrators);
      if (isInvalid(obis)) {
        return {};
      }
      auto energy = real{};
      for (auto b = obis->first; b != obis->second; ++b) {
        const auto obi = p.getBehaviourIntegrator(ctx, m, b);
        if (isInvalid(obi)) {
          return {};
        }
        const auto oe = computeDissipatedEnergy(ctx, *obi);
        if (isInvalid(oe)) {
          return {};
        }
        energy += *oe;
      }
      energies.push_back(energy);
    }
    return energies;
  }  // end of computeEnergies

  template <bool parallel>
  DissipatedEnergyPostProcessing<parallel>::~DissipatedEnergyPostProcessing() =
      default;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_IXX */
