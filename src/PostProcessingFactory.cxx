/*!
 * \file   src/PostProcessingFactory.cxx
 * \brief  This file implements the `PostProcessingFactory` class
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"
#include "MFEMMGIS/ParaviewExportResults.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/MeanThermodynamicForces.hxx"
#include "MFEMMGIS/EnergyPostProcessings.hxx"
#include "MFEMMGIS/ComputeResultantForceOnBoundary.hxx"
#include "MFEMMGIS/PostProcessingFactory.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  PostProcessingFactory<true>& PostProcessingFactory<true>::getFactory() {
    static PostProcessingFactory<true> factory;
    return factory;
  }  // end of getFactory

  void PostProcessingFactory<true>::add(std::string_view n, Generator g) {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      std::string msg("PostProcessingFactory<true>::add: ");
      msg += "a post-processing called '";
      msg += n;
      msg += "' has already been declared";
      raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<true>>
  PostProcessingFactory<true>::generate(
      Context& ctx,
      std::string_view n,
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& params) const noexcept {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      return ctx.registerErrorMessage(
          "PostProcessingFactory<true>::generate: no post-processing called '" +
          std::string{n} + "' declared");
    }
    return pg->second(ctx, p, params);
  }  // end of generate

  PostProcessingFactory<true>::PostProcessingFactory() {
    this->add("ParaviewExportResults",
              [](Context& ctx, NonLinearEvolutionProblemImplementation<true>& p,
                 const Parameters& params) noexcept {
                return make_unique<ParaviewExportResults<true>>(ctx, p, params);
              });
    this->add(
        "ParaviewExportIntegrationPointResultsAtNodes",
        [](Context& ctx, NonLinearEvolutionProblemImplementation<true>& p,
           const Parameters& params) noexcept {
          return make_unique<
              ParaviewExportIntegrationPointResultsAtNodesImplementation<true>>(
              ctx, p, params);
        });
    this->add("ComputeResultantForceOnBoundary",
              [](Context& ctx, NonLinearEvolutionProblemImplementation<true>& p,
                 const Parameters& params) noexcept {
                return make_unique<ComputeResultantForceOnBoundary<true>>(
                    ctx, p, params);
              });
    this->add("MeanThermodynamicForces",
              [](Context& ctx, NonLinearEvolutionProblemImplementation<true>& p,
                 const Parameters& params) noexcept {
                return make_unique<MeanThermodynamicForces<true>>(ctx, p,
                                                                  params);
              });
    this->add("StoredEnergy",
              [](Context& ctx, NonLinearEvolutionProblemImplementation<true>& p,
                 const Parameters& params) noexcept {
                return make_unique<StoredEnergyPostProcessing<true>>(ctx, p,
                                                                     params);
              });
    this->add("DissipatedEnergy",
              [](Context& ctx, NonLinearEvolutionProblemImplementation<true>& p,
                 const Parameters& params) noexcept {
                return make_unique<DissipatedEnergyPostProcessing<true>>(
                    ctx, p, params);
              });
  }  // end of PostProcessingFactory

  PostProcessingFactory<true>::~PostProcessingFactory() = default;

#endif /* MFEM_USE_MPI */

  PostProcessingFactory<false>& PostProcessingFactory<false>::getFactory() {
    static PostProcessingFactory<false> factory;
    return factory;
  }  // end of getFactory

  void PostProcessingFactory<false>::add(std::string_view n, Generator g) {
    const auto pg = this->generators.find(n);
    if (pg != this->generators.end()) {
      std::string msg("PostProcessingFactory<false>::add: ");
      msg += "a post-processing called '";
      msg += n;
      msg += "' has already been declared";
      raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<false>>
  PostProcessingFactory<false>::generate(
      Context& ctx,
      std::string_view n,
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& params) const noexcept {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      return ctx.registerErrorMessage(
          "PostProcessingFactory<false>::generate: no post-processing called "
          "'" +
          std::string{n} + "' declared");
    }
    return pg->second(ctx, p, params);
  }  // end of generate

  PostProcessingFactory<false>::PostProcessingFactory() {
    this->add(
        "ParaviewExportResults",
        [](Context& ctx, NonLinearEvolutionProblemImplementation<false>& p,
           const Parameters& params) noexcept {
          return make_unique<ParaviewExportResults<false>>(ctx, p, params);
        });
    this->add(
        "ParaviewExportIntegrationPointResultsAtNodes",
        [](Context& ctx, NonLinearEvolutionProblemImplementation<false>& p,
           const Parameters& params) noexcept {
          return make_unique<
              ParaviewExportIntegrationPointResultsAtNodesImplementation<
                  false>>(ctx, p, params);
        });
    this->add(
        "ComputeResultantForceOnBoundary",
        [](Context& ctx, NonLinearEvolutionProblemImplementation<false>& p,
           const Parameters& params) noexcept {
          return make_unique<ComputeResultantForceOnBoundary<false>>(ctx, p,
                                                                     params);
        });
    this->add(
        "MeanThermodynamicForces",
        [](Context& ctx, NonLinearEvolutionProblemImplementation<false>& p,
           const Parameters& params) noexcept {
          return make_unique<MeanThermodynamicForces<false>>(ctx, p, params);
        });
    this->add(
        "StoredEnergy",
        [](Context& ctx, NonLinearEvolutionProblemImplementation<false>& p,
           const Parameters& params) noexcept {
          return make_unique<StoredEnergyPostProcessing<false>>(ctx, p, params);
        });
    this->add(
        "DissipatedEnergy",
        [](Context& ctx, NonLinearEvolutionProblemImplementation<false>& p,
           const Parameters& params) noexcept {
          return make_unique<DissipatedEnergyPostProcessing<false>>(ctx, p,
                                                                    params);
        });
  }  // end of PostProcessingFactory

  PostProcessingFactory<false>::~PostProcessingFactory() = default;

}  // end of namespace mfem_mgis
