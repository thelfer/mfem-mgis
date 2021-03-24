/*!
 * \file   src/PostProcessingFactory.cxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#include <utility>
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/PostProcessing.hxx"
#include "MFEMMGIS/ParaviewExportResults.hxx"
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
      mgis::raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  std::unique_ptr<PostProcessing<true>> PostProcessingFactory<true>::generate(
      std::string_view n,
      NonLinearEvolutionProblemImplementation<true>& p,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("PostProcessingFactory<true>::generate: ");
      msg += "no post-processing called '";
      msg += n;
      msg += "' declared";
      mgis::raise(msg);
    }
    const auto& g = pg->second;
    return g(p, params);
  }  // end of generate

  PostProcessingFactory<true>::PostProcessingFactory() {
    this->add("ParaviewExportResults",
              [](NonLinearEvolutionProblemImplementation<true>& p,
                 const Parameters& params) {
                return std::make_unique<ParaviewExportResults<true>>(p, params);
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
      mgis::raise(msg);
    }
    this->generators.insert({std::string(n), std::move(g)});
  }  // end of add

  std::unique_ptr<PostProcessing<false>> PostProcessingFactory<false>::generate(
      std::string_view n,
      NonLinearEvolutionProblemImplementation<false>& p,
      const Parameters& params) const {
    const auto pg = this->generators.find(n);
    if (pg == this->generators.end()) {
      std::string msg("PostProcessingFactory<false>::generate: ");
      msg += "no post-processing called '";
      msg += n;
      msg += "' declared";
      mgis::raise(msg);
    }
    const auto& g = pg->second;
    return g(p, params);
  }  // end of generate

  PostProcessingFactory<false>::PostProcessingFactory(){
    this->add("ParaviewExportResults",
              [](NonLinearEvolutionProblemImplementation<false>& p,
                 const Parameters& params) {
                return std::make_unique<ParaviewExportResults<false>>(p,
                                                                      params);
              });
  }  // end of PostProcessingFactory

  PostProcessingFactory<false>::~PostProcessingFactory() = default;

}  // end of namespace mfem_mgis
