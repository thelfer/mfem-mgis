/*!
 * \file   src/Behaviour.cxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include "MGIS/Behaviour/FiniteStrainBehaviourOptions.hxx"
#include "MFEMMGIS/Behaviour.hxx"

namespace mfem_mgis {

  std::unique_ptr<Behaviour> load(Context& ctx,
                                  const std::string& l,
                                  const std::string& b,
                                  const Hypothesis h) noexcept {
    using namespace mgis::behaviour;
    try {
      if (isStandardFiniteStrainBehaviour(l, b)) {
        auto opts = FiniteStrainBehaviourOptions{};
        opts.stress_measure = FiniteStrainBehaviourOptions::PK1;
        opts.tangent_operator = FiniteStrainBehaviourOptions::DPK1_DF;
        return std::make_unique<Behaviour>(
            mgis::behaviour::load(opts, l, b, h));
      }
      return std::make_unique<Behaviour>(mgis::behaviour::load(l, b, h));
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of load

}  // end of namespace mfem_mgis
