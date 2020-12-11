/*!
 * \file   src/Behaviour.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#include "MGIS/Behaviour/FiniteStrainBehaviourOptions.hxx"
#include "MFEMMGIS/Behaviour.hxx"

namespace mfem_mgis{

  std::shared_ptr<Behaviour> load(const std::string& l,
                                  const std::string& b,
                                  const Hypothesis h) {
    using namespace mgis::behaviour;
    if (isStandardFiniteStrainBehaviour(l, b)) {
      auto opts = FiniteStrainBehaviourOptions{};
      opts.stress_measure = FiniteStrainBehaviourOptions::PK1;
      opts.tangent_operator = FiniteStrainBehaviourOptions::DSIG_DF;
      return std::make_shared<Behaviour>(mgis::behaviour::load(opts, l, b, h));
    }
    return std::make_shared<Behaviour>(mgis::behaviour::load(l, b, h));
  }  // end of load

} // end of namespace mfem_mgis
