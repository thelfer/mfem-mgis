/*!
 * \file   include/MFEMMGIS/BehaviourIntegratorFactory.hxx
 * \brief
 * \author Thomas Helfer
 * \date   13/10/2020
 */

#ifndef LIB_MFEMMGIS_BEHAVIOURINTEGRATORFACTORY_HXX
#define LIB_MFEMMGIS_BEHAVIOURINTEGRATORFACTORY_HXX

#include <map>
#include <memory>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;

  /*!
   * \brief an abstract factory for behaviour integrators
   */
  struct MFEM_MGIS_EXPORT BehaviourIntegratorFactory {
    //! a simple alias
    using Generator =
        std::function<std::unique_ptr<AbstractBehaviourIntegrator>(
            const FiniteElementDiscretization&,
            const size_type,
            std::unique_ptr<const Behaviour>)>;
    /*!
     * \return the unique instance of this class for the given hypothesis
     * \param[in] h: modelling hypothesis
     */
    static BehaviourIntegratorFactory& get(const Hypothesis);
    //! \brief constructor
    BehaviourIntegratorFactory();
    //! \brief move constructor
    BehaviourIntegratorFactory(BehaviourIntegratorFactory&&);
    //! \brief copy constructor
    BehaviourIntegratorFactory(const BehaviourIntegratorFactory&);
    //! \brief move assignement
    BehaviourIntegratorFactory& operator=(BehaviourIntegratorFactory&&);
    //! \brief standard assignement
    BehaviourIntegratorFactory& operator=(const BehaviourIntegratorFactory&);
    /*!
     * \brief register a new behaviour integrator
     * \param[in] n: name
     * \param[in] g: generator
     */
    void addGenerator(const std::string&, const Generator);
    /*!
     * \return a newly created behaviour integrator
     * \param[in] n: name
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b: behaviour
     */
    std::unique_ptr<AbstractBehaviourIntegrator> generate(
        const std::string&,
        const FiniteElementDiscretization&,
        const size_type,
        std::unique_ptr<const Behaviour>) const;
    //! \brief destructor
    ~BehaviourIntegratorFactory() noexcept;

   private:
    /*!
     * \return default initialised factories for the supported modelling
     * hypotheses.
     */
    static std::map<Hypothesis, BehaviourIntegratorFactory> buildFactories();
    //! \brief registred generators
    std::map<std::string, Generator> generators;

  };  // end of BehaviourIntegratorFactory

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_BEHAVIOURINTEGRATORFACTORY_HXX */
