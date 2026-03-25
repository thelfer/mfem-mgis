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
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractBehaviourIntegrator.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  struct BehaviourIntegratorFactory;
  /*!
   * \brief fill the given factory for the given hypothesis with some
   * predeclared generators.
   *
   * \tparam H: modelling hypothesis
   * \param[in] f: factory
   *
   * \note this function is meant to be overriden to declare behaviour
   * integrators which are only valid for this modelling hypothesis. Behaviour
   * integrators defined for all modelling hypotheses are declared by the
   * `fillWithDefaultBehaviourIntegrators`.
   */
  template <Hypothesis H>
  void buildFactory(BehaviourIntegratorFactory&);
  /*!
   * \brief an abstract factory for behaviour integrators
   */
  struct MFEM_MGIS_EXPORT BehaviourIntegratorFactory {
    //! a simple alias
    using Generator =
        std::function<std::unique_ptr<AbstractBehaviourIntegrator>(
            Context&,
            const FiniteElementDiscretization&,
            const size_type,
            std::unique_ptr<const Behaviour>,
            const Parameters&)>;
    //! a simple alias, kept for backward compatibility
    using DeprecatedGeneratorType =
        std::function<std::unique_ptr<AbstractBehaviourIntegrator>(
            const FiniteElementDiscretization&,
            const size_type,
            std::unique_ptr<const Behaviour>)>;
    /*!
     * \return the unique instance of this class for the given hypothesis
     * \param[in] ctx: execution context
     * \param[in] h: modelling hypothesis
     */
    static OptionalReference<BehaviourIntegratorFactory> get(
        Context&, const Hypothesis) noexcept;
    /*!
     * \return the unique instance of this class for the given hypothesis
     * \param[in] h: modelling hypothesis
     */
    static BehaviourIntegratorFactory& get(const Hypothesis);
    /*!
     * \brief register a new behaviour integrator
     * \param[in] ctx: execution context
     * \param[in] n: name
     * \param[in] g: generator
     */
    [[nodiscard]] bool addGenerator(Context&,
                                    const std::string&,
                                    const Generator) noexcept;
    /*!
     * \return a newly created behaviour integrator
     * \param[in] ctx: execution context
     * \param[in] n: name
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b: behaviour
     * \param[in] params: additional parameters
     */
    [[nodiscard, deprecated]] std::unique_ptr<AbstractBehaviourIntegrator>
    generate(Context&,
             std::string_view,
             const FiniteElementDiscretization&,
             const size_type,
             std::unique_ptr<const Behaviour>,
             const Parameters& = Parameters{}) const noexcept;
    /*!
     * \return a newly created behaviour integrator
     * \param[in] n: name
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b: behaviour
     */
    [[nodiscard, deprecated]] std::unique_ptr<AbstractBehaviourIntegrator>
    generate(std::string_view,
             const FiniteElementDiscretization&,
             const size_type,
             std::unique_ptr<const Behaviour>) const;
    //! \brief destructor
    ~BehaviourIntegratorFactory() noexcept;

   private:
    //! \brief constructor
    BehaviourIntegratorFactory();
    //! \brief move constructor
    BehaviourIntegratorFactory(BehaviourIntegratorFactory&&);
    //! \brief copy constructor
    BehaviourIntegratorFactory(const BehaviourIntegratorFactory&) = delete;
    //! \brief move assignement
    BehaviourIntegratorFactory& operator=(BehaviourIntegratorFactory&&) =
        delete;
    //! \brief standard assignement
    BehaviourIntegratorFactory& operator=(const BehaviourIntegratorFactory&) =
        delete;
    /*!
     * \brief register a new behaviour integrator
     * \param[in] a: attribute expliciting that this method may abort on error
     * \param[in] n: name
     * \param[in] g: generator
     */
    void addGenerator(attributes::MayAbort, std::string_view, const Generator);
    /*!
     * \brief register a new behaviour integrator
     * \param[in] a: attribute expliciting that this method may abort on error
     * \param[in] n: name
     * \param[in] g: generator
     */
    void addGenerator(attributes::MayAbort,
                      std::string_view,
                      const DeprecatedGeneratorType);
    //
    template <Hypothesis H>
    friend void buildFactory(BehaviourIntegratorFactory&);
    /*!
     * \brief an helper function which add the behaviour integrators for the
     * given hypothesis \tparam H: modelling hypothesis \param[in, out]
     * factories: list of factories per modelling hypotheses
     */
    template <Hypothesis H>
    static void addFactory(
        std::map<Hypothesis, std::unique_ptr<BehaviourIntegratorFactory>>&);
    /*!
     * \return default initialised factories for the supported modelling
     * hypotheses.
     */
    static std::map<Hypothesis, std::unique_ptr<BehaviourIntegratorFactory>>
    buildFactories();
    //! \brief registred generators
    std::map<std::string, Generator, std::less<>> generators;

  };  // end of BehaviourIntegratorFactory

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_BEHAVIOURINTEGRATORFACTORY_HXX */
