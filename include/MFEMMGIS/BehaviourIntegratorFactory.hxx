/*!
 * \file   include/MFEMMGIS/BehaviourIntegratorFactory.hxx
 * \brief  This file declares the `BehaviourIntegratorFactory` class
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
   * \param[in, out] f: factory
   *
   * \note this function is meant to be specialised to declare behaviour
   * integrators which are only valid for this modelling hypothesis. Behaviour
   * integrators defined for all modelling hypotheses are declared by the
   * `fillWithDefaultBehaviourIntegrators`.
   */
  template <Hypothesis H>
  void buildFactory(BehaviourIntegratorFactory& f);
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
     * \brief return the unique instance of this class for the given hypothesis
     * \param[in, out] ctx: execution context
     * \param[in] h: modelling hypothesis
     * \return the factory associated with the given hypothesis
     */
    static OptionalReference<BehaviourIntegratorFactory> get(
        Context& ctx, const Hypothesis h) noexcept;
    /*!
     * \brief return the unique instance of this class for the given hypothesis
     * \param[in] h: modelling hypothesis
     * \return the factory associated with the given hypothesis
     */
    static BehaviourIntegratorFactory& get(const Hypothesis h);
    /*!
     * \brief register a new behaviour integrator
     * \param[in, out] ctx: execution context
     * \param[in] n: name
     * \param[in] g: generator
     * \return true on success
     */
    [[nodiscard]] bool addGenerator(Context& ctx,
                                    const std::string& n,
                                    const Generator g) noexcept;
    /*!
     * \brief create a new behaviour integrator
     * \param[in, out] ctx: execution context
     * \param[in] n: name
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b: behaviour
     * \param[in] params: additional parameters
     * \return a newly created behaviour integrator
     */
    [[nodiscard, deprecated]] std::unique_ptr<AbstractBehaviourIntegrator>
    generate(Context& ctx,
             std::string_view n,
             const FiniteElementDiscretization& fed,
             const size_type m,
             std::unique_ptr<const Behaviour> b,
             const Parameters& params = Parameters{}) const noexcept;
    /*!
     * \brief create a new behaviour integrator
     * \param[in] n: name
     * \param[in] fed: finite element discretization.
     * \param[in] m: material attribute.
     * \param[in] b: behaviour
     * \return a newly created behaviour integrator
     */
    [[nodiscard, deprecated]] std::unique_ptr<AbstractBehaviourIntegrator>
    generate(std::string_view n,
             const FiniteElementDiscretization& fed,
             const size_type m,
             std::unique_ptr<const Behaviour> b) const;
    //! \brief destructor
    ~BehaviourIntegratorFactory() noexcept;

   private:
    //! \brief constructor
    BehaviourIntegratorFactory();
    //! \brief move constructor
    BehaviourIntegratorFactory(BehaviourIntegratorFactory&&);
    //! \brief copy constructor
    BehaviourIntegratorFactory(const BehaviourIntegratorFactory&) = delete;
    //! \brief move assignment
    BehaviourIntegratorFactory& operator=(BehaviourIntegratorFactory&&) =
        delete;
    //! \brief standard assignment
    BehaviourIntegratorFactory& operator=(const BehaviourIntegratorFactory&) =
        delete;
    /*!
     * \brief register a new behaviour integrator
     * \param[in] a: attribute stating that this method may abort on error
     * \param[in] n: name
     * \param[in] g: generator
     */
    void addGenerator(attributes::MayAbort a,
                      std::string_view n,
                      const Generator g);
    /*!
     * \brief register a new behaviour integrator
     * \param[in] a: attribute stating that this method may abort on error
     * \param[in] n: name
     * \param[in] g: generator
     */
    void addGenerator(attributes::MayAbort a,
                      std::string_view n,
                      const DeprecatedGeneratorType g);
    //
    template <Hypothesis H>
    friend void buildFactory(BehaviourIntegratorFactory& f);
    /*!
     * \brief build the factory associated with the given hypothesis and add
     * it to the given map
     * \tparam H: modelling hypothesis
     * \param[in, out] factories: list of factories per modelling hypotheses
     */
    template <Hypothesis H>
    static void addFactory(
        std::map<Hypothesis, std::unique_ptr<BehaviourIntegratorFactory>>&
            factories);
    /*!
     * \return default initialised factories for the supported modelling
     * hypotheses.
     */
    static std::map<Hypothesis, std::unique_ptr<BehaviourIntegratorFactory>>
    buildFactories();
    //! \brief registered generators
    std::map<std::string, Generator, std::less<>> generators;

  };  // end of BehaviourIntegratorFactory

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_BEHAVIOURINTEGRATORFACTORY_HXX */
