/*!
 * \file   include/MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#ifndef LIB_MFEM_MGIS_MULTIMATERIALNONLINEARINTEGRATOR_HXX
#define LIB_MFEM_MGIS_MULTIMATERIALNONLINEARINTEGRATOR_HXX

#include <memory>
#include <vector>
#include "mfem/fem/nonlininteg.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  // forward declaration
  enum struct IntegrationType;
  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  struct AbstractBehaviourIntegrator;

  /*!
   * \brief base class for non linear integrators based on an MGIS' behaviours.
   * This class manages an mapping associating a material and its identifier
   */
  struct MFEM_MGIS_EXPORT [[nodiscard]] MultiMaterialNonLinearIntegrator final
      : public NonlinearFormIntegrator {
    //! \brief a simple alias
    using Behaviour = mgis::behaviour::Behaviour;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretisation
     * \param[in] h: modelling hypothesis
     */
    MultiMaterialNonLinearIntegrator(
        std::shared_ptr<const FiniteElementDiscretization> fed,
        const Hypothesis h);
    // MFEM API
    void AssembleElementVector(const mfem::FiniteElement& e,
                               mfem::ElementTransformation& tr,
                               const mfem::Vector& U,
                               mfem::Vector& F) override;

    void AssembleElementGrad(const mfem::FiniteElement& e,
                             mfem::ElementTransformation& tr,
                             const mfem::Vector& U,
                             mfem::DenseMatrix& K) override;
    /*!
     * \brief integrate the behaviour for the current estimate of the unknowns
     * at the end of the time step.
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     */
    [[nodiscard]] bool integrate(const mfem::FiniteElement& e,
                                 mfem::ElementTransformation& tr,
                                 const mfem::Vector& U,
                                 const IntegrationType it);
    //! \return the current time increment
    [[nodiscard]] real getTimeIncrement() const noexcept;
    /*!
     * \brief set the value of the time increment
     * \param[in] dt: time increment
     */
    void setTimeIncrement(const real dt);
    /*!
     * \brief method called before each resolution
     *
     * \param[in, out] ctx: execution context
     *
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    [[nodiscard]] bool setup(Context& ctx,
                             const real t,
                             const real dt) noexcept;
    /*!
     * \brief add a new behaviour integrator
     * \return the behaviour integrator identifier
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material id
     * \param[in] l: library name
     * \param[in] b: behaviour name
     * \param[in] params: additional parameters
     */
    [[nodiscard]] std::optional<size_type> addBehaviourIntegrator(
        Context& ctx,
        const std::string& n,
        const size_type m,
        const std::string& l,
        const std::string& b,
        const Parameters& params = {}) noexcept;
    /*!
     * \return the material with the given id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] OptionalReference<const Material> getMaterial(
        Context& ctx, const size_type m, const size_type b) const noexcept;
    /*!
     * \return the material with the given id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] OptionalReference<Material> getMaterial(
        Context& ctx, const size_type m, const size_type b) noexcept;
    /*!
     * \return the number of behaviour integrators associated with the given
     * material id
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     */
    [[nodiscard]] std::optional<size_type> getNumberOfBehaviourIntegrators(
        Context& ctx, const size_type m) const noexcept;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] OptionalReference<const AbstractBehaviourIntegrator>
    getBehaviourIntegrator(Context& ctx,
                           const size_type m,
                           const size_type b) const noexcept;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] OptionalReference<AbstractBehaviourIntegrator>
    getBehaviourIntegrator(Context& ctx,
                           const size_type m,
                           const size_type b) noexcept;
    /*!
     * \brief revert the internal state variables.
     *
     * The values of the internal state variables at the beginning of the time
     * step are copied on the values of the internal state variables at
     * end of the time step.
     */
    void revert();
    /*!
     * \brief update the internal state variables.
     *
     * The values of the internal state variables at the end of the time step
     * are copied on the values of the internal state variables at beginning of
     * the time step.
     */
    void update();
    /*!
     * \brief set the macroscropic gradients
     * \param[in] g: macroscopic gradients
     */
    void setMacroscopicGradients(std::span<const real> g);
    /*!
     * \return the list of material identifiers for which a behaviour
     * integrator has been defined.
     */
    std::vector<size_type> getAssignedMaterialsIdentifiers() const;
    /*!
     * \return linearized operators
     * \param[in] u: current estimate of the unknowns
     * \note: those linearized operators used the consistent tangent operators
     * and thermodynamic forces computed by the integration. The user is
     * responible for calling the behaviour integration before using those
     * operators
     */
    [[nodiscard]] LinearizedOperators getLinearizedOperators(
        const mfem::Vector& u);
    //! \brief destructor
    ~MultiMaterialNonLinearIntegrator() override;

   protected:
    //! \brief underlying finite element space
    const std::shared_ptr<const FiniteElementDiscretization> fe_discretization;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
    /*!
     * \brief mapping between the material identifiers and the behaviour
     * integrators.
     */
    std::vector<std::vector<std::unique_ptr<AbstractBehaviourIntegrator>>>
        behaviour_integrators;
  };  // end of MultiMaterialNonLinearIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MULTIMATERIALNONLINEARINTEGRATOR_HXX */
