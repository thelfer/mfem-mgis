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
#include "MFEMMGIS/Utilities/OptionalReference.hxx"
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
        std::shared_ptr<const FiniteElementDiscretization>, const Hypothesis);
    // MFEM API
    void AssembleElementVector(const mfem::FiniteElement &,
                               mfem::ElementTransformation &,
                               const mfem::Vector &,
                               mfem::Vector &) override;

    void AssembleElementGrad(const mfem::FiniteElement &,
                             mfem::ElementTransformation &,
                             const mfem::Vector &,
                             mfem::DenseMatrix &) override;
    /*!
     * \brief integrate the behaviour for the current estimate of the unknowns
     * at the end of the time step.
     * \param[in] e: finite element
     * \param[in] tr: finite element transformation
     * \param[in] u: current estimate of the unknowns
     * \param[in] it: integration type
     */
    virtual bool integrate(const mfem::FiniteElement &,
                           mfem::ElementTransformation &,
                           const mfem::Vector &,
                           const IntegrationType);
    /*!
     * \brief set the value of the time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real);
    /*!
     * \brief method called before each resolution
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void setup(const real, const real);
    /*!
     * \brief add a new behaviour integrator
     * \return the behaviour integrator identifier
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material id
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    [[nodiscard]] virtual size_type addBehaviourIntegrator(const std::string &,
                                                           const size_type,
                                                           const std::string &,
                                                           const std::string &);
    /*!
     * \return the material with the given id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] virtual OptionalReference<const Material> getMaterial(
        Context &, const size_type, const size_type) const noexcept;
    /*!
     * \return the material with the given id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] virtual OptionalReference<Material> getMaterial(
        Context &, const size_type, const size_type) noexcept;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] virtual OptionalReference<const AbstractBehaviourIntegrator>
    getBehaviourIntegrator(Context &,
                           const size_type,
                           const size_type) const noexcept;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in, out] ctx: execution context
     * \param[in] m: material id
     * \param[in] b: behaviour id
     */
    [[nodiscard]] virtual OptionalReference<AbstractBehaviourIntegrator>
    getBehaviourIntegrator(Context &,
                           const size_type,
                           const size_type) noexcept;
    /*!
     * \brief revert the internal state variables.
     *
     * The values of the internal state variables at the beginning of the time
     * step are copied on the values of the internal state variables at
     * end of the time step.
     */
    virtual void revert();
    /*!
     * \brief update the internal state variables.
     *
     * The values of the internal state variables at the end of the time step
     * are copied on the values of the internal state variables at beginning of
     * the time step.
     */
    virtual void update();
    /*!
     * \brief set the macroscropic gradients
     * \param[in] g: macroscopic gradients
     */
    virtual void setMacroscopicGradients(std::span<const real>);
    /*!
     * \return the list of material identifiers for which a behaviour
     * integrator has been defined.
     */
    virtual std::vector<size_type> getAssignedMaterialsIdentifiers() const;
    /*!
     * \return linearized operators
     * \param[in] u: current estimate of the unknowns
     * \note: those linearized operators used the consistent tangent operators
     * and thermodynamic forces computed by the integration. The user is
     * responible for calling the behaviour integration before using those
     * operators
     */
    [[nodiscard]] virtual LinearizedOperators getLinearizedOperators(
        const mfem::Vector &);
    /*!
     * \return the material with the given id for the first behaviour integrator
     * \param[in] m: material id
     */
    [[nodiscard, deprecated]] virtual const Material &getMaterial(
        const size_type) const;
    /*!
     * \return the material with the given id for the first behaviour integrator
     * \param[in] m: material id
     */
    [[nodiscard, deprecated]] virtual Material &getMaterial(const size_type);
    /*!
     * \return the first behaviour integrator with the given material id
     * \param[in] m: material id
     */
    [[nodiscard, deprecated]] virtual const AbstractBehaviourIntegrator &
    getBehaviourIntegrator(const size_type) const;
    /*!
     * \return the first behaviour integrator with the given material id
     * \param[in] m: material id
     */
    [[nodiscard, deprecated]] virtual AbstractBehaviourIntegrator &
    getBehaviourIntegrator(const size_type);
    //! \brief destructor
    virtual ~MultiMaterialNonLinearIntegrator();

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
