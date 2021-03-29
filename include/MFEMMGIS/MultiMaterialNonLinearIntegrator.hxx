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

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  struct BehaviourIntegrator;

  /*!
   * \brief base class for non linear integrators based on an MGIS' behaviours.
   * This class manages an mapping associating a material and its identifier
   */
  struct MFEM_MGIS_EXPORT MultiMaterialNonLinearIntegrator final
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
     */
    virtual bool integrate(const mfem::FiniteElement &,
                           mfem::ElementTransformation &,
                           const mfem::Vector &);
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
     * \brief add a new material
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material id
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    virtual void addBehaviourIntegrator(const std::string &,
                                        const size_type,
                                        const std::string &,
                                        const std::string &);
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual const Material &getMaterial(const size_type) const;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual Material &getMaterial(const size_type);
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in] m: material id
     */
    virtual const BehaviourIntegrator &getBehaviourIntegrator(
        const size_type) const;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in] m: material id
     */
    virtual BehaviourIntegrator &getBehaviourIntegrator(const size_type);
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
    virtual void setMacroscopicGradients(mgis::span<const real>);
    /*!
     * \return the list of material identifiers for which a behaviour
     * integrator has been defined.
     */
    virtual std::vector<size_type> getMaterialIdentifiers() const;

    //! \brief destructor
    virtual ~MultiMaterialNonLinearIntegrator();

   protected:
    //! \brief underlying finite element space
    const std::shared_ptr<const FiniteElementDiscretization> fe_discretization;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
    /*!
     * \brief mapping between the material integrator and the behaviour
     * integrator.
     */
    std::vector<std::unique_ptr<BehaviourIntegrator>> behaviour_integrators;
  };  // end of MultiMaterialNonLinearIntegrator

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MULTIMATERIALNONLINEARINTEGRATOR_HXX */
