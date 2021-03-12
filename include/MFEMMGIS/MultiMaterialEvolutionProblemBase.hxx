/*!
 * \file   MultiMaterialEvolutionProblemBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#ifndef LIB_MFEM_MGIS_MULTIMATERIALEVOLUTIONPROBLEMBASE_HXX
#define LIB_MFEM_MGIS_MULTIMATERIALEVOLUTIONPROBLEMBASE_HXX

#include <memory>
#include "MGIS/Span.hxx"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  struct Material;
  // forward declaration
  struct MultiMaterialNonLinearIntegrator;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT MultiMaterialEvolutionProblemBase {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     */
    MultiMaterialEvolutionProblemBase(
        std::shared_ptr<FiniteElementDiscretization>, const Hypothesis);
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
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual const Material& getMaterial(const size_type) const;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual Material& getMaterial(const size_type);
    /*!
     * \brief add a new material
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material id
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    virtual void addBehaviourIntegrator(const std::string&,
                                        const size_type,
                                        const std::string&,
                                        const std::string&);
    //! \brief destructor
    virtual ~MultiMaterialEvolutionProblemBase();

   protected:
    /*!
     * \brief revert the state of the materials at the beginning of the time
     * step, i.e. the state at the end of the time step is reset using the
     * state at the beginning of the time step.
     */
    void revert();
    /*!
     * \brief update the state of the materials, i.e. the state at the end of
     * the time step becomes the state at the beginning of the next time step.
     */
    void update();
    /*!
     * \brief set the time increment for the current time step
     * \param[in] dt: time increment
     */
    void setTimeIncrement(const real);
    /*!
     * \brief method called at the beginning of each resolution.
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    void setup(const real, const real);
    //! \brief pointer to the underlying domain integrator
    MultiMaterialNonLinearIntegrator* const mgis_integrator;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
  };  // namespace mfem_mgis

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MULTIMATERIALEVOLUTIONPROBLEMBASE_HXX */
