/*!
 * \file   MultiMaterialEvolutionProblemBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#ifndef LIB_MFEM_MGIS_MULTIMATERIALEVOLUTIONPROBLEMBASE_HXX
#define LIB_MFEM_MGIS_MULTIMATERIALEVOLUTIONPROBLEMBASE_HXX

#include <memory>
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  struct Material;
  // forward declaration
  struct MultiMaterialNonLinearIntegratorBase;
  // forward declaration
  template <bool parallel>
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
     * \return the material with the given id
     * \param[in] m: material id
     */
    const Material& getMaterial(const size_type) const;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    Material& getMaterial(const size_type);
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
     * \brief return the underlying base class for the active multi-material
     * integrator.
     */
    MultiMaterialNonLinearIntegratorBase&
    getMultiMaterialNonLinearIntegratorBase();

    /*!
     * \brief return the underlying base class for the active multi-material
     * integrator.
     */
    const MultiMaterialNonLinearIntegratorBase&
    getMultiMaterialNonLinearIntegratorBase() const;
    //
    void revert();
    void update();
    void setTimeIncrement(const real);
    void setup();

#ifdef MFEM_USE_MPI
    //! \brief pointer to the underlying domain integrator
    MultiMaterialNonLinearIntegrator<true>* const parallel_mgis_integrator;
#endif /* MFEM_USE_MPI */
    //! \brief pointer to the underlying domain integrator
    MultiMaterialNonLinearIntegrator<false>* const sequential_mgis_integrator;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
  };  // namespace mfem_mgis

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MULTIMATERIALEVOLUTIONPROBLEMBASE_HXX */
