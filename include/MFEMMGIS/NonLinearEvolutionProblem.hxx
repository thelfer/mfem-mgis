/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_EVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_EVOLUTIONPROBLEM_HXX

#include <memory>
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemBase.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Material;
  // forward declaration
  struct MultiMaterialNonLinearIntegrator;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblem
      : public NonLinearEvolutionProblemBase {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     */
    NonLinearEvolutionProblem(std::shared_ptr<FiniteElementDiscretization>,
                              const Hypothesis);
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
    //
    void revert() override;
    void update() override;
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   private:
    void setTimeIncrement(const real) override;
    //! \brief pointer to the underlying domain integrator
    MultiMaterialNonLinearIntegrator* const mgis_integrator;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
  };  // end of struct NonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_EVOLUTIONPROBLEM */
