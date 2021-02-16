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
#include "MFEMMGIS/MultiMaterialEvolutionProblemBase.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief class for solving non linear evolution problems
   */
  template <bool parallel>
  struct NonLinearEvolutionProblem;

#ifdef MFEM_USE_MPI

  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblem<true>
      : public NonLinearEvolutionProblemBase<true>,
        public MultiMaterialEvolutionProblemBase {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     */
    NonLinearEvolutionProblem(std::shared_ptr<FiniteElementDiscretization>,
                              const Hypothesis);
    //
    void revert() override;
    void update() override;
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   private:
    void setTimeIncrement(const real) override;
    void setup() override;
  };  // end of struct NonLinearEvolutionProblem

#endif /* MFEM_USE_MPI */

  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblem<false>
      : public NonLinearEvolutionProblemBase<false>,
        public MultiMaterialEvolutionProblemBase {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     */
    NonLinearEvolutionProblem(std::shared_ptr<FiniteElementDiscretization>,
                              const Hypothesis);
    //
    void revert() override;
    void update() override;
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   private:
    void setTimeIncrement(const real) override;
    void setup() override;
  };  // end of struct NonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_EVOLUTIONPROBLEM */
