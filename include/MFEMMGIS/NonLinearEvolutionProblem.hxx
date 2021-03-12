/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEM_HXX

#include <memory>
#include <vector>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MultiMaterialEvolutionProblemBase.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemBase.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct PostProcessing;

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
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(std::unique_ptr<PostProcessing<true>>);
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(
        const std::function<void(const real, const real)>&);
    /*!
     * \brief execute the registred postprocessings
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void executePostProcessings(const real, const real);
    //
    void revert() override;
    void update() override;
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   protected:
    void setup(const real, const real) override;

   private:
    void setTimeIncrement(const real) override;
    //! \brief registred post-processings
    std::vector<std::unique_ptr<PostProcessing<true>>> postprocessings;
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
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(std::unique_ptr<PostProcessing<false>>);
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(
        const std::function<void(const real, const real)>&);
    /*!
     * \brief execute the registred postprocessings
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void executePostProcessings(const real, const real);
    //
    void revert() override;
    void update() override;
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   protected:
    void setup(const real, const real) override;

   private:
    void setTimeIncrement(const real) override;
    //! \brief registred post-processings
    std::vector<std::unique_ptr<PostProcessing<false>>> postprocessings;
  };  // end of struct NonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEM */
