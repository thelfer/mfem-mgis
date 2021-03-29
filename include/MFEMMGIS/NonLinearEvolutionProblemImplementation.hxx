/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_HXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION_HXX

#include <memory>
#include <vector>
#include "mfem/fem/nonlinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pnonlinearform.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementationBase.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct PostProcessing;

  /*!
   * \brief class for solving non linear evolution problems
   */
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

#ifdef MFEM_USE_MPI

  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemImplementation<true>
      : public NonLinearEvolutionProblemImplementationBase,
        public NonlinearForm<true> {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblemImplementation(
        std::shared_ptr<FiniteElementDiscretization>,
        const Hypothesis,
        const Parameters&);
    //! \return the finite element space
    FiniteElementSpace<true>& getFiniteElementSpace();
    //! \return the finite element space
    const FiniteElementSpace<true>& getFiniteElementSpace() const;
    // ! \brief return the underlying newton solver
    virtual NewtonSolver& getSolver();
    /*!
     * \brief set the linear solver
     * \param[in] s: linear solver
     */
    virtual void setLinearSolver(std::unique_ptr<LinearSolver>);
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(std::unique_ptr<PostProcessing<true>>);
    //
    void setSolverParameters(const Parameters&) override;
    void setLinearSolver(std::string_view, const Parameters&) override;
    void addPostProcessing(
        const std::function<void(const real, const real)>&) override;
    void addPostProcessing(std::string_view, const Parameters&) override;
    void executePostProcessings(const real, const real) override;
    void solve(const real, const real) override;
    //! \brief destructor
    ~NonLinearEvolutionProblemImplementation() override;

   protected:
    //
    void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type>) override;

   private:
    //! \brief newton solver
    NewtonSolver solver;
    //! \brief linear solver
    std::unique_ptr<LinearSolver> linear_solver;
    //! \brief registred post-processings
    std::vector<std::unique_ptr<PostProcessing<true>>> postprocessings;
  };  // end of struct NonLinearEvolutionProblemImplementation

#endif /* MFEM_USE_MPI */

  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemImplementation<false>
      : public NonLinearEvolutionProblemImplementationBase,
        public NonlinearForm<false> {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblemImplementation(
        std::shared_ptr<FiniteElementDiscretization>,
        const Hypothesis,
        const Parameters&);
    //! \return the finite element space
    FiniteElementSpace<false>& getFiniteElementSpace();
    //! \return the finite element space
    const FiniteElementSpace<false>& getFiniteElementSpace() const;
    // ! \brief return the underlying newton solver
    virtual NewtonSolver& getSolver();
    /*!
     * \brief set the linear solver
     * \param[in] s: linear solver
     */
    virtual void setLinearSolver(std::unique_ptr<LinearSolver>);
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(std::unique_ptr<PostProcessing<false>>);
    //
    void setSolverParameters(const Parameters&) override;
    void setLinearSolver(std::string_view, const Parameters&) override;
    void addPostProcessing(
        const std::function<void(const real, const real)>&) override;
    void addPostProcessing(std::string_view, const Parameters&) override;
    void executePostProcessings(const real, const real) override;
    void solve(const real, const real) override;
    //! \brief destructor
    ~NonLinearEvolutionProblemImplementation() override;

   protected:
    //
    void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type>) override;

   private:
    //! \brief newton solver
    NewtonSolver solver;
    //! \brief linear solver
    std::unique_ptr<LinearSolver> linear_solver;
    //! \brief registred post-processings
    std::vector<std::unique_ptr<PostProcessing<false>>> postprocessings;
  };  // end of struct NonLinearEvolutionProblemImplementation

  /*!
   * \brief return the resultant of the inner forces on the given boundary
   * \param[out] F: resultant
   * \param[in] p: non linear evolution problem
   * \param[in] faces: description of the boundary by a vector of pair
   * associating for each face its identifier and the identifier of the adjacent
   * element.
   *
   * \note in parallel, the resultant is only the contribution of the given
   * process
   */
  template <bool parallel>
  void computeResultantForceOnBoundary(
      mfem::Vector&,
      NonLinearEvolutionProblemImplementation<parallel>&,
      const std::vector<std::pair<size_type, size_type>>&);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx"

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION */