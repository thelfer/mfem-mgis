/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx
 * \brief  This file declares the `NonLinearEvolutionProblemImplementation`
 * class
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
  struct AbstractNonLinearEvolutionProblemPostProcessing;

  /*!
   * \brief class for solving non linear evolution problems
   */
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

#ifdef MFEM_USE_MPI

  //! \brief parallel specialisation
  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemImplementation<true>
      : public NonLinearEvolutionProblemImplementationBase,
        public NonlinearForm<true> {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblemImplementation(
        Context& ctx,
        std::shared_ptr<FiniteElementDiscretization> fed,
        const Hypothesis h,
        const Parameters& p);
    //! \return the mesh
    [[nodiscard]] Mesh<true>& getMesh();
    //! \return the mesh
    [[nodiscard]] const Mesh<true>& getMesh() const;
    //! \return the finite element space
    [[nodiscard]] FiniteElementSpace<true>& getFiniteElementSpace();
    //! \return the finite element space
    [[nodiscard]] const FiniteElementSpace<true>& getFiniteElementSpace() const;
    /*!
     * \brief compute the residual
     * \param[in] u: current estimate of the unknowns
     * \param[out] r: residual
     */
    void Mult(const mfem::Vector& u, mfem::Vector& r) const override;
    //
    [[nodiscard]] bool addBoundaryCondition(
        Context& ctx,
        std::unique_ptr<AbstractDirichletBoundaryCondition> bc) noexcept
        override;
    [[nodiscard]] bool addBoundaryCondition(
        Context& ctx,
        std::unique_ptr<AbstractBoundaryCondition> f) noexcept override;
    /*!
     * \brief add a new post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: post-processing
     * \return true on success
     */
    [[nodiscard]] virtual bool addPostProcessing(
        Context& ctx,
        std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<true>>
            p) noexcept;
    //
    [[nodiscard]] bool integrate(const mfem::Vector& u,
                                 const IntegrationType it,
                                 const std::optional<real> odt) override;
    [[nodiscard]] bool setLinearSolver(Context& ctx,
                                       LinearSolverHandler s) noexcept override;
    [[nodiscard]] bool setLinearSolver(
        Context& ctx,
        std::string_view n,
        const Parameters& params) noexcept override;
    [[nodiscard]] bool addPostProcessing(
        Context& ctx,
        const std::function<void(const real, const real)>& p) noexcept override;
    [[nodiscard]] bool addPostProcessing(Context& ctx,
                                         std::string_view n,
                                         const Parameters& p) noexcept override;
    [[nodiscard]] bool executeInitialPostProcessings(
        Context& ctx, const real t) noexcept override;
    [[nodiscard]] bool executePostProcessings(Context&,
                                              const real t,
                                              const real dt) noexcept override;
    //! \brief destructor
    ~NonLinearEvolutionProblemImplementation() override;

   protected:
    //
    [[nodiscard]] std::optional<real> computePrediction(
        Context& ctx, const real t, const real dt) noexcept override;
    void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type> dofs) override;
    //! \brief registered post-processings
    std::vector<
        std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<true>>>
        postprocessings;
  };  // end of struct NonLinearEvolutionProblemImplementation

#endif /* MFEM_USE_MPI */

  //! \brief sequential specialisation
  template <>
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemImplementation<false>
      : public NonLinearEvolutionProblemImplementationBase,
        public NonlinearForm<false> {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblemImplementation(
        Context& ctx,
        std::shared_ptr<FiniteElementDiscretization> fed,
        const Hypothesis h,
        const Parameters& p);
    //! \return the mesh
    [[nodiscard]] Mesh<false>& getMesh() noexcept;
    //! \return the mesh
    [[nodiscard]] const Mesh<false>& getMesh() const noexcept;
    //! \return the finite element space
    [[nodiscard]] FiniteElementSpace<false>& getFiniteElementSpace() noexcept;
    //! \return the finite element space
    [[nodiscard]] const FiniteElementSpace<false>& getFiniteElementSpace()
        const noexcept;
    //
    [[nodiscard]] bool addBoundaryCondition(
        Context& ctx,
        std::unique_ptr<AbstractDirichletBoundaryCondition> bc) noexcept
        override;
    [[nodiscard]] bool addBoundaryCondition(
        Context& ctx,
        std::unique_ptr<AbstractBoundaryCondition> f) noexcept override;
    /*!
     * \brief add a new post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: post-processing
     * \return true on success
     */
    [[nodiscard]] virtual bool addPostProcessing(
        Context& ctx,
        std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<false>>
            p) noexcept;
    //
    [[nodiscard]] bool setLinearSolver(Context& ctx,
                                       LinearSolverHandler s) noexcept override;
    [[nodiscard]] bool setLinearSolver(
        Context& ctx,
        std::string_view n,
        const Parameters& params) noexcept override;
    [[nodiscard]] bool integrate(const mfem::Vector& u,
                                 const IntegrationType it,
                                 const std::optional<real> odt) override;
    [[nodiscard]] bool addPostProcessing(
        Context& ctx,
        const std::function<void(const real, const real)>& p) noexcept override;
    [[nodiscard]] bool addPostProcessing(Context& ctx,
                                         std::string_view n,
                                         const Parameters& p) noexcept override;
    [[nodiscard]] bool executeInitialPostProcessings(
        Context& ctx, const real t) noexcept override;
    [[nodiscard]] bool executePostProcessings(Context&,
                                              const real t,
                                              const real dt) noexcept override;
    //! \brief destructor
    ~NonLinearEvolutionProblemImplementation() override;

   protected:
    //
    [[nodiscard]] std::optional<real> computePrediction(
        Context& ctx, const real t, const real dt) noexcept override;
    void markDegreesOfFreedomHandledByDirichletBoundaryConditions(
        std::vector<size_type> dofs) override;
    //! \brief registered post-processings
    std::vector<
        std::unique_ptr<AbstractNonLinearEvolutionProblemPostProcessing<false>>>
        postprocessings;
  };  // end of struct NonLinearEvolutionProblemImplementation

  /*!
   * \brief compute the resultant of the inner forces on the given boundary
   * \param[in, out] ctx: execution context
   * \param[out] F: resultant
   * \param[in] p: non linear evolution problem
   * \param[in] elements: a structure which gives for each element having at
   * least one node on the boundary the list of the nodes of this element on the
   * boundary.
   * \return true on success
   *
   * \note in parallel, the resultant is only the contribution of the given
   * process
   */
  template <bool parallel>
  bool computeResultantForceOnBoundary(
      Context& ctx,
      mfem::Vector& F,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const std::vector<
          std::pair<size_type, std::vector<std::vector<size_type>>>>&
          elements) noexcept;

  /*!
   * \brief compute the integral of the thermodynamic forces and the volume of
   * each material
   * \return the integral of the thermodynamic forces at the end of the time
   * step and the volume of each material.
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \note in parallel, the returned value is only the contribution of the given
   * process
   */
  template <bool parallel>
  std::optional<std::pair<std::vector<std::vector<real>>, std::vector<real>>>
  computeMeanThermodynamicForcesValues(
      Context& ctx,
      NonLinearEvolutionProblemImplementation<parallel>& p) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.ixx"

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMIMPLEMENTATION */
