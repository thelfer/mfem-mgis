/*!
 * \file   include/MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx
 * \brief  This file declares the `PeriodicNonLinearEvolutionProblem` class
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX

#include <span>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  //! \brief coordinate minimized to select the node whose displacement is fixed
  enum BoundaryConditionType {
    FIX_XMIN = 0,  //!< x coordinate
    FIX_YMIN = 1,  //!< y coordinate
    FIX_ZMIN = 2   //!< z coordinate
  };

#ifdef MFEM_USE_MPI

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in, out] ctx: execution context
   * \param[in, out] p: problem
   * \param[in] corner1: one corner of the computation domain
   * \param[in] corner2: second corner of the computation domain
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      mgis::Context& ctx,
      NonLinearEvolutionProblemImplementation<true>& p,
      const std::span<const real>& corner1,
      const std::span<const real>& corner2);

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in, out] ctx: execution context
   * \param[in, out] p: problem
   * \param[in] bct: impose the zero value on displacement field at xmin, or
   * ymin, or zmin
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      mgis::Context& ctx,
      NonLinearEvolutionProblemImplementation<true>& p,
      const mfem_mgis::BoundaryConditionType bct = mfem_mgis::FIX_XMIN);

#endif /* MFEM_USE_MPI */

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in, out] ctx: execution context
   * \param[in, out] p: problem
   * \param[in] corner1: one corner of the computation domain
   * \param[in] corner2: second corner of the computation domain
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      mgis::Context& ctx,
      NonLinearEvolutionProblemImplementation<false>& p,
      const std::span<const real>& corner1,
      const std::span<const real>& corner2);

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in, out] ctx: execution context
   * \param[in, out] p: problem
   * \param[in] bct: impose the zero value on displacement field at xmin, or
   * ymin, or zmin
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      mgis::Context& ctx,
      NonLinearEvolutionProblemImplementation<false>& p,
      const mfem_mgis::BoundaryConditionType bct = mfem_mgis::FIX_XMIN);

  /*!
   * \brief compute the squared distance from the node identified in vector
   * `nodes` at index `index` to the closest corner of the box defined by
   * `corner1` and `corner2`
   * \param[in] nodes: nodes coordinates
   * \param[in] reorder_space: true if the coordinates are ordered by
   * `mfem::Ordering::byNODES`
   * \param[in] dim: space dimension
   * \param[in] index: index of the node
   * \param[in] size: number of nodes
   * \param[in] corner1: one corner of the computation domain
   * \param[in] corner2: second corner of the computation domain
   * \return the squared distance
   */
  MFEM_MGIS_EXPORT real getNodesDistance(const mfem::GridFunction& nodes,
                                         const bool reorder_space,
                                         const size_t dim,
                                         const int index,
                                         const int size,
                                         const std::span<const real>& corner1,
                                         const std::span<const real>& corner2);

  /*!
   * \brief a base class handling the evolution of the macroscopic gradients
   */
  struct MFEM_MGIS_EXPORT PeriodicNonLinearEvolutionProblem
      : NonLinearEvolutionProblem {
    /*!
     * \brief constructor with profiling support
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] corner1: one corner of the computation domain
     * \param[in] corner2: second corner of the computation domain
     */
    PeriodicNonLinearEvolutionProblem(
        mgis::Context& ctx,
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::span<const real>& corner1,
        const std::span<const real>& corner2);

    /*!
     * \brief constructor with profiling support
     * \param[in, out] ctx: execution context
     * \param[in] fed: finite element discretization
     * \param[in] bct: impose the zero value on displacement field at xmin, or
     * ymin, or zmin
     */
    PeriodicNonLinearEvolutionProblem(
        mgis::Context& ctx,
        std::shared_ptr<FiniteElementDiscretization> fed,
        const mfem_mgis::BoundaryConditionType bct = mfem_mgis::FIX_XMIN);
    //! \brief move constructor
    PeriodicNonLinearEvolutionProblem(
        PeriodicNonLinearEvolutionProblem&&) noexcept = default;
    //! \brief deleted copy constructor
    PeriodicNonLinearEvolutionProblem(
        const PeriodicNonLinearEvolutionProblem&) noexcept = delete;
    [[nodiscard]] bool addBoundaryCondition(
        Context& ctx,
        std::unique_ptr<AbstractBoundaryCondition> f) noexcept override;
    /*!
     * \brief always fails: Dirichlet boundary conditions are not allowed
     * \param[in, out] ctx: execution context
     * \param[in] bc: boundary condition
     * \return false
     */
    [[nodiscard]] bool addBoundaryCondition(
        Context& ctx,
        std::unique_ptr<AbstractDirichletBoundaryCondition> bc) noexcept
        override;
    /*!
     * \brief set the evolution of the macroscopic gradients
     * \param[in] ev: function
     */
    virtual void setMacroscopicGradientsEvolution(
        const std::function<std::vector<real>(const real)>& ev);
    /*!
     * \brief return the value of the macroscopic gradients
     * \return the value of the macroscopic gradients at the end of the time
     * step.
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual std::vector<real> getMacroscopicGradients(const real t,
                                                      const real dt) const;
    //! \brief destructor
    ~PeriodicNonLinearEvolutionProblem() override;

   protected:
    //
    [[nodiscard]] bool setup(Context& ctx,
                             const real t,
                             const real dt) noexcept override;
    //! \brief a function describing the evolution of the macroscopic gradients
    std::function<std::vector<real>(const real)>
        macroscopic_gradients_evolution;
  };  // end of PeriodicNonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX */
