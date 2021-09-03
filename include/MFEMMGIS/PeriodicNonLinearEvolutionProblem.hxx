/*!
 * \file   include/MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX

#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

#ifdef MFEM_USE_MPI

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in] p: problem
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<true>&,
      const mgis::span<const real>&, const mgis::span<const real>&);

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in] p: problem
   * \param[in] bct: impose the zero value on displacement field at xmin, or ymin, or zmin
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<true>&,
      const mfem_mgis::BoundaryConditionType = mfem_mgis::FIX_XMIN);
#endif /* MFEM_USE_MPI */

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in] p: problem
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<false>&,
      const mgis::span<const real>&, const mgis::span<const real>&);

  /*!
   * \brief set the boundary conditions specific to periodic problems
   * \param[in] p: problem
   * \param[in] bct: impose the zero value on displacement field at xmin, or ymin, or zmin
   */
  MFEM_MGIS_EXPORT void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<false>&,
      const mfem_mgis::BoundaryConditionType = mfem_mgis::FIX_XMIN);

  /*!
   * \brief compute minimal distance from corners to point 
   *        identified in vector `nodes` at index `index`
   */
  real getNodesDistance(const mfem::GridFunction& nodes,
			const bool reorder_space,
			const size_t dim,
			const int index,
			const int size,
			const mgis::span<const real>& corner1,
			const mgis::span<const real>& corner2);
  
  
  /*!
   * \brief a base class handling the evolution of the macroscopic gradients
   */
  struct MFEM_MGIS_EXPORT PeriodicNonLinearEvolutionProblem
      : NonLinearEvolutionProblem {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] corner1: one corner of the computation domain
     * \param[in] corner2: second corner of the computation domain
     */
    PeriodicNonLinearEvolutionProblem(
        std::shared_ptr<FiniteElementDiscretization>,
	const mgis::span<const real>&, const mgis::span<const real>&);

    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] bct: impose the zero value on displacement field at xmin, or ymin, or zmin
     */
    PeriodicNonLinearEvolutionProblem(
        std::shared_ptr<FiniteElementDiscretization>,
	const mfem_mgis::BoundaryConditionType = mfem_mgis::FIX_XMIN);
    
    // disable adding boundary conditions
    [[noreturn]] void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) override;
    /*!
     * \brief set the evolution of the macroscopic gradients
     * \param[in] e : function
     */
    virtual void setMacroscopicGradientsEvolution(
        const std::function<std::vector<real>(const real)>&);
    /*!
     * \return the value of the macroscopic gradients at the end of the time
     * step.
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual std::vector<real> getMacroscopicGradients(const real,
                                                      const real) const;
    //! \brief destructor
    virtual ~PeriodicNonLinearEvolutionProblem();

   protected:
    //
    void setup(const real, const real) override;
    //! \brief a function describing the evolution of the macroscopic gradients
    std::function<std::vector<real>(const real)>
        macroscopic_gradients_evolution;
  };  // end of PeriodicNonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM */
