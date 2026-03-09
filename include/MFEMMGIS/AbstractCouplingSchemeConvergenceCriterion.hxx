/*!
 * \file   MFEMMGIS/AbstractCouplingSchemeConvergenceCriterion.hxx
 * \brief  This file declares the `AbstractCouplingSchemeConvergenceCriterion`
 * class
 * \date   05/12/2022
 */

#ifndef LIB_MFEM_MGIS_ABSTRACT_COUPLING_SCHEME_CONVERGENCE_CRITERION_HXX
#define LIB_MFEM_MGIS_ABSTRACT_COUPLING_SCHEME_CONVERGENCE_CRITERION_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractCouplingItem.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  struct ComputeNextStateOutput;

  //! \brief interface of all coupling schemes
  struct MFEM_MGIS_EXPORT AbstractCouplingSchemeConvergenceCriterion {
    /*!
     * \brief method called at the beginning of a time step
     *
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     *
     * \note if required, the time step can be retrieved from the clock held by
     * the physical system
     */
    virtual bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept = 0;
    /*!
     * \return if the criterion is satisfied
     * \param[in] ctx: execution context
     * \param[in] o: output of all items of the
     * coupling scheme
     */
    virtual std::optional<bool> check(
        Context &, const ComputeNextStateOutput &) const noexcept = 0;
    /*!
     * \brief update the state of the system for the next time step
     *
     * This method is called at the end of the time step once the state of
     * the system at the end of the time step is known. This state is generally
     * defined by a set of unknown fields (usually defined as node fields or
     * element fields) and a set of internal state variables (usually defined as
     * fields at integration points).
     *
     * This method is first meant to copy the fields at the end of the time on
     * the fields at the beginning of the time step.
     */
    virtual bool update(Context &) noexcept = 0;
    /*!
     * \brief revert the state of the system at the beginning of the time step
     *
     * This method is typically called in case of non convergence of the
     * coupling scheme to restart the computation with a smaller the time step.
     *
     * This method is first meant to copy the fields at the end of the time on
     * the fields at the beginning of the time step.
     */
    virtual bool revert(Context &) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractCouplingSchemeConvergenceCriterion() noexcept;
  };  // end of AbstractCouplingSchemeConvergenceCriterion

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACT_COUPLING_SCHEME_CONVERGENCE_CRITERION_HXX */
