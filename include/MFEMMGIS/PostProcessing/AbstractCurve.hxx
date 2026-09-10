/*!
 * \file   MFEMMGIS/PostProcessing/AbstractCurve.hxx
 * \brief  This file declares the `AbstractCurve` class.
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_ABSTRACTCURVE_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_ABSTRACTCURVE_HXX

#include <string>
#include <vector>
#include <optional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;

  /*!
   * \brief base struct for curves
   *
   * A curve is defined as an abstraction to get compute a set of values
   * for the post-processing of a simulation, such as:
   *
   * - the value of a nodal field at a given point
   * - the integral value of a field on integration points
   * - the mean value of a field on integration points
   * - the value of an uniform evaluator
   * - etc...
   *
   * Most curves imply a reduction across MPI processes.
   */
  struct MFEM_MGIS_EXPORT AbstractCurve {
    //! \return a description of each value returned by the curve
    [[nodiscard]] virtual std::vector<std::string> getDescriptions()
        const noexcept = 0;
    /*!
     * \brief return the values of the curve at the given time
     * \param[in] ctx: execution context
     * \param[in] ts: time step stage
     *
     * \note all MPI processes are synchronized after this call
     */
    [[nodiscard]] virtual std::optional<std::vector<real>> getValues(
        Context &, const TimeStepStage) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractCurve() noexcept;
  };  // end of AbstractCurve

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_ABSTRACTCURVE_HXX */
