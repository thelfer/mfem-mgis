/*!
 * \file   MFEMMGIS/PostProcessing/CurvesWriter.hxx
 * \brief  This file declares the `CurvesWriter` classs
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_CURVESWRITER_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_CURVESWRITER_HXX

#include <fstream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/Utilities/DataFileUtilities.hxx"
#include "MFEMMGIS/PostProcessing/MultipleCurves.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;
  struct MeshDiscretization;
  struct PhysicalSystem;

  /*!
   * \brief helper class meant to write the results of curves to
   * an output file
   */
  struct MFEM_MGIS_EXPORT CurvesWriter {
    //! \return a description of each parameter
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh discretization
     * \param[in] parameters: parameters
     */
    CurvesWriter(Context& ctx,
                 const MeshDiscretization& m,
                 const Parameters& parameters);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] parameters: parameters
     */
    CurvesWriter(Context& ctx,
                 const PhysicalSystem& ps,
                 const Parameters& parameters);
    /*!
     * \brief add a new curve
     * \param[in, out] ctx: execution context
     * \param[in] c: curve
     * \return true on success
     */
    [[nodiscard]] bool addCurve(Context& ctx,
                                std::shared_ptr<const AbstractCurve> c);
    /*!
     * \brief write the file headers
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    bool writeFileHeader(Context& ctx);
    /*!
     * \brief write the values of the curves
     * \param[in, out] ctx: execution context
     * \param[in] ts: time step
     * \param[in] tss: time step stage
     * \return true on success
     */
    bool writeValues(Context& ctx,
                     const TimeStep& ts,
                     const TimeStepStage& tss);

   private:
    //! \brief list of registered curves packed into a `MultipleCurves`
    MultipleCurves curves;
    //! \brief output file
    std::ofstream out;
    //! \brief data file format
    DataFileFormat fileFormat = DataFileFormat::TXT;
    //! \brief boolean stating if new curves are allowed
    bool allowNewCurves = true;
    /*!
     * \brief  boolean stating if the current MPI process is the main of the
     * group
     */
    const bool isMainProcess;
  };  // end of CurvesWriter

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_CURVESWRITER_HXX */
