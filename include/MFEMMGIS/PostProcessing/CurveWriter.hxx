/*!
 * \file   MFEMMGIS/PostProcessing/CurveWriter.hxx
 * \brief  This file declares the `CurveWriter` classs
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_CURVEWRITER_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_CURVEWRITER_HXX

#include <fstream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/Utilities/DataFileUtilities.hxx"
#include "MFEMMGIS/PostProcessing/MultipleCurves.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;
  struct MeshDiscretization;
  struct PhysicalSystem;

  /*!
   * \brief helper class meant to write the results of a curve to
   * an output file
   */
  struct MFEM_MGIS_EXPORT CurveWriter {
    //! \return a description of each parameters
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] ctx: execution context
     * \param[in] m: mesh discretization
     * \param[in] params: parameters
     */
    CurveWriter(Context &, const MeshDiscretization &, const Parameters &);
    /*!
     * \brief constructor
     * \param[in] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] params: parameters
     */
    CurveWriter(Context &, const PhysicalSystem &, const Parameters &);
    /*!
     * \brief add a new curve
     * \param[in] ctx: execution context
     * \param[in] c: curve
     */
    [[nodiscard]] bool addCurve(Context &,
                                std::shared_ptr<const AbstractCurve>);
    /*!
     * \brief write the file headers
     * \param[in] ctx: execution context
     */
    bool writeFileHeader(Context &);
    /*!
     * \brief write the file headers
     * \param[in] ctx: execution context
     * \param[in] t: current time
     * \param[in] ts: time step stage
     */
    bool writeValues(Context &, const real, const TimeStepStage &);

   private:
    //! \brief list of registred curves packed into a `MultipleCurves`
    MultipleCurves curves;
    //! \brief output file
    std::ofstream out;
    //! \brief data file format
    DataFileFormat fileFormat = DataFileFormat::TXT;
    //! \brief boolean stating if new curves are allowed
    bool allowNewCuves = true;
    /*!
     * \brief  boolean stating if the current MPI process is the main of the
     * group
     */
    const bool isMainProcess;
  };  // end of CurveWriter

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_CURVEWRITER_HXX */
