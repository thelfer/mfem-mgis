/*!
 * \file   MFEMMGIS/PostProcessing/PointsSetCurvesWriter.hxx
 * \brief  This file declares the `PointsSetCurvesWriter` classs
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESWRITER_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESWRITER_HXX

#include <fstream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStep.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/Utilities/DataFileUtilities.hxx"
#include "MFEMMGIS/PostProcessing/PointsSetCurves.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;
  struct PhysicalSystem;

  /*!
   * \brief helper class meant to write the results of a curve to
   * an output file
   */
  struct MFEM_MGIS_EXPORT PointsSetCurvesWriter {
    //! \return a description of each parameters
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] m: mesh discretization
     * \param[in] params: parameters
     */
    PointsSetCurvesWriter(const MeshDiscretization &, const Parameters &);
    /*!
     * \brief constructor
     * \param[in] ps: physical system
     * \param[in] params: parameters
     */
    PointsSetCurvesWriter(const PhysicalSystem &, const Parameters &);
#ifdef MFEM_USE_MPI
    /*!
     * \brief add a grid function  (parallel version)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     */
    [[nodiscard]] bool add(Context &,
                           std::string_view,
                           const GridFunction<true> &) noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add a grid function (sequential version)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the grid function
     * \param[in] f: grid function
     */
    [[nodiscard]] bool add(Context &,
                           std::string_view,
                           const GridFunction<false> &) noexcept;
    /*!
     * \brief write the file headers
     * \param[in] ctx: execution context
     */
    bool writeFileHeader(Context &);
    /*!
     * \brief write the file headers
     * \param[in] ctx: execution context
     * \param[in] ts: time step
     * \param[in] tss: time step stage
     */
    bool writeValues(Context &, const TimeStep &, const TimeStepStage &);

   private:
    //! \brief mesh
    MeshDiscretization mesh;
    //! \brief list of registred curves packed into a `MultiplePointsSetCurvess`
    PointsSetCurves curves;
    //! \brief output file
    std::ofstream out;
    //! \brief boolean stating if new curves are allowed
    bool allowNewPointsSetCurves = true;
  };  // end of PointsSetCurvesWriter

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_POINTSSETCURVESWRITER_HXX */
