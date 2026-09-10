/*!
 * \file   MFEMMGIS/PostProcessing/MultipleCurves.hxx
 * \brief  This file declares the `MultipleCurves` class
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#include <memory>
#include <vector>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PostProcessing/AbstractCurve.hxx"

#ifndef LIB_MFEMMGIS_POSTPROCESSING_CURVE_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_CURVE_HXX

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  struct PhysicalSystem;

  //! \brief class gathering the results of various curves
  struct MFEM_MGIS_EXPORT MultipleCurves : AbstractCurve {
    //! \brief default constructor
    MultipleCurves();
    /*!
     * \brief constructor
     * \param[in] ctx: execution context
     * \param[in] ps: physical system
     * \param[in] parameters: parameters
     */
    MultipleCurves(Context &, PhysicalSystem &, const Parameters &);
    /*!
     * \brief add a new curve
     * \param[in] ctx: execution context
     * \param[in] c: curve
     */
    [[nodiscard]] bool addCurve(Context &,
                                std::shared_ptr<const AbstractCurve>);
    //
    [[nodiscard]] std::vector<std::string> getDescriptions()
        const noexcept override;
    [[nodiscard]] std::optional<std::vector<real>> getValues(
        Context &ctx, const TimeStepStage) const noexcept override;
    //! \brief destructor
    ~MultipleCurves() noexcept override;

   private:
    //! \brief underlying curves
    std::vector<std::shared_ptr<const AbstractCurve>> curves;
  };  // end of MultipleCurves

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_CURVE_HXX */
