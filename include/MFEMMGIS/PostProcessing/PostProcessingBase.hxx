/*!
 * \file   MFEMMGIS/PostProcessing/PostProcessingBase.hxx
 * \brief  This file declares the `PostProcessingBase` class
 * \author Thomas Helfer
 * \date   29/09/2023
 */

#ifndef LIB_MFEMMGIS_POSTPROCESSING_POSTPROCESSINGBASE_HXX
#define LIB_MFEMMGIS_POSTPROCESSING_POSTPROCESSINGBASE_HXX

#include <string>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PostProcessing/AbstractPostProcessing.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  //! \brief a base class suitable for most post-processings
  struct MFEM_MGIS_EXPORT PostProcessingBase : public AbstractPostProcessing {
    //! \return a description of each parameters of this class
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] ps: physical system
     * \param[in] parameters: parameters
     * \param[in] defaultAllTimeStepsValue: default value for the `allTimeSteps`
     * parameter
     */
    PostProcessingBase(PhysicalSystem &, const Parameters &, const bool);
    //
    [[nodiscard]] PhysicalSystem &getPhysicalSystem() noexcept override;
    [[nodiscard]] const PhysicalSystem &getPhysicalSystem()
        const noexcept override;
    //! \brief destructor
    ~PostProcessingBase() noexcept override;

   protected:
    //! \brief the underlying physical system
    PhysicalSystem &physicalSystem;
    /*!
     * \brief boolean stating if the post-processing must be executed at the end
     * of all time steps or only at time explicitely marked by the user as a
     * *post-processing time*
     */
    const bool allTimeSteps;
  };  // end of PostProcessingBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POSTPROCESSING_POSTPROCESSINGBASE_HXX */
