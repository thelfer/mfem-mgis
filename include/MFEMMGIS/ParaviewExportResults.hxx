/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing to export the results to paraview
   */
  template <bool parallel>
  struct ParaviewExportResults final : public PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportResults(NonLinearEvolutionProblemImplementation<parallel>&,
                          const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ParaviewExportResults() override;

   private:
    const int dim, sdim;
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! exported displacement grid function
    mfem_mgis::GridFunction<parallel> displacement;
    //! list of materials
    std::vector<mfem_mgis::Material *> mgis_materials;
    //! nb_material
    int nb_materials;
    //! exported stress grid function
    std::vector<mfem_mgis::GridFunction<parallel>> stress;
    //! exported strain grid function
    std::vector<mfem_mgis::GridFunction<parallel>> strain;
    //!
    size_type cycle;
  };  // end of struct ParaviewExportResults

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportResults.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX */
