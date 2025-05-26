/*!
 * \file   include/MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX

#include <memory>
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing to export integration points results to paraview
   * after projecting them to nodes
   */
  template <bool parallel>
  struct ParaviewExportIntegrationPointResultsAtNodes final
      : public PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        NonLinearEvolutionProblemImplementation<parallel>&, const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ParaviewExportIntegrationPointResultsAtNodes() override;

   private:
    struct IntegrationPointResult {
      //! \brief enumeration of the kind of results that can be post-processed.
      enum Category {
        GRADIENTS,
        THERMODYNAMIC_FORCES,
        INTERNAL_STATE_VARIABLES
      };
      //! \brief name of the result
      std::string name;
      //! \brief number of components
      size_type number_of_components;
      //! \brief kind of results treated
      Category category;
      //! \brief finite element space
      std::unique_ptr<FiniteElementSpace<parallel>> fespace;
      //! \brief grid function
      std::unique_ptr<GridFunction<parallel>> f;
    };
    /*!
     * \brief get information about the given result (number of components and
     * category)
     */
    void getResultDescription(
        IntegrationPointResult&,
        const NonLinearEvolutionProblemImplementation<parallel>&);
    /*!
     * \return the functions associated with the given result
     * \param[in] p: non linear evolution problem
     * \param[in] r: results considered
     */
    std::vector<ImmutablePartialQuadratureFunctionView>
    getPartialQuadratureFunctionViews(
        const NonLinearEvolutionProblemImplementation<parallel>&,
        const IntegrationPointResult&);
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief list of material' identifiers
    std::vector<size_type> materials_identifiers;
    //! \brief list of results
    std::vector<IntegrationPointResult> results;
    //! \brief number of records
    size_type cycle;
  };  // end of struct ParaviewExportIntegrationPointResultsAtNodes

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX */
