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
   *
   * The functions to be exported can be either:
   *
   * 1. automatically extracted from the materials defined in a nonlinear
   *    evolution problem
   * 2. explicitely given by the user
   */
  template <bool parallel>
  struct ParaviewExportIntegrationPointResultsAtNodes final
      : public PostProcessing<parallel> {
    //! \brief a simple structure to describe the functions to be exported.
    struct ExportedFunctionsDescription {
      //! \brief name of the exported function in the Paraview's file
      const std::string name;
      //! \brief list of exported function name
      const std::vector<ImmutablePartialQuadratureFunctionView> functions;
    };
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        NonLinearEvolutionProblemImplementation<parallel>&, const Parameters&);
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] d: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        NonLinearEvolutionProblemImplementation<parallel> &,
        const ExportedFunctionsDescription &,
        const std::string &);
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] d: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        NonLinearEvolutionProblemImplementation<parallel> &,
        const std::vector<ExportedFunctionsDescription> &,
        const std::string &);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ParaviewExportIntegrationPointResultsAtNodes() override;

   private:
    struct MaterialIntegrationPointResult {
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
        MaterialIntegrationPointResult&,
        const NonLinearEvolutionProblemImplementation<parallel>&);
    /*!
     * \return the functions associated with the given result
     * \param[in] p: non linear evolution problem
     * \param[in] r: results considered
     */
    std::vector<ImmutablePartialQuadratureFunctionView>
    getPartialQuadratureFunctionViews(
        const NonLinearEvolutionProblemImplementation<parallel>&,
        const MaterialIntegrationPointResult&);
    //! \brief submesh defined when exporting data for domain or boundary
    //! attributes
    std::shared_ptr<mfem_mgis::SubMesh<parallel>> submesh;
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief list of material' identifiers
    std::vector<size_type> materials_identifiers;
    //! \brief list of results defined through parameters
    std::vector<MaterialIntegrationPointResult> results;
    /*!
     * \brief a small structure gathering information about fields to be exported
     */
    struct ExportedFunctions {
      ExportedFunctions() = default;
      ExportedFunctions(ExportedFunctions &&) = default;
      //! \brief name of the exported function in the Paraview's file
      std::string name;
      //! \brief list of exported function name
      std::vector<ImmutablePartialQuadratureFunctionView> functions;
      //! \brief finite element space used to define the exported grid function
      std::unique_ptr<FiniteElementSpace<parallel>>
          grid_function_fespace;
      //! \brief exported grid functions corresponding to the exported functions
      std::unique_ptr<GridFunction<parallel>> grid_function;
    };
    std::vector<std::unique_ptr<ExportedFunctions>> exported_functions;
    //! \brief number of records
    size_type cycle;
  };  // end of struct ParaviewExportIntegrationPointResultsAtNodes

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX */
