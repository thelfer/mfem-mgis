/*!
 * \file   include/MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX

#include <map>
#include <string>
#include <memory>
#include <variant>
#include <string_view>
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#ifdef MGIS_FUNCTION_SUPPORT
#include "MFEMMGIS/PartialQuadratureFunctionsSet.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

namespace mfem_mgis {

#ifdef MGIS_FUNCTION_SUPPORT
  // forward declaration
  struct PartialQuadratureFunctionsSet;
#endif /* MGIS_FUNCTION_SUPPORT */

  /*!
   * \brief a base class to factorize methods between sequential and parallel
   * versions
   */
  struct MFEM_MGIS_EXPORT ParaviewExportIntegrationPointResultsAtNodesBase {
    //! \brief a simple structure to describe the functions to be exported.
    struct ExportedFunctionsDescription {
      //! \brief name of the exported function in the Paraview's file
      const std::string name;
      //! \brief list of exported function name
      const std::vector<ImmutablePartialQuadratureFunctionView> functions;
    };
    /*!
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodesBase(const std::string &);

   protected:
    //
    struct MaterialIntegrationPointResultBase {
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
    };
    /*!
     * \brief extract the material identifiers from the description of
     * of the exported functions
     * \param[in] ds: functions to be exported
     */
    void extractMaterialIdentifiers(
        const std::vector<ExportedFunctionsDescription> &);
    /*!
     * \brief get information about the given result (number of components and
     * category)
     */
    void getResultDescription(
        MaterialIntegrationPointResultBase &,
        const NonLinearEvolutionProblemImplementationBase &);
    /*!
     * \return the functions associated with the given result
     * \param[in] p: non linear evolution problem
     * \param[in] r: results considered
     */
    std::vector<ImmutablePartialQuadratureFunctionView>
    getPartialQuadratureFunctionViews(
        const NonLinearEvolutionProblemImplementationBase &,
        const MaterialIntegrationPointResultBase &);
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief list of material' identifiers
    std::vector<size_type> materials_identifiers;
    //! \brief number of records
    size_type cycle;
  };

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
  struct ParaviewExportIntegrationPointResultsAtNodesImplementation final
      : public PostProcessing<parallel>,
        public ParaviewExportIntegrationPointResultsAtNodesBase {
    /*!
     * \brief constructor
     * \param[in,out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        Context &,
        NonLinearEvolutionProblemImplementation<parallel> &,
        const Parameters &);
    /*!
     * \brief constructor
     * \param[in,out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] d: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        Context &,
        NonLinearEvolutionProblemImplementation<parallel> &,
        const ExportedFunctionsDescription &,
        const std::string &);
    /*!
     * \brief constructor
     * \param[in,out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] ds: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        Context &,
        NonLinearEvolutionProblemImplementation<parallel> &,
        const std::vector<ExportedFunctionsDescription> &,
        const std::string &);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel> &,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ParaviewExportIntegrationPointResultsAtNodesImplementation() override;

   private:
    struct MaterialIntegrationPointResult
        : public MaterialIntegrationPointResultBase {
      //! \brief finite element space
      std::unique_ptr<FiniteElementSpace<parallel>> fespace;
      //! \brief grid function
      std::unique_ptr<GridFunction<parallel>> f;
    };
    /*!
     * \brief create the sub mesh once the material identifiers are known
     * \param[in] p: non linear problem
     */
    void createSubMesh(NonLinearEvolutionProblemImplementation<parallel> &);
    //! \brief submesh defined when exporting data
    std::shared_ptr<mfem_mgis::SubMesh<parallel>> submesh;
    //! \brief list of results defined through parameters
    std::vector<MaterialIntegrationPointResult> results;
    /*!
     * \brief a small structure gathering information about fields to be
     * exported
     */
    struct ExportedFunctions {
      ExportedFunctions() = default;
      ExportedFunctions(ExportedFunctions &&) = default;
      //! \brief name of the exported function in the Paraview's file
      std::string name;
      //! \brief list of exported function name
      std::vector<ImmutablePartialQuadratureFunctionView> functions;
      //! \brief finite element space used to define the exported grid function
      std::unique_ptr<FiniteElementSpace<parallel>> grid_function_fespace;
      //! \brief exported grid functions corresponding to the exported functions
      std::unique_ptr<GridFunction<parallel>> grid_function;
    };
    std::vector<std::unique_ptr<ExportedFunctions>> exported_functions;
  };  // end of struct
      // ParaviewExportIntegrationPointResultsAtNodesImplementation

  /*!
   * \brief a facade to export quadrature function in sequential and parallel.
   */
  struct ParaviewExportIntegrationPointResultsAtNodes {
    //
    using ExportedFunctionsDescription =
        ParaviewExportIntegrationPointResultsAtNodesBase::
            ExportedFunctionsDescription;
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportIntegrationPointResultsAtNodes(NonLinearEvolutionProblem &,
                                                 const Parameters &);
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] d: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        NonLinearEvolutionProblem &,
        const ExportedFunctionsDescription &,
        const std::string &);
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] ds: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        NonLinearEvolutionProblem &,
        const std::vector<ExportedFunctionsDescription> &,
        const std::string &);
    //
    void execute(NonLinearEvolutionProblem &, const real, const real);
    //! \brief destructor
    ~ParaviewExportIntegrationPointResultsAtNodes();

   private:
    std::variant<
        std::monostate,
        ParaviewExportIntegrationPointResultsAtNodesImplementation<true>,
        ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>
        implementations;
  };

#ifdef MGIS_FUNCTION_SUPPORT
  /*!
   * \brief generates a description of functions to be exported from  a set of
   * partial quadrature functions.
   * \param[in] n: name of the exported functions
   * \param[in] f: partial quadrature functions setx
   */
  MFEM_MGIS_EXPORT ParaviewExportIntegrationPointResultsAtNodesBase::
      ExportedFunctionsDescription
      makeExportedFunctionsDescription(std::string_view,
                                       const PartialQuadratureFunctionsSet &);
  /*!
   *
   */
  MFEM_MGIS_EXPORT
  std::vector<ParaviewExportIntegrationPointResultsAtNodesBase::
                  ExportedFunctionsDescription>
  makeExportedFunctionsDescriptions(
      const std::map<std::string, const PartialQuadratureFunctionsSet &> &);

  MFEM_MGIS_EXPORT PartialQuadratureFunctionsSet
  buildPartialQuadratureFunctionsSet(
      const NonLinearEvolutionProblemImplementationBase &,
      const std::vector<size_type> &,
      const size_type);

  template <bool parallel>
  struct ParaviewExportIntegrationPointPostProcessingsResultsAtNodes
      : public PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] n: name of the field
     * \param[in] mids: material identifier
     * \param[in] nc: number of components
     * \param[in] fct: function used to compute the exported fields
     * \param[in] d: output directory
     */
    ParaviewExportIntegrationPointPostProcessingsResultsAtNodes(
        NonLinearEvolutionProblemImplementation<true> &,
        std::string_view,
        const std::vector<size_type>,
        const size_type,
        std::function<bool(Context &, PartialQuadratureFunction &)>,
        std::string_view);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel> &p,
                 const real t,
                 const real dt) override;

   private:
    //! \brief exported functions
    PartialQuadratureFunctionsSet functions;
    //! \brief update function
    std::function<bool(Context &, PartialQuadratureFunction &)> update_function;
    //! \brief exporter
    ParaviewExportIntegrationPointResultsAtNodesImplementation<parallel>
        exporter;
  };

#endif /* MGIS_FUNCTION_SUPPORT */

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX */
