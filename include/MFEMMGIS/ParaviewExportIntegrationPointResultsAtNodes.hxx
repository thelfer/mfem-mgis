/*!
 * \file   include/MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx
 * \brief  This file declares the `ParaviewExportIntegrationPointResultsAtNodes`
 * class
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTINTEGRATIONPOINTRESULTSATNODES_HXX

#include <map>
#include <string>
#include <memory>
#include <vector>
#include <variant>
#include <string_view>
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#ifdef MGIS_FUNCTION_SUPPORT
#include "MFEMMGIS/PartialQuadratureFunctionsSet.hxx"
#endif /* MGIS_FUNCTION_SUPPORT */

namespace mfem_mgis {

  // forward declaration
  struct AbstractBehaviourIntegrator;
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
      //! \brief exported functions, one per material
      const std::vector<ImmutablePartialQuadratureFunctionView> functions;
    };
    /*!
     * \brief constructor
     * \param[in] d: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodesBase(const std::string &d);
    //! \brief move constructor
    ParaviewExportIntegrationPointResultsAtNodesBase(
        ParaviewExportIntegrationPointResultsAtNodesBase &&) noexcept = default;
    ParaviewExportIntegrationPointResultsAtNodesBase(
        const ParaviewExportIntegrationPointResultsAtNodesBase &) = delete;

   protected:
    //! \brief description of a result defined at the integration points
    struct MaterialIntegrationPointResultBase {
      //! \brief enumeration of the kind of results that can be post-processed.
      enum Category {
        GRADIENTS,                //!< gradients
        THERMODYNAMIC_FORCES,     //!< thermodynamic forces
        INTERNAL_STATE_VARIABLES  //!< internal state variables
      };
      //! \brief name of the result
      std::string name;
      //! \brief number of components
      size_type number_of_components;
      //! \brief kind of results treated
      Category category;
      /*!
       * \brief behaviour integrators providing the result, in the same order
       * as the material identifiers
       */
      std::vector<const AbstractBehaviourIntegrator *> behaviour_integrators;
    };
    /*!
     * \brief extract the material identifiers from the description of
     * the exported functions
     * \param[in] ds: functions to be exported
     */
    void extractMaterialIdentifiers(
        const std::vector<ExportedFunctionsDescription> &ds);
    /*!
     * \brief get information about the given result (number of components,
     * category and behaviour integrators providing the result)
     *
     * On each material, exactly one behaviour integrator must provide the
     * result, either as a gradient, a thermodynamic force or an internal
     * state variable.
     *
     * \param[in] throwing: dummy attribute to indicate that this function may
     * throw an exception
     * \param[in, out] r: result considered. The name of the result must be
     * set.
     * \param[in] p: non linear evolution problem
     */
    void getResultDescription(
        attributes::Throwing throwing,
        MaterialIntegrationPointResultBase &r,
        const NonLinearEvolutionProblemImplementationBase &p);
    /*!
     * \brief get the functions associated with the given result
     * \return the functions associated with the given result
     * \param[in, out] ctx: execution context
     * \param[in] r: result considered
     * \param[in] s: time step stage
     */
    std::optional<std::vector<ImmutablePartialQuadratureFunctionView>>
    getPartialQuadratureFunctionViews(
        Context &ctx,
        const MaterialIntegrationPointResultBase &r,
        const TimeStepStage s = ets) noexcept;
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief list of material identifiers
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
   * evolution problem
   * 2. explicitly given by the user
   */
  template <bool parallel>
  struct ParaviewExportIntegrationPointResultsAtNodesImplementation final
      : public AbstractNonLinearEvolutionProblemPostProcessing<parallel>,
        public ParaviewExportIntegrationPointResultsAtNodesBase {
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        const Parameters &params);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] d: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        const ExportedFunctionsDescription &d,
        const std::string &n);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] ds: functions to be exported
     * \param[in] n: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        const std::vector<ExportedFunctionsDescription> &ds,
        const std::string &n);
    //! \brief move constructor
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        ParaviewExportIntegrationPointResultsAtNodesImplementation
            &&) noexcept = default;
    ParaviewExportIntegrationPointResultsAtNodesImplementation(
        const ParaviewExportIntegrationPointResultsAtNodesImplementation &) =
        delete;
    /*!
     * \brief execute the post-processing at the initial time
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessing(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        const real t) noexcept override;
    /*!
     * \brief execute the post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool execute(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        const real t,
        const real dt) noexcept override;
    //! \brief destructor
    ~ParaviewExportIntegrationPointResultsAtNodesImplementation() override;

   private:
    //! \brief description of a result and of the grid function exporting it
    struct MaterialIntegrationPointResult
        : public MaterialIntegrationPointResultBase {
      //! \brief grid function
      std::unique_ptr<GridFunction<parallel>> f;
    };
    /*!
     * \brief create the sub mesh once the material identifiers are known
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     */
    void createSubMesh(Context &ctx,
                       NonLinearEvolutionProblemImplementation<parallel> &p);
    /*!
     * \brief update the grid functions and export them
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] t: time
     * \param[in] s: time step stage
     * \return true on success
     */
    bool exportResults(Context &ctx,
                       NonLinearEvolutionProblemImplementation<parallel> &p,
                       const real t,
                       const TimeStepStage s) noexcept;
    //! \brief submesh defined when exporting data
    std::shared_ptr<mfem_mgis::SubMesh<parallel>> submesh;
    //! \brief list of results defined through parameters
    std::vector<MaterialIntegrationPointResult> results;
    /*!
     * \brief a small structure gathering information about fields to be
     * exported
     */
    struct ExportedFunctions {
      //! \brief default constructor
      ExportedFunctions() = default;
      //! \brief move constructor
      ExportedFunctions(ExportedFunctions &&) = default;
      //! \brief name of the exported function in the Paraview's file
      std::string name;
      //! \brief exported functions, one per material
      std::vector<ImmutablePartialQuadratureFunctionView> functions;
      //! \brief exported grid functions corresponding to the exported functions
      std::unique_ptr<GridFunction<parallel>> grid_function;
    };
    //! \brief exported functions
    std::vector<std::unique_ptr<ExportedFunctions>> exported_functions;
    /*!
     * \brief boolean stating if the results shall be exported at the initial
     * time of the simulation
     */
    const bool shallExecuteInitialPostProcessing;
  };  // end of struct
      // ParaviewExportIntegrationPointResultsAtNodesImplementation

  /*!
   * \brief a facade to export quadrature functions in sequential and parallel.
   */
  struct ParaviewExportIntegrationPointResultsAtNodes {
    //! \brief a simple alias
    using ExportedFunctionsDescription =
        ParaviewExportIntegrationPointResultsAtNodesBase::
            ExportedFunctionsDescription;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportIntegrationPointResultsAtNodes(Context &ctx,
                                                 NonLinearEvolutionProblem &p,
                                                 const Parameters &params);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] efcts: functions to be exported
     * \param[in] d: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        Context &ctx,
        NonLinearEvolutionProblem &p,
        const ExportedFunctionsDescription &efcts,
        const std::string &d);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] ds: functions to be exported
     * \param[in] d: output directory name
     */
    ParaviewExportIntegrationPointResultsAtNodes(
        Context &ctx,
        NonLinearEvolutionProblem &p,
        const std::vector<ExportedFunctionsDescription> &ds,
        const std::string &d);
    //! \brief move constructor
    ParaviewExportIntegrationPointResultsAtNodes(
        ParaviewExportIntegrationPointResultsAtNodes &&) = default;
    ParaviewExportIntegrationPointResultsAtNodes(
        const ParaviewExportIntegrationPointResultsAtNodes &) = delete;
    /*!
     * \brief execute the export at the initial time of the simulation
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessing(
        Context &ctx, NonLinearEvolutionProblem &p, const real t) noexcept;
    /*!
     * \brief execute the export
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool execute(Context &ctx,
                               NonLinearEvolutionProblem &p,
                               const real t,
                               const real dt) noexcept;
    //! \brief destructor
    ~ParaviewExportIntegrationPointResultsAtNodes();

   private:
    //! \brief underlying implementation, parallel or sequential
    std::variant<
        std::monostate,
#ifdef MFEM_USE_MPI
        ParaviewExportIntegrationPointResultsAtNodesImplementation<true>,
#endif /* MFEM_USE_MPI */
        ParaviewExportIntegrationPointResultsAtNodesImplementation<false>>
        implementations;
  };

#ifdef MGIS_FUNCTION_SUPPORT
  /*!
   * \brief generates a description of functions to be exported from a set of
   * partial quadrature functions.
   * \param[in] n: name of the exported functions
   * \param[in] f: partial quadrature functions set
   * \return the description of the exported functions
   */
  MFEM_MGIS_EXPORT ParaviewExportIntegrationPointResultsAtNodesBase::
      ExportedFunctionsDescription
      makeExportedFunctionsDescription(std::string_view n,
                                       const PartialQuadratureFunctionsSet &f);
  /*!
   * \brief generates the descriptions of functions to be exported from
   * named sets of partial quadrature functions.
   * \param[in] fcts: sets of partial quadrature functions sorted by name
   * \return the descriptions of the exported functions
   */
  MFEM_MGIS_EXPORT
  std::vector<ParaviewExportIntegrationPointResultsAtNodesBase::
                  ExportedFunctionsDescription>
  makeExportedFunctionsDescriptions(
      const std::map<std::string, const PartialQuadratureFunctionsSet &> &fcts);

  /*!
   * \brief build a set of partial quadrature functions
   * \return a set of partial quadrature functions with the given number of
   * components defined on the given materials
   * \param[in, out] ctx: execution context
   * \param[in] p: non linear evolution problem
   * \param[in] mids: material identifiers
   * \param[in] nc: number of components
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<PartialQuadratureFunctionsSet>
  buildPartialQuadratureFunctionsSet(
      Context &ctx,
      const NonLinearEvolutionProblemImplementationBase &p,
      const std::vector<size_type> &mids,
      const size_type nc) noexcept;

  /*!
   * \brief a post-processing which updates partial quadrature functions and
   * exports them at nodes
   */
  template <bool parallel>
  struct ParaviewExportIntegrationPointPostProcessingsResultsAtNodes
      : public AbstractNonLinearEvolutionProblemPostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] n: name of the field
     * \param[in] mids: material identifiers
     * \param[in] nc: number of components
     * \param[in] fct: function used to compute the exported fields
     * \param[in] d: output directory
     */
    ParaviewExportIntegrationPointPostProcessingsResultsAtNodes(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        std::string_view n,
        const std::vector<size_type> mids,
        const size_type nc,
        std::function<bool(Context &, PartialQuadratureFunction &)> fct,
        std::string_view d);
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessing(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        const real t) noexcept override;
    /*!
     * \brief execute the post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool execute(
        Context &ctx,
        NonLinearEvolutionProblemImplementation<parallel> &p,
        const real t,
        const real dt) noexcept override;

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
