/*!
 * \file   MFEMMGIS/StateManager.hxx
 * \brief  This file declares the `StateManager` class.
 * \author Thomas Helfer
 * \date   29/03/2026
 */

#ifndef LIB_MFEMMGIS_STATEMANAGER_HXX
#define LIB_MFEMMGIS_STATEMANAGER_HXX

#include <map>
#include <string>
#include <memory>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/PartialQuadratureSpaceIdentifiersManager.hxx"

namespace mfem_mgis {

  // forward declaration
  struct AbstractNonLinearEvolutionProblem;

  /*!
   * \brief structure containing the state of a physical system or a set of non
   * linear evolution problems.
   */
  struct MFEM_MGIS_EXPORT StateManager
      : PartialQuadratureSpaceIdentifiersManager {
    /*!
     * \brief constructor from a mesh discretization
     * \param[in] m: mesh discretization
     */
    StateManager(const MeshDiscretization&) noexcept;
    /*!
     * \brief register a partial quadrature function
     *
     * \param[in] ctx: execution context
     * \param[in] n: name of the partial quadrature function
     * \param[in] f: partial quadrature function
     * \param[in] ts: time step stage
     *
     * \note the callee is responsible for keeping the view alive
     */
    [[nodiscard]] bool add(Context&,
                           std::string_view,
                           ImmutablePartialQuadratureFunctionView,
                           const TimeStepStage) noexcept;
    /*!
     * \brief register a partial quadrature function
     *
     * \param[in] ctx: execution context
     * \param[in] n: name of the partial quadrature function
     * \param[in] f: partial quadrature function
     * \param[in] ts: time step stage
     *
     * \note the shared pointer is stored internally
     */
    [[nodiscard]] bool add(Context&,
                           std::string_view,
                           const std::shared_ptr<PartialQuadratureFunction>&,
                           const TimeStepStage) noexcept;
    /*!
     * \return if a partial quadrature function with the given name has been
     * registred
     *
     * \param[in] ctx: execution context
     * \param[in] s: partial quadrature space
     * \param[in] f: function name
     * \param[in] ts: time step stage
     */
    [[nodiscard]] std::optional<bool> contains(
        Context&,
        const std::shared_ptr<const PartialQuadratureSpace>&,
        std::string_view,
        const TimeStepStage) const noexcept;
    /*!
     * \return the registred partial quadrature function
     *
     * \param[in] ctx: execution context
     * \param[in] s: partial quadrature space
     * \param[in] f: function name
     * \param[in] ts: time step stage
     */
    [[nodiscard]] std::optional<ImmutablePartialQuadratureFunctionView> get(
        Context&,
        const std::shared_ptr<const PartialQuadratureSpace>&,
        std::string_view,
        const TimeStepStage) const noexcept;

    //! \brief destructor
    ~StateManager() noexcept;

   private:
    //
    struct PartialQuadratureFunctionManager {
      /*!
       * \return the registred partial quadrature function
       *
       * \param[in] ctx: execution context
       * \param[in] n: partial quadrature function name
       * \param[in] f: partial quadrature function
       * \param[in] ts: time step stage
       */
      [[nodiscard]] bool add(Context&,
                             std::string_view,
                             ImmutablePartialQuadratureFunctionView,
                             const TimeStepStage) noexcept;
      /*!
       * \return if a partial quadrature function with the given name has been
       * registred
       *
       * \param[in] f: function name
       * \param[in] ts: time step stage
       */
      [[nodiscard]] bool contains(std::string_view,
                                  const TimeStepStage) const noexcept;
      /*!
       * \return the registred partial quadrature function
       *
       * \param[in] ctx: execution context
       * \param[in] f: function name
       * \param[in] ts: time step stage
       */
      [[nodiscard]] std::optional<ImmutablePartialQuadratureFunctionView> get(
          Context&, std::string_view, const TimeStepStage) const noexcept;

     protected:
      /*!
       * \brief registred partial quadrature functions at the beginning of the
       * time step
       */
      std::map<std::string, ImmutablePartialQuadratureFunctionView, std::less<>>
          qfunctions_bts;
      /*!
       * \brief registred partial quadrature functions at the end of the
       * time step
       */
      std::map<std::string, ImmutablePartialQuadratureFunctionView, std::less<>>
          qfunctions_ets;
    };
    /*!
     * \brief list of partial quadrature function managers, sorted by quadrature
     * space identifiers
     */
    std::map<std::pair<size_type, size_type>,
             std::unique_ptr<PartialQuadratureFunctionManager>>
        qfunctions;
    //
    std::vector<std::shared_ptr<PartialQuadratureFunction>> qfunctions_pointers;
  };  // end of StateManager

  /*!
   * \brief declare all the partial quadrature functions defined by the
   * behaviour integrators of a nonlinear evolution problem in the given state
   * manager
   *
   * \param[in] ctx: execution context
   * \param[in] sm: state manager
   * \param[in] p: nonlinear evolution problem
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addPartialQuadratureFunctions(
      Context&,
      StateManager&,
      const AbstractNonLinearEvolutionProblem&) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_STATEMANAGER_HXX */
