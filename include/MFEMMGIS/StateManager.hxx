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
    StateManager(const MeshDiscretization& m) noexcept;
    /*!
     * \brief register a partial quadrature function
     *
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the partial quadrature function
     * \param[in] f: partial quadrature function
     * \param[in] ts: time step stage
     * \return true on success
     *
     * \note the caller is responsible for keeping the viewed function alive
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           ImmutablePartialQuadratureFunctionView f,
                           const TimeStepStage ts) noexcept;
    /*!
     * \brief register a partial quadrature function
     *
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the partial quadrature function
     * \param[in] f: partial quadrature function
     * \param[in] ts: time step stage
     * \return true on success
     *
     * \note the shared pointer is stored internally
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           const std::shared_ptr<PartialQuadratureFunction>& f,
                           const TimeStepStage ts) noexcept;
    /*!
     * \return if a partial quadrature function with the given name has been
     * registered
     *
     * \param[in, out] ctx: execution context
     * \param[in] qspace: partial quadrature space
     * \param[in] n: function name
     * \param[in] ts: time step stage
     */
    [[nodiscard]] std::optional<bool> contains(
        Context& ctx,
        const std::shared_ptr<const PartialQuadratureSpace>& qspace,
        std::string_view n,
        const TimeStepStage ts) const noexcept;
    /*!
     * \return the registered partial quadrature function
     *
     * \param[in, out] ctx: execution context
     * \param[in] qspace: partial quadrature space
     * \param[in] n: function name
     * \param[in] ts: time step stage
     */
    [[nodiscard]] std::optional<ImmutablePartialQuadratureFunctionView> get(
        Context& ctx,
        const std::shared_ptr<const PartialQuadratureSpace>& qspace,
        std::string_view n,
        const TimeStepStage ts) const noexcept;

    //! \brief destructor
    ~StateManager() noexcept;

   private:
    //! \brief partial quadrature functions registered on a quadrature space
    struct PartialQuadratureFunctionManager {
      /*!
       * \brief register a partial quadrature function
       *
       * \param[in, out] ctx: execution context
       * \param[in] n: partial quadrature function name
       * \param[in] f: partial quadrature function
       * \param[in] ts: time step stage
       * \return true on success
       */
      [[nodiscard]] bool add(Context& ctx,
                             std::string_view n,
                             ImmutablePartialQuadratureFunctionView f,
                             const TimeStepStage ts) noexcept;
      /*!
       * \return if a partial quadrature function with the given name has been
       * registered
       *
       * \param[in] n: function name
       * \param[in] ts: time step stage
       */
      [[nodiscard]] bool contains(std::string_view n,
                                  const TimeStepStage ts) const noexcept;
      /*!
       * \return the registered partial quadrature function
       *
       * \param[in, out] ctx: execution context
       * \param[in] n: function name
       * \param[in] ts: time step stage
       */
      [[nodiscard]] std::optional<ImmutablePartialQuadratureFunctionView> get(
          Context& ctx,
          std::string_view n,
          const TimeStepStage ts) const noexcept;

     protected:
      /*!
       * \brief registered partial quadrature functions at the beginning of the
       * time step
       */
      std::map<std::string, ImmutablePartialQuadratureFunctionView, std::less<>>
          qfunctions_bts;
      /*!
       * \brief registered partial quadrature functions at the end of the
       * time step
       */
      std::map<std::string, ImmutablePartialQuadratureFunctionView, std::less<>>
          qfunctions_ets;
    };
    /*!
     * \brief list of partial quadrature function managers, sorted by quadrature
     * space identifiers
     */
    std::map<std::pair<LocationIdentifier, size_type>,
             std::unique_ptr<PartialQuadratureFunctionManager>>
        qfunctions;
    //! \brief functions kept alive by the state manager
    std::vector<std::shared_ptr<PartialQuadratureFunction>> qfunctions_pointers;
  };  // end of StateManager

  /*!
   * \brief declare all the partial quadrature functions defined by the
   * behaviour integrators of a nonlinear evolution problem in the given state
   * manager
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] s: state manager
   * \param[in] p: nonlinear evolution problem
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool addPartialQuadratureFunctions(
      Context& ctx,
      StateManager& s,
      const AbstractNonLinearEvolutionProblem& p) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_STATEMANAGER_HXX */
