/*!
 * \file   MFEMMGIS/LoopCouplingScheme.hxx
 * \brief  This file declares the `LoopCouplingScheme` class
 * \date   05/12/2022
 */

#ifndef LIB_MFEM_MGIS_LOOP_COUPLING_SCHEME_HXX
#define LIB_MFEM_MGIS_LOOP_COUPLING_SCHEME_HXX

#include <map>
#include <string>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/CouplingSchemeBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief the simplest coupling scheme: all coupling items are called a fixed
   * number of times
   */
  struct MFEM_MGIS_EXPORT LoopCouplingScheme : CouplingSchemeBase {
    //! \return a description of this scheme
    static std::string getDescription() noexcept;
    //! \return a description of the parameters of this scheme
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    LoopCouplingScheme(Context &ctx, const MeshDiscretization &m);
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] params: parameters
     */
    LoopCouplingScheme(Context &ctx,
                       const MeshDiscretization &m,
                       const Parameters &params);
    /*!
     * \brief set the number of iterations
     * \param[in, out] ctx: execution context
     * \param[in] n: number of iterations
     * \return true on success
     */
    [[nodiscard]] bool setNumberOfIterations(Context &ctx,
                                             const size_type n) noexcept;
    //! \return the name of the scheme, `LoopCouplingScheme` by default
    [[nodiscard]] std::string getName() const noexcept override;
    /*!
     * \return a description of the coupling scheme
     * \param[in, out] ctx: execution context
     * \param[in] b: boolean being the default value for information requests
     * \param[in] parameters: information requests. Supported requests are
     * `ShortDescription`, `NumericalParameters` and `CouplingItems`.
     */
    [[nodiscard]] std::optional<std::string> describe(
        Context &ctx,
        const bool b,
        const Parameters &parameters) const noexcept override;
    //     [[nodiscard]] bool addConvergenceCriterion(
    //         Context &,
    //         std::string_view,
    //         const Parameters &) noexcept override final;
    /*!
     * \brief report an error: this scheme does not support convergence criteria
     * \param[in, out] ctx: execution context
     * \param[in] c: convergence criterion
     * \return false
     */
    [[nodiscard]] bool addConvergenceCriterion(
        Context &ctx,
        std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion> c) noexcept
        override final;
    /*!
     * \brief call all coupling items a fixed number of times
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return the exit status and, on success, the outputs of the scheme
     */
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &ctx, const TimeStep &ts) noexcept override;
    //! \brief destructor
    ~LoopCouplingScheme() noexcept override;

   private:
    //! \brief number of iterations
    size_type number_of_iterations = 1;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LOOP_COUPLING_SCHEME_HXX */
