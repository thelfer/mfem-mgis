/*!
 * \file   MFEMMGIS/NonLinearModel.hxx
 * \brief  This file declares the `NonLinearModel` class
 * \author Thomas Helfer
 * \date   05/03/2026
 */

#ifndef LIB_MFEM_MGIS_NONLINEARMODEL_HXX
#define LIB_MFEM_MGIS_NONLINEARMODEL_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ModelBase.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  //! \brief a model based on a nonlinear evolution problem
  struct MFEM_MGIS_EXPORT NonLinearModel : ModelBase {
    /*!
     * \brief constructor from a mesh and parameters
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] parameters: parameters of the model and of the non linear
     * evolution problem
     */
    NonLinearModel(Context &ctx,
                   MeshDiscretization &m,
                   const Parameters &parameters);
    /*!
     * \brief constructor from an existing problem
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     */
    NonLinearModel(Context &ctx, std::shared_ptr<NonLinearEvolutionProblem> p);
    //! \return the underlying problem
    [[nodiscard]] NonLinearEvolutionProblem &getProblem() noexcept;
    //! \return the underlying problem
    [[nodiscard]] const NonLinearEvolutionProblem &getProblem() const noexcept;
    //! \return the name of the model, `NonLinearModel` by default
    [[nodiscard]] std::string getName() const noexcept override;
    /*!
     * \brief execute the initial post-processings of the underlying problem
     * \param[in, out] ctx: execution context
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &ctx, const real t) noexcept override;
    /*!
     * \brief call the `ModelBase` implementation
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return true on success
     */
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \brief execute the registered post-processings and, if `b` is true, the
     * post-processings of the underlying problem
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \param[in] b: boolean stating if the time at the end of the time
     * step is a post-processing time
     * \return true on success
     */
    [[nodiscard]] bool executePostProcessingTasks(
        Context &ctx, const TimeStep &ts, const bool b) noexcept override;
    /*!
     * \brief solve the underlying problem over the time step
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return `ExitStatus::recoverableError` if the resolution fails,
     * `ExitStatus::success` and the output of the resolution otherwise
     */
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \brief update the underlying problem
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    /*!
     * \brief revert the underlying problem
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    //! \brief destructor
    ~NonLinearModel() noexcept override;

   private:
    //! \brief underlying non linear evolution problem
    std::shared_ptr<NonLinearEvolutionProblem> problem;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARMODEL_HXX */
