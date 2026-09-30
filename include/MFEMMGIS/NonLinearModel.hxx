/*!
 * \file   MFEMMGIS/NonLinearModel.hxx
 * \brief
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
    NonLinearModel(Context &ctx,
                   MeshDiscretization &m,
                   const Parameters &parameters);
    NonLinearModel(Context &ctx, std::shared_ptr<NonLinearEvolutionProblem> p);
    //
    [[nodiscard]] NonLinearEvolutionProblem &getProblem() noexcept;
    [[nodiscard]] const NonLinearEvolutionProblem &getProblem() const noexcept;
    //
    [[nodiscard]] std::string getName() const noexcept override;
    [[nodiscard]] bool executeInitialPostProcessingTasks(
        Context &ctx, const real t) noexcept override;
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &ctx, const TimeStep &ts) noexcept override;
    [[nodiscard]] bool executePostProcessingTasks(
        Context &ctx, const TimeStep &ts, const bool b) noexcept override;
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &ctx, const TimeStep &ts) noexcept override;
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    //! \brief destructor
    ~NonLinearModel() noexcept override;

   private:
    std::shared_ptr<NonLinearEvolutionProblem> problem;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARMODEL_HXX */
