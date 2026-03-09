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
    NonLinearModel(MeshDiscretization &, const Parameters &);
    NonLinearModel(std::shared_ptr<NonLinearEvolutionProblem>);
    //
    [[nodiscard]] NonLinearEvolutionProblem &getProblem() noexcept;
    [[nodiscard]] const NonLinearEvolutionProblem &getProblem() const noexcept;
    //
    [[nodiscard]] std::string getName() const noexcept override;
    [[nodiscard]] bool performInitializationTaksAtTheBeginningOfTheTimeStep(
        Context &, const TimeStep &) noexcept override;
    [[nodiscard]] bool executePostProcessingTasks(Context &,
                                                  const TimeStep &,
                                                  const bool) noexcept override;
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &, const TimeStep &) noexcept override;
    [[nodiscard]] bool update(Context &) noexcept override;
    [[nodiscard]] bool revert(Context &) noexcept override;
    //! \brief destructor
    ~NonLinearModel();

   private:
    std::shared_ptr<NonLinearEvolutionProblem> problem;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEARMODEL_HXX */
