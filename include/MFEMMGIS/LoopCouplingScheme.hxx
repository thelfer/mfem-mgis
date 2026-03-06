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

  //! \brief the simpliest coupling scheme: all declared models are called once
  struct MFEM_MGIS_EXPORT LoopCouplingScheme : CouplingSchemeBase {
    //! \return a description of this scheme
    static std::string getDescription() noexcept;
    //! \return a description of the parameters of this scheme
    static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] m: mesh
     */
    LoopCouplingScheme(const MeshDiscretization &);
    /*!
     * \brief set the number of iterations
     * \param[in] ctx: execution context
     * \param[in] n: number of iterations
     */
    [[nodiscard]] bool setNumberOfIterations(Context &,
                                             const size_type) noexcept;
    //
    [[nodiscard]] std::string getName() const noexcept override;
    [[nodiscard]] std::optional<std::string> describe(
        Context &, const bool, const Parameters &) const noexcept override;
    //     [[nodiscard]] bool addConvergenceCriterion(
    //         Context &,
    //         std::string_view,
    //         const Parameters &) noexcept override final;
    [[nodiscard]] bool addConvergenceCriterion(
        Context &,
        std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion>) noexcept
        override final;
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &, const TimeStep &) noexcept override;
    //! \brief destructor
    ~LoopCouplingScheme() noexcept override;

   private:
    //! \brief number of iterations
    size_type number_of_iterations = 1;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LOOP_COUPLING_SCHEME_HXX */