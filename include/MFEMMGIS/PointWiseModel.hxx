/*!
 * \file   PointWiseModel.hxx
 * \brief  This file declares the `PointWiseModel` class
 * \author Thomas Helfer
 * \date   04/05/2026
 */

#ifndef LIB_MFEMMGIS_POINTWISEMODEL_HXX
#define LIB_MFEMMGIS_POINTWISEMODEL_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ModelBase.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  //! \brief a model integrating a MFront model at each integration point
  struct PointWiseModel : public ModelBase, protected Material {
    //! \return a description of the parameters of this model
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] qspace: partial quadrature space
     * \param[in] parameters: parameters
     */
    PointWiseModel(Context &ctx,
                   std::shared_ptr<const PartialQuadratureSpace> qspace,
                   const Parameters &parameters);
    //! \return the underlying material
    Material &getMaterial() noexcept;
    //! \return the underlying material
    const Material &getMaterial() const noexcept;
    //! \return the name of the MFront model
    [[nodiscard]] std::string getName() const noexcept override;
    /*!
     * \brief integrate the model at each integration point over the time step
     * \param[in, out] ctx: execution context
     * \param[in] ts: description of the time step
     * \return `ExitStatus::recoverableError` if the integration fails,
     * `ExitStatus::unreliableResults` if the model reports unreliable results,
     * `ExitStatus::success` otherwise, and an empty output
     */
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &ctx, const TimeStep &ts) noexcept override;
    /*!
     * \brief copy the state at the end of the time step on the state at the
     * beginning of the time step
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool update(Context &ctx) noexcept override;
    /*!
     * \brief copy the state at the beginning of the time step on the state at
     * the end of the time step
     * \param[in, out] ctx: execution context
     * \return true on success
     */
    [[nodiscard]] bool revert(Context &ctx) noexcept override;
    //! \brief destructor
    ~PointWiseModel() override;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POINTWISEMODEL_HXX */
