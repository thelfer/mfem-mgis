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

  struct PointWiseModel : public ModelBase, protected Material {
    //! \return a description of the parameters of this model
    [[nodiscard]] static std::map<std::string, std::string>
    getParametersDescription() noexcept;
    /*!
     * \brief constructor
     * \param[in] qspace: partial quadrature space
     * \param[in] parameters: parameters
     */
    PointWiseModel(std::shared_ptr<const PartialQuadratureSpace>,
                   const Parameters &);
    //! \brief return the underlying material
    Material &getMaterial() noexcept;
    //! \brief return the underlying material
    const Material &getMaterial() const noexcept;
    //
    [[nodiscard]] std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
    computeNextState(Context &, const TimeStep &) noexcept override;
    [[nodiscard]] bool update(Context &) noexcept override;
    [[nodiscard]] bool revert(Context &) noexcept override;
    //! \brief destructor
    ~PointWiseModel() override;
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_POINTWISEMODEL_HXX */
