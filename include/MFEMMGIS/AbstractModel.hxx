/*!
 * \file   MFEMMGIS/AbstractModel.hxx
 * \brief  This file declares the `AbstractModel` class
 * \date   15/11/2022
 */

#ifndef LIB_MFEM_MGIS_ABSTRACT_MODEL_HXX
#define LIB_MFEM_MGIS_ABSTRACT_MODEL_HXX

#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Provider.hxx"
#include "MFEMMGIS/AbstractCouplingItem.hxx"

namespace mfem_mgis {

  /*!
   * \brief an abstract class for models
   *
   * Most models are built on balance equations
   * (heat transfer, mechanics, diffusion, etc..) which require to solve
   * partial differential equations (PDEs).
   *
   * Models may also define internal state variables which are generally
   * defined at integration points.
   *
   * Unknown fields and internal state variables are generally sufficient to
   * define the state of the model.
   */
  struct MFEM_MGIS_EXPORT AbstractModel : AbstractCouplingItem, Provider {
    /*!
     * \brief add a new post-processing
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the post-processing
     * \param[in] params: parameters defining the post-processing
     * \return true on success
     */
    [[nodiscard]] virtual bool addPostProcessing(
        Context& ctx,
        std::string_view n,
        const Parameters& params) noexcept = 0;
    //! \return the list of available post-processings
    [[nodiscard]] virtual std::vector<std::string> getAvailablePostProcessings()
        const noexcept = 0;
    //! \brief destructor
    ~AbstractModel() noexcept override;
  };  // end of  class AbstractModel

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACT_MODEL_HXX */
