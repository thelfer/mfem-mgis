/*!
 * \file   MFEMMGIS/AbstractPartialQuadratureFunctionEvaluator.hxx
 * \brief
 * \author Thomas HElfer
 * \date   06/03/2026
 */

#ifndef LIB_MFEMMGIS_ABSTRACTPARTIALQUADRATUREFUNCTIONEVALUATOR_HXX
#define LIB_MFEMMGIS_ABSTRACTPARTIALQUADRATUREFUNCTIONEVALUATOR_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct PartialQuadratureSpace;
  struct PartialQuadratureFunction;

  struct MFEM_MGIS_EXPORT AbstractPartialQuadratureFunctionEvaluator {
    [[nodiscard]] virtual const PartialQuadratureSpace& getQuadratureSpace()
        const noexcept = 0;
    [[nodiscard]] virtual std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const noexcept = 0;
    //! \return the number of components
    [[nodiscard]] virtual size_type getNumberOfComponents() const noexcept;
    //!
    [[nodiscard]] virtual bool evaluate(
        Context&, PartialQuadratureFunction&) const noexcept = 0;
    //! \brief destructor
    virtual ~AbstractPartialQuadratureFunctionEvaluator() noexcept;
  };  // end of AbstractPartialQuadratureFunctionEvaluator

  MFEM_MGIS_EXPORT [[nodiscard]] std::shared_ptr<PartialQuadratureFunction>
  evaluate(Context&,
           const AbstractPartialQuadratureFunctionEvaluator&) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTPARTIALQUADRATUREFUNCTIONEVALUATOR_HXX */
