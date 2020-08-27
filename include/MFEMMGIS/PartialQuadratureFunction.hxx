/*!
 * \file   PartialQuadratureFunction.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX
#define LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX

#include <memory>
#include <vector>
#include "MGIS/Span.hxx"
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct PartialQuadratureSpace;

  /*!
   * \brief quadrature function defined on a partial quadrature space.
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureFunction {
    /*!
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] size: size of the data stored per integration points.
     */
    PartialQuadratureFunction(std::shared_ptr<const PartialQuadratureSpace>,
                              const size_type = 1);
    /*!
     * \brief return the data associated with integration point
     * \param[in] e: element number
     * \param[in] i: integration point number in the element
     */
    mgis::span<real> getIntegrationPointValues(const size_type, const size_type) ;
    /*!
     * \brief return the data associated with integration point
     * \param[in] e: element number
     * \param[in] i: integration point number in the element
     */
    mgis::span<const real> getIntegrationPointValues(const size_type, const size_type) const;

    //! \brief destructor
    ~PartialQuadratureFunction();

   private:
    //! underlying values
    std::vector<double> values;
    //! \brief underlying finite element space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
    //! \brief data size
    size_type data_size;
  };  // end of struct QuadratureFunction

}  // end of mfem_mgis

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX */
