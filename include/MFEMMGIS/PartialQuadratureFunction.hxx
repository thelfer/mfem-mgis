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
#include <limits>
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
     * \brief constructor
     * \param[in] s: quadrature space.
     * \param[in] v: values
     * \param[in] db: start of the view inside the given data
     * \param[in] ds: size of the view
     */
    PartialQuadratureFunction(
        std::shared_ptr<const PartialQuadratureSpace>,
        mgis::span<real>,
        const size_type = 0,
        const size_type = std::numeric_limits<size_type>::max());
    /*!
     * \brief return the value associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    const real& getIntegrationPointValue(const size_type) const;
    /*!
     * \brief return the value associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    real& getIntegrationPointValue(const size_type, const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     * \note this method is only meaningful when the quadrature function is
     * scalar
     */
    const real& getIntegrationPointValue(const size_type,
                                         const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type N>
    mgis::span<real, N> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    template <size_type N>
    mgis::span<const real, N> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    mgis::span<real> getIntegrationPointValues(const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] o: offset associated with the integration point
     */
    mgis::span<const real> getIntegrationPointValues(const size_type) const;
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    mgis::span<real> getIntegrationPointValues(const size_type,
                                               const size_type);
    /*!
     * \brief return the data associated with an integration point
     * \param[in] e: global element number
     * \param[in] i: integration point number in the element
     */
    mgis::span<const real> getIntegrationPointValues(const size_type,
                                                     const size_type) const;

    //! \brief destructor
    ~PartialQuadratureFunction();

   protected:
    /*!
     * \return the data offset associated with the given integration point.
     * \param[in] o: offset associated with the integration point
     */
    size_type getDataOffset(const size_type) const;
    //! \brief underlying finite element space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
    //! \brief underlying values
    mgis::span<real> values;
    //! \brief storage for the values when the partial function holds the
    //! values
    std::vector<real> values_storage;
    //! \brief data stride
    size_type data_stride;
    //! \brief begin of the data
    size_type data_begin;
    //! \brief data size
    size_type data_size;
    };  // end of struct PartialQuadratureFunction

}  // namespace mfem_mgis

#include "MFEMMGIS/PartialQuadratureFunction.ixx"

#endif /* LIB_MFEM_MGIS_PARTIALQUADRATUREFUNCTION_HXX */
