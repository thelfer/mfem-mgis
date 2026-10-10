/*!
 * \file   MFEMMGIS/PartialQuadratureFunctionsSet.hxx
 * \brief  This file declares the `PartialQuadratureFunctionsSet` class
 * \author Thomas Helfer
 * \date   02/06/2025
 */

#ifndef LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONSSET_HXX
#define LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONSSET_HXX

#ifndef MGIS_FUNCTION_SUPPORT
#error "PartialQuadratureFunctionsSet requires mgis/function"
#endif /* MGIS_FUNCTION_SUPPORT */

#include <vector>
#include <memory>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  /*!
   * \brief a structure grouping a set of partial quadrature functions.
   */
  struct MFEM_MGIS_EXPORT PartialQuadratureFunctionsSet
      : protected std::vector<std::shared_ptr<PartialQuadratureFunction>> {
    //! \brief prototype of a function able to update the functions of the set
    using UpdateFunction =
        std::function<bool(Context&, PartialQuadratureFunction&)>;
    //! \brief prototype of a function able to update the functions of the set
    using UpdateFunction2 = std::function<void(PartialQuadratureFunction&)>;
    /*!
     * \brief create the partial quadrature functions set using
     * the list of partial quadrature spaces and
     * the number of components
     *
     * \param[in] qspaces: partial quadrature spaces
     * \param[in] n: number of components
     */
    PartialQuadratureFunctionsSet(
        const std::vector<std::shared_ptr<const PartialQuadratureSpace>>&
            qspaces,
        const mfem_mgis::size_type n = 1);
    /*!
     * \brief create the partial quadrature functions set using
     * the given partial quadrature functions
     *
     * \param[in] functions: list of partial quadrature functions
     */
    PartialQuadratureFunctionsSet(
        const std::vector<std::shared_ptr<PartialQuadratureFunction>>&
            functions);
    //! \return the functions of the set
    std::vector<std::shared_ptr<const PartialQuadratureFunction>> getFunctions()
        const;
    //! \return the functions of the set
    const std::vector<std::shared_ptr<PartialQuadratureFunction>>&
    getFunctions();
    //! \return the list of locations
    std::vector<LocationIdentifier> getLocations() const noexcept;
    /*!
     * \brief return the partial quadrature function associated with the given
     * material identifier
     * \param[in, out] ctx: execution context
     * \param[in] l: location identifier
     * \return the partial quadrature function
     *
     * \note if no function associated with this identifier is found, a nullptr
     * is returned.
     */
    std::shared_ptr<PartialQuadratureFunction> get(
        Context& ctx, const LocationIdentifier l) noexcept;
    /*!
     * \brief return the partial quadrature function associated with the given
     * material identifier
     * \param[in, out] ctx: execution context
     * \param[in] l: location identifier
     * \return the partial quadrature function
     *
     * \note if no function associated with this identifier is found, a nullptr
     * is returned.
     */
    std::shared_ptr<const PartialQuadratureFunction> get(
        Context& ctx, const LocationIdentifier l) const noexcept;
    /*!
     * \brief update the set using an external function
     * \param[in, out] ctx: execution context
     * \param[in] f: function applied to each function of the set
     * \return true on success
     */
    bool update(Context& ctx, UpdateFunction& f);
    /*!
     * \brief update the set using an external function
     * \param[in] f: function applied to each function of the set
     */
    void update(UpdateFunction2& f);
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONSSET_HXX */
