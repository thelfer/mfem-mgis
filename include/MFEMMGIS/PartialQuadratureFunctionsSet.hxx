/*!
 * \file   MFEMMGIS/PartialQuadratureFunctionsSet.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/06/2025
 */

#ifndef LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONSSET_HXX
#define LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONSSET_HXX

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
     * \param[in] nc: number of components
     */
    PartialQuadratureFunctionsSet(
        const std::vector<std::shared_ptr<const PartialQuadratureSpace>>&,
        const mfem_mgis::size_type = 1);
    /*!
     * \brief create the partial quadrature functions set using
     * the give partial quadrature functions
     *
     * \param[in] functions: list of quadrature
     * \param[in] nc: number of components
     */
    PartialQuadratureFunctionsSet(
        const std::vector<std::shared_ptr<PartialQuadratureFunction>>&);
    //! \brief return the functions of the set
    std::vector<std::shared_ptr<const PartialQuadratureFunction>> getFunctions()
        const;
    //! \brief return the functions of the set
    const std::vector<std::shared_ptr<PartialQuadratureFunction>>&
    getFunctions();
    //! \brief return the list of material identifiers
    std::vector<mfem_mgis::size_type> getMaterialIdentifiers() const;
    /*!
     * \brief return the partial quadrature function associated with the given
     * material identifier
     * \param[in] ctx: execution context
     * \param[in] m: material identifier
     *
     * \note if no function associated with this identifier is found, a nullptr
     * is returned.
     */
    std::shared_ptr<PartialQuadratureFunction> get(Context&,
                                                   const mfem_mgis::size_type);
    /*!
     * \brief return the partial quadrature function associated with the given
     * material identifier
     * \param[in] ctx: execution context
     * \param[in] m: material identifier
     *
     * \note if no function associated with this identifier is found, a nullptr
     * is returned.
     */
    std::shared_ptr<const PartialQuadratureFunction> get(
        Context&, const mfem_mgis::size_type) const;
    //! \brief update the set using an external function
    bool update(Context&, UpdateFunction&);
    //! \brief update the set using an external function
    void update(UpdateFunction2&);
  };

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARTIALQUADRATUREFUNCTIONSSET_HXX */
