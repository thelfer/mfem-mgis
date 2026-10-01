/*!
 * \file   MFEMMGIS/AbstractCouplingScheme.hxx
 * \brief  This file declares the `AbstractCouplingScheme` class
 * \date   05/12/2022
 */

#ifndef LIB_MFEM_MGIS_ABSTRACT_COUPLING_SCHEME_HXX
#define LIB_MFEM_MGIS_ABSTRACT_COUPLING_SCHEME_HXX

#include <vector>
#include <memory>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractCouplingItem.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  struct Provider;
  struct AbstractModel;
  struct NonLinearEvolutionProblem;
  struct AbstractCouplingSchemeConvergenceCriterion;

  //! \brief interface of all coupling schemes
  struct MFEM_MGIS_EXPORT AbstractCouplingScheme : AbstractCouplingItem {
    //! \return the list of providers handled by the coupling scheme
    [[nodiscard]] virtual std::vector<const Provider*>
    getProviders() noexcept = 0;
    //     /*
    //      * \brief add a new coupling item
    //      * \param[in, out] ctx: execution context
    //      * \param[in] parameters: parameters describing the item
    //      * \return true on success
    //      *
    //      * \note the parameters must contain a dictionary with one element.
    //      * The key gives the type of the item inserted, for example `Model`
    //      * or `model`. The value must be a dictionary describing the item.
    //      * This description must be a dictionary with one element. Its key
    //      * is the name of the item and its value is a dictionary describing
    //      * the parameters of the item to be created.
    //      *
    //      * A typical example, in JSON notation, is:
    //      * {"Model": {"MFrontPointWiseModel": {
    //      *     "Material": "mesh", "Model": "ChemicalReaction5"}}}
    //      */
    //     [[nodiscard]] virtual bool add(Context &, const Parameters &)
    //     noexcept = 0;
    //     /*
    //      * \brief add a new coupling item
    //      *
    //      * \param[in, out] ctx: execution context
    //      * \param[in] t: type of the coupling item, one of `model`, `Model`,
    //      * `post-processing` and `PostProcessing`
    //      * \param[in] n: name of the coupling item
    //      * \param[in] parameters: parameters passed to the constructor of the
    //      * coupling item
    //      * \return true on success
    //      */
    //     [[nodiscard]] virtual bool addCouplingItem(Context &,
    //                                                std::string_view,
    //                                                std::string_view,
    //                                                const Parameters &)
    //                                                noexcept = 0;
    /*!
     * \brief add a new coupling item
     * \param[in, out] ctx: execution context
     * \param[in] i: coupling item
     * \return true on success
     */
    [[nodiscard]] virtual bool addCouplingItem(
        Context& ctx, std::shared_ptr<AbstractCouplingItem> i) noexcept = 0;
    //     /*
    //      * \brief add a new model
    //      * \param[in, out] ctx: execution context
    //      * \param[in] n: name of the model
    //      * \param[in] parameters: parameters passed to the model
    //      * \return true on success
    //      */
    //     [[nodiscard]] virtual bool addModel(Context &,
    //                                         std::string_view,
    //                                         const Parameters &) noexcept = 0;
    /*!
     * \brief add a new model
     * \param[in, out] ctx: execution context
     * \param[in] m: model
     * \return true on success
     */
    [[nodiscard]] virtual bool addModel(
        Context& ctx, std::shared_ptr<AbstractModel> m) noexcept = 0;
    /*!
     * \brief add a new model
     * \param[in, out] ctx: execution context
     * \param[in] m: non linear evolution problem
     * \return true on success
     */
    [[nodiscard]] virtual bool addModel(
        Context& ctx,
        std::shared_ptr<NonLinearEvolutionProblem> m) noexcept = 0;
    /*!
     * \brief add a new convergence criterion
     * \param[in, out] ctx: execution context
     * \param[in] c: convergence criterion
     * \return true on success
     */
    [[nodiscard]] virtual bool addConvergenceCriterion(
        Context& ctx,
        std::shared_ptr<AbstractCouplingSchemeConvergenceCriterion>
            c) noexcept = 0;
    //     /*
    //      * \brief add a new convergence criterion
    //      * \param[in, out] ctx: execution context
    //      * \param[in] n: name of the convergence criterion
    //      * \param[in] parameters: parameters passed to the convergence
    //      * criterion
    //      * \return true on success
    //      */
    //     [[nodiscard]] virtual bool addConvergenceCriterion(
    //         Context &, std::string_view, const Parameters &) noexcept = 0;
    //! \brief destructor
    ~AbstractCouplingScheme() noexcept override;
  };  // end of AbstractCouplingScheme

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_ABSTRACT_COUPLING_SCHEME_HXX */
