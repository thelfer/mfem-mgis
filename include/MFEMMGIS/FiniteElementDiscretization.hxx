/*!
 * \file   include/MFEMMGIS/FiniteElementDiscretization.hxx
 * \brief
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#ifndef LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX
#define LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX

#include <map>
#include <string>
#include <vector>
#include <memory>
#include "MFEMMGIS/Info.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameter;
  struct Parameters;

  /*!
   * \brief a simple class used to:
   * - handle the life time of the mesh and the finite element collection.
   * - create and handle a finite element space
   */
  struct MFEM_MGIS_EXPORT FiniteElementDiscretization : MeshDiscretization {
    //! \brief string associated to the `FiniteElementFamily` parameter
    static const char* const FiniteElementFamily;
    //! \brief string associated to the `FiniteElementOrder` parameter
    static const char* const FiniteElementOrder;
    //! \brief string associated to the `UnknownsSize` parameter
    static const char* const UnknownsSize;
    //!
    [[noreturn]] static void reportInvalidParallelFiniteElementSpace();
    //!
    [[noreturn]] static void reportInvalidSequentialFiniteElementSpace();
    //!
    static std::vector<std::string> getParametersList();
    /*!
     * \brief constructor
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - parameters describing the mesh (see the `MeshDiscretization` class for
     *   details)
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. Supported families are:
     *   - `H1`
     *   The default value if `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `UnknownsSize` (int): number of components of the unknows
     */
    FiniteElementDiscretization(const Parameters&);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. Supported families are:
     *   - `H1`:
     *   The default value if `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `UnknownsSize` (int): number of components of the unknows
     */
    FiniteElementDiscretization(const MeshDiscretization&, const Parameters&);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. Supported families are:
     *   - `H1`:
     *   The default value if `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `UnknownsSize` (int): number of components of the unknows
     */
    FiniteElementDiscretization(std::shared_ptr<Mesh<true>>, const Parameters&);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. Supported families are:
     *   - `H1`:
     *   The default value if `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `UnknownsSize` (int): number of components of the unknows
     */
    FiniteElementDiscretization(std::shared_ptr<Mesh<false>>,
                                const Parameters&);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] c: finite element collection
     * \param[in] d: size of the unknowns
     *
     * \note this methods creates the finite element space.
     */
    FiniteElementDiscretization(std::shared_ptr<Mesh<true>>,
                                std::shared_ptr<const FiniteElementCollection>,
                                const size_type);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] c: collection
     * \param[in] s: finite element space
     */
    FiniteElementDiscretization(std::shared_ptr<Mesh<true>>,
                                std::shared_ptr<const FiniteElementCollection>,
                                std::unique_ptr<FiniteElementSpace<true>>);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] c: collection
     * \param[in] d: size of the unknowns
     *
     * \note this methods creates the finite element space.
     */
    FiniteElementDiscretization(std::shared_ptr<Mesh<false>>,
                                std::shared_ptr<const FiniteElementCollection>,
                                const size_type);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     * \param[in] c: collection
     * \param[in] s: finite element space
     */
    FiniteElementDiscretization(std::shared_ptr<Mesh<false>>,
                                std::shared_ptr<const FiniteElementCollection>,
                                std::unique_ptr<FiniteElementSpace<false>>);
    //! \return the finite element space
    template <bool parallel>
    [[nodiscard]] FiniteElementSpace<parallel>& getFiniteElementSpace();
    //! \return the finite element space
    template <bool parallel>
    [[nodiscard]] const FiniteElementSpace<parallel>& getFiniteElementSpace()
        const;
    //! \return the finite element collection
    [[nodiscard]] const FiniteElementCollection& getFiniteElementCollection()
        const noexcept;
    //! \return the finite element collection
    [[nodiscard]] std::shared_ptr<const FiniteElementCollection>
    getFiniteElementCollectionPointer() const noexcept;
    //
    // expose MeshDiscretization's methods, even deprecated ones for backward
    // compatibility
    //
    using MeshDiscretization::getBoundariesIdentifiers;
    using MeshDiscretization::getBoundariesNames;
    using MeshDiscretization::getBoundaryIdentifier;
    using MeshDiscretization::getBoundaryName;
    using MeshDiscretization::getMaterialIdentifier;
    using MeshDiscretization::getMaterialName;
    using MeshDiscretization::getMaterialsIdentifiers;
    using MeshDiscretization::getMaterialsNames;
    using MeshDiscretization::setBoundariesNames;
    using MeshDiscretization::setMaterialsNames;
    //! \brief destructor
    ~FiniteElementDiscretization();

   private:
    //! \brief finite element collection
    std::shared_ptr<const FiniteElementCollection> fec;
    //! \brief finite element space
#ifdef MFEM_USE_MPI
    std::unique_ptr<FiniteElementSpace<true>> parallel_fe_space;
#endif /* MFEM_USE_MPI */
    std::unique_ptr<FiniteElementSpace<false>> sequential_fe_space;
  };  // end of FiniteElementDiscretization

  /*!
   * \brief return the total number of unknowns
   * \param[in] fed: finite element discretization
   */
  MFEM_MGIS_EXPORT size_type getTrueVSize(const FiniteElementDiscretization&);

  /*!
   * \brief display information about a finite element discretization
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] fed: finite element discretization
   */
  template <>
  MFEM_MGIS_EXPORT bool getInformation<FiniteElementDiscretization>(
      Context&, std::ostream&, const FiniteElementDiscretization&) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/FiniteElementDiscretization.ixx"

#endif /* LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX */
