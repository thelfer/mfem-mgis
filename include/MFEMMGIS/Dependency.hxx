/*!
 * \file   MFEMMGIS/Dependency.hxx
 * \brief  This file declares the DependencyBase and QPDependency classes
 * \author Thomas Helfer
 * \date   01/04/2026
 */

#ifndef LIB_MFEMMGIS_DEPENDENCY_HXX
#define LIB_MFEMMGIS_DEPENDENCY_HXX

#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct PartialQuadratureSpace;

  //! \brief a class describing a dependency at integration points
  struct MFEM_MGIS_EXPORT DependencyBase {
    //! \brief enumeration stating if the
    enum DependencyStatus {
      REQUIRED,
      OPTIONAL
    };
    //! \returns the name of the dependency
    [[nodiscard]] const std::string &getName() const noexcept;
    //! \return a string representation of the specifications
    [[nodiscard]] std::string getSpecificationsAsString() const noexcept;
    //! \return if the size have been specified
    [[nodiscard]] bool hasConcreteSpecifications() const noexcept;
    /*!
     * \brief check that this dependency is consistent with the given
     * specifications.
     *
     * \note if one of the dependency's specification is not
     * still optional, the specification is matched. If the specification has a
     * value, its value must match with the given value.
     *
     * \param[in] ctx: execution context
     * \param[in] nc: number of components
     */
    [[nodiscard]] bool matchesSpecifications(Context &ctx,
                                             const size_type nc) const noexcept;
    /*!
     * \brief check that this dependency is consistent with the given
     * specifications.
     *
     * \note this methods assumes that
     * `hasConcreteSpecifications` returns true
     *
     * \param[in] ctx: execution context
     * \param[in] nc: number of components
     */
    [[nodiscard]] bool checkSpecifications(Context &ctx,
                                           const size_type nc) const noexcept;
    /*!
     * \brief set the provider of the dependency
     * param[in] ctx: execution context
     * \param[in] p: provider
     */
    [[nodiscard]] bool setProvider(Context &ctx, const Provider &p) noexcept;
    //! \return if a this dependency has a provider
    [[nodiscard]] bool hasProvider() const noexcept;
    //! \return the provider of the dependency
    [[nodiscard]] OptionalReference<const Provider> getProvider(
        Context &ctx) const noexcept;
    //! \return if the dependency is required
    [[nodiscard]] bool isRequired() const noexcept;
    //! \return if the dependency is optional
    [[nodiscard]] bool isOptional() const noexcept;
    //! \return a description of the dependency for error reporting
    [[nodiscard]] virtual std::string getDescription() const noexcept = 0;

   protected:
    /*!
     * \brief constructor
     * \param[in] n: name
     * \param[in] s: status
     */
    DependencyBase(std::string_view n, const DependencyStatus s) noexcept;
    //! \brief copy constructor
    DependencyBase(const DependencyBase &) noexcept;
    //! \brief move constructor
    DependencyBase(DependencyBase &&) noexcept;
    //! \returns the number of components expected to be computed by the
    //! dependency
    [[nodiscard]] std::optional<size_type> getNumberOfComponents()
        const noexcept;
    /*!
     * \brief check that this dependency is consistent with other dependency and
     * update specifications (size) if the other dependency is more specific.
     *
     * \param[in] ctx: execution context
     * \param[in] d: dependency
     */
    [[nodiscard]] bool checkAndUpdate(Context &ctx,
                                      const DependencyBase &d) noexcept;
    /*!
     * \brief check that this dependency is consistent with the given
     * specifications and eventually update the specifications.
     *
     * \param[in] ctx: execution context
     * \param[in] n: name
     * \param[in] nc: number of components
     */
    [[nodiscard]] bool checkAndUpdate(Context &ctx,
                                      std::string_view n,
                                      const size_type nc) noexcept;
    /*!
     * \brief check that this dependency is consistent with the given
     * specifications and eventually update the specifications.
     *
     * \param[in] ctx: execution context
     * \param[in] n: number of components
     */
    [[nodiscard]] bool checkAndUpdate(Context &ctx,
                                      const size_type nc) noexcept;
    //! \brief destructor
    virtual ~DependencyBase() noexcept;
    //! \brief dependency status
    const DependencyStatus status;
    //! \brief name of the dependency
    const std::string name;
    //! \brief expected number of components
    std::optional<size_type> number_of_components;
    //! \brief provider
    const Provider *provider = nullptr;
  };  // end of DependencyBase

  //! \brief a class describing a dependency at quadrature points
  struct MFEM_MGIS_EXPORT QPDependency final : public DependencyBase {
    /*!
     * \brief report that the given provider requires the quadrature space to be
     * defined.
     *
     * \param[in] ctx: execution context
     * \param[in] d: dependency
     * \param[in] n: name of the provider
     */
    [[nodiscard]] static bool reportProviderRequiresQuadratureIdToBeDefined(
        Context &ctx, const QPDependency &d, const std::string &n) noexcept;
    /*!
     * \brief constructor
     * \param[in] s: partial quadrature space
     * \param[in] n: name of the dependency
     * \param[in] ds: status
     */
    QPDependency(
        std::shared_ptr<const PartialQuadratureSpace> s,
        std::string_view n,
        const DependencyStatus ds = DependencyStatus::REQUIRED) noexcept;
    /*!
     * \brief constructor
     * \param[in] m: material identifier
     * \param[in] n: name of the dependency
     * \param[in] s: status
     */
    QPDependency(
        const size_type m,
        std::string_view n,
        const DependencyStatus s = DependencyStatus::REQUIRED) noexcept;
    //! \brief copy constructor
    QPDependency(const QPDependency &) noexcept;
    //! \brief move constructor
    QPDependency(QPDependency &&) noexcept;
    //! \return the material identifier
    [[nodiscard]] size_type getMaterialIdentifier() const noexcept;
    //! \return the quadrature space
    [[nodiscard]] OptionalReference<const PartialQuadratureSpace>
    getPartialQuadratureSpace() const noexcept;
    //! \return the quadrature space
    [[nodiscard]] std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const noexcept;
    /*!
     * \brief set the quadrature id
     * \param[in] ctx: execution context
     * \param[in] qspace: partial quadrature space
     */
    [[nodiscard]] bool setPartialQuadratureSpace(
        Context &ctx, std::shared_ptr<const PartialQuadratureSpace> s) noexcept;
    /*!
     * \brief set this dependency as duplicate
     * \param[in] ctx: execution context
     */
    [[nodiscard]] bool setDuplicate(Context &ctx) noexcept;
    //! \return if this dependency is the duplicate of another one
    [[nodiscard]] bool isDuplicate() const noexcept;
    /*!
     * \brief check that this dependency is consistent with the given
     * specification (number of components) and update specifications if
     * required.
     *
     * \param[in] ctx: execution context
     * \param[in] nc: number of components
     */
    [[nodiscard]] bool checkAndUpdate(Context &ctx,
                                      const size_type nc) noexcept;
    /*!
     * \brief check that this dependency is consistent with other dependency
     * and update specifications (number of components) if the other
     * dependency is more specific.
     *
     * \param[in] ctx: execution context
     * \param[in] d: dependency
     */
    [[nodiscard]] bool checkAndUpdate(Context &ctx,
                                      const QPDependency &d) noexcept;
    //
    [[nodiscard]] std::string getDescription() const noexcept override;
    //! \brief destructor
    ~QPDependency() noexcept override;

   private:
    //! \brief material identifier
    const size_type material_identifier;
    //! \brief partial quadrature space
    std::shared_ptr<const PartialQuadratureSpace> qspace;
    /*!
     * \brief boolean stating if this dependency is the duplicate of another
     * dependency. This may happen if this dependency did not have
     * any initial quadrature id.
     */
    bool is_duplicate = false;
  };  // end of QPDependency

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_DEPENDENCY_HXX */