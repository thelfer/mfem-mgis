/*!
 * \file   include/MFEMMGIS/MeshDiscretization.hxx
 * \brief
 * \author Thomas Helfer
 * \date   06/03/2026
 */

#ifndef LIB_MFEM_MGIS_MESHDISCRETIZATION_HXX
#define LIB_MFEM_MGIS_MESHDISCRETIZATION_HXX

#include <map>
#include <string>
#include <vector>
#include <memory>
#include "MFEMMGIS/Info.hxx"
#include "MFEMMGIS/Config.hxx"
#ifdef MGIS_HAVE_TFEL
#include "MFEMMGIS/Geometry.hxx"
#endif /* MGIS_HAVE_TFEL */

namespace mfem_mgis {

  // forward declarations
  struct Parameter;
  struct Parameters;

  //! \brief a simple class used to handle the life time of the mesh
  struct MFEM_MGIS_EXPORT [[nodiscard]] MeshDiscretization {
    //! \brief structure holding a list of attributes to be used as a key.
    struct AttributesList {
      /*!
       * \brief constructor
       * \param[ids] ids: list of attributes
       */
      AttributesList(const std::vector<size_type>& ids) : attributes(ids) {
        std::sort(this->attributes.begin(), this->attributes.end());
      }  // end of AttributesList
      //! \brief comparison operator
      [[nodiscard]] bool operator<(const AttributesList& rhs) const noexcept {
        return this->attributes < rhs.attributes;
      }  // end of operator<
      [[nodiscard]] const std::vector<size_type>& getAttributes()
          const noexcept {
        return this->attributes;
      }

     private:
      //! \brief list of attributes
      std::vector<size_type> attributes;
    };
    //! \brief location on which submeshes can be defined
    enum struct Location { ON_MATERIALS, ON_BOUNDARIES };
    /*!
     * \brief an helper structure used to distinguish the identifiers of a
     * single material (associated with a unique attribute) and indentifiers of
     * a single boundary (associated with a unique boundary attribute)
     */
    template <MeshDiscretization::Location>
    struct RawLocationIdentifier {
      //! \brief attribute associated with the material or the boundary
      size_type id;
      //! \brief comparisons operator
      constexpr auto operator<=>(const RawLocationIdentifier&) const noexcept =
          default;
    };
    /*!
     * \brief a simple alias for the identifier of a single
     * material (associated with a unique attribute)
     */
    using MaterialIdentifier =
        RawLocationIdentifier<MeshDiscretization::Location::ON_MATERIALS>;
    /*!
     * \brief a simple alias for the identifier of a single
     * boundary (associated with a unique attribute)
     */
    using BoundaryIdentifier =
        RawLocationIdentifier<MeshDiscretization::Location::ON_BOUNDARIES>;
    /*!
     * \brief a simple structure to store either a material identifier or a
     * boundary identifier.
     *
     * \note This identifier must refer to the main mesh. See
     * the getLocationIdentifier in `MeshDiscretization` for details.
     *
     * \note This structure is sortable, and can be used as key in standard
     * associative containers.
     *
     * \note this structure may be invalid and works with MGIS's error handling
     * scheme.
     */
    struct LocationIdentifier {
      //! \brief identifier associated with a material
      std::optional<MaterialIdentifier> material_identifier;
      //! \brief identifier associated with a boundary
      std::optional<BoundaryIdentifier> boundary_identifier;
      //! \brief comparisons operator
      constexpr auto operator<=>(const LocationIdentifier&) const noexcept =
          default;
    };  // end of LocationIdentifier
    //! \brief string associated to the `Parallel` parameter
    static const char* const Parallel;
    //! \brief string associated to the `MeshFileName` parameter
    static const char* const MeshFileName;
    //! \brief string associated to the `MeshReadMode` parameter
    static const char* const MeshReadMode;
    //! \brief string associated to the `Materials` parameter
    static const char* const Materials;
    //! \brief string associated to the `Boundaries` parameter
    static const char* const Boundaries;
    //! \brief string associated to the `Points` parameter
    static const char* const Points;
    //! \brief string associated to the `PointsSets` parameter
    static const char* const PointsSets;
    //! \brief string associated to the `NumberOfUniformRefinements` parameter
    static const char* const NumberOfUniformRefinements;
    //! \brief string associated to the `VerbosityLevel` parameter
    static const char* const GeneralVerbosityLevel;
    //!
    [[nodiscard]] static std::vector<std::string> getParametersList() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context used for profiling
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - `Parallel` (boolean): if true, a parallel computation is to be be
     *    performed. This value if assumed to be false by default.
     * - `MeshFileName` (string): mesh file.
     * - `NumberOfUniformRefinements` (int): number of uniform refinements
     *   applied to the mesh
     * - `GeneralVerbosityLevel` (int): with large positive numbers, expect more
     * verbosity
     */
    MeshDiscretization(mgis::Context&, const Parameters&);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     */
    MeshDiscretization(std::shared_ptr<Mesh<true>>);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     */
    MeshDiscretization(std::shared_ptr<Mesh<false>>);
    //! \brief move constructor
    MeshDiscretization(MeshDiscretization&&) noexcept;
    //! \brief copy constructor
    MeshDiscretization(const MeshDiscretization&) noexcept;
    /*!
     * \return if the given mesh is managed by this mesh discretization
     * \param[in] m: parallel mesh
     */
    bool manages(const Mesh<true>&) const noexcept;
    /*!
     * \return if the given mesh is managed by this mesh discretization
     * \param[in] m: sequential mesh
     */
    bool manages(const Mesh<false>&) const noexcept;
    /*!
     * \return if the given mesh is defined on (a subset of) the
     * materials of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this methods fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnMaterials(
        Context& ctx, const Mesh<true>&) const noexcept;
    /*!
     * \return if the given mesh is defined on (a subset of) the
     * materials of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this methods fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnMaterials(
        Context& ctx, const Mesh<false>&) const noexcept;
    /*!
     * \return if the given mesh is defined on (a subset of) the
     * boundaries of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this methods fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnBoundaries(
        Context& ctx, const Mesh<true>&) const noexcept;
    /*!
     * \return if the given mesh is defined on (a subset of) the
     * boundaries of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this methods fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnBoundaries(
        Context& ctx, const Mesh<false>&) const noexcept;
    /*!
     * \return a pointer to the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    template <bool parallel>
    std::shared_ptr<SubMesh<parallel>> getMutableSubMeshPointer(
        Context&, const Parameter&, const Location) const noexcept;
    /*!
     * \return a pointer to the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    template <bool parallel>
    std::shared_ptr<SubMesh<parallel>> getSubMeshPointer(
        Context&, const Parameter&, const Location) noexcept;
    /*!
     * \return a pointer to the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    template <bool parallel>
    std::shared_ptr<const SubMesh<parallel>> getSubMeshPointer(
        Context&, const Parameter&, const Location) const noexcept;
    /*!
     * \return the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    template <bool parallel>
    OptionalReference<SubMesh<parallel>> getMutableSubMeshReference(
        Context&, const Parameter&, const Location) const noexcept;
    /*!
     * \return the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    template <bool parallel>
    OptionalReference<SubMesh<parallel>> getSubMesh(Context&,
                                                    const Parameter&,
                                                    const Location) noexcept;
    /*!
     * \return the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    template <bool parallel>
    OptionalReference<const SubMesh<parallel>> getSubMesh(
        Context&, const Parameter&, const Location) const noexcept;
    /*!
     * \return the shared pointer associated with the given mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    template <bool parallel>
    std::shared_ptr<Mesh<parallel>> getMutableMeshPointer(
        Context&, const Mesh<parallel>&) const noexcept;
    /*!
     * \return the shared pointer associated with the given mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    template <bool parallel>
    std::shared_ptr<Mesh<parallel>> getMeshPointer(
        Context&, const Mesh<parallel>&) noexcept;
    /*!
     * \return the shared pointer associated with the given mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    template <bool parallel>
    std::shared_ptr<const Mesh<parallel>> getMeshPointer(
        Context&, const Mesh<parallel>&) const noexcept;
#ifdef MFEM_USE_MPI
    /*!
     * \return the location identifier in the main mesh
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] id: attribute in the mesh
     *
     * \note This rationale behind this method is that submesh may be created on
     * boundaries. In this case, the boundary attributes used to create
     * the boundaries becomes standard attributes of the submesh. This may lead
     * to ambiguity when its comes to determine where a partial quadrature space
     * is defined for instance. The returned location identifier does not have
     * such ambiguity.
     */
    [[nodiscard]] std::optional<LocationIdentifier> getLocationIdentifier(
        Context&, const Mesh<true>&, const size_type) const noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \return the location identifier in the main mesh
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] id: attribute in the mesh
     *
     * \note This rationale behind this method is that submesh may be created on
     * boundaries. In this case, the boundary attributes used to create
     * the boundaries becomes standard attributes of the submesh. This may lead
     * to ambiguity when its comes to determine where a partial quadrature space
     * is defined for instance. The returned location identifier does not have
     * such ambiguity.
     */
    [[nodiscard]] std::optional<LocationIdentifier> getLocationIdentifier(
        Context&, const Mesh<false>&, const size_type) const noexcept;
    /*!
     * \brief set material names
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     */
    [[nodiscard]] bool setMaterialsNames(
        Context&, const std::map<size_type, std::string>&) noexcept;
    /*!
     * \brief set material names
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     */
    [[nodiscard]] bool setBoundariesNames(
        Context&, const std::map<size_type, std::string>&) noexcept;
    /*!
     * \return the name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \param[in, out] ctx: execution context
     * \param[in] id: location identifier
     * \note the method only fails if the identifier is invalid or not defined
     * in the mesh
     */
    [[nodiscard]] std::optional<std::string> getLocationName(
        Context&, const LocationIdentifier&) const noexcept;
    /*!
     * \return the material name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \param[in, out] ctx: execution context
     * \param[in] id: material identifier
     * \note the method only fails if the material identifier is not defined in
     * the mesh
     */
    [[nodiscard]] std::optional<std::string> getMaterialName(
        Context&, const size_type) const noexcept;
    /*!
     * \return the boundary name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \param[in, out] ctx: execution context
     * \param[in] id: boundary identifier
     * \note the method only fails if the boundary identifier is not defined in
     * the mesh
     */
    [[nodiscard]] std::optional<std::string> getBoundaryName(
        Context&, const size_type) const noexcept;
    /*!
     * \return the material identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] std::optional<size_type> getMaterialIdentifier(
        Context&, const Parameter&) const noexcept;
    /*!
     * \return the boundary identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] std::optional<size_type> getBoundaryIdentifier(
        Context&, const Parameter&) const noexcept;
    /*!
     * \return the list of materials identifiers described by the given
     * parameter.
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings and integers.
     *
     * Integers are directly intepreted as materials identifiers.
     *
     * Strings are intepreted as regular expressions which allows the selection
     * of materials by names.
     */
    [[nodiscard]] std::optional<std::vector<size_type>> getMaterialsIdentifiers(
        Context&, const Parameter&) const noexcept;
    /*!
     * \return the list of boundaries identifiers described by the given
     * parameter.
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings and integers.
     *
     * Integers are directly intepreted as boundaries identifiers.
     *
     * Strings are intepreted as regular expressions which allows the selection
     * of boundaries by names.
     */
    [[nodiscard]] std::optional<std::vector<size_type>>
    getBoundariesIdentifiers(Context&, const Parameter&) const noexcept;
    /*!
     * \return the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     */
    template <bool parallel>
    [[nodiscard]] Mesh<parallel>& getMesh() noexcept;
    /*!
     * \return the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     */
    template <bool parallel>
    [[nodiscard]] const Mesh<parallel>& getMesh() const noexcept;
    /*!
     * \return a pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<Mesh<parallel>> getMeshPointer() noexcept;
    /*!
     * \return a pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<const Mesh<parallel>> getMeshPointer()
        const noexcept;
    /*!
     * \return a mutable pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<Mesh<parallel>> getMutableMeshPointer()
        const noexcept;
    //! \return if this object is built to run parallel computations
    [[nodiscard]] bool describesAParallelComputation() const noexcept;
    /*!
     * \return the names of the materials (and their mapping with their
     * identifiers
     */
    [[nodiscard]] std::map<size_type, std::string> getMaterialsNames()
        const noexcept;
    /*!
     * \return the names of the boundaries (and their mapping with their
     * identifiers
     */
    [[nodiscard]] std::map<size_type, std::string> getBoundariesNames()
        const noexcept;
#ifdef MGIS_HAVE_TFEL
    /*!
     * \brief add point
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \param[in] pt: coordinates of the point
     */
    [[nodiscard]] bool addPoint(Context&,
                                std::string_view,
                                const Point<2>&) noexcept;
    /*!
     * \brief add point
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \param[in] pt: coordinates of the point
     */
    [[nodiscard]] bool addPoint(Context&,
                                std::string_view,
                                const Point<3>&) noexcept;
    /*!
     * \brief add point set
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point set
     * \param[in] pts: list of points
     */
    [[nodiscard]] bool addPointsSet(Context&,
                                    std::string_view,
                                    const std::vector<Point<2>>&) noexcept;
    /*!
     * \brief add point set
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point set
     * \param[in] pts: list of points
     */
    [[nodiscard]] bool addPointsSet(Context&,
                                    std::string_view,
                                    const std::vector<Point<3>>&) noexcept;
    /*!
     * \return the point with the given name
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] std::optional<Point<N>> getPoint(
            Context&, std::string_view) const noexcept;
    /*!
     * \return the registered points
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] OptionalReference<
            const std::map<std::string, Point<N>, std::less<>>>  //
        getPoints(Context&) const noexcept;
    /*!
     * \return the registered points sets
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] OptionalReference<
            const std::map<std::string, std::vector<Point<N>>, std::less<>>>  //
        getPointsSets(Context&) const noexcept;
    /*!
     * \return the points set with the given name
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     */
    template <size_type N>
    requires((N == 2) || (N == 3))                                    //
        [[nodiscard]] OptionalReference<const std::vector<Point<N>>>  //
        getPointsSet(Context&, std::string_view) const noexcept;
#endif /* MGIS_HAVE_TFEL */

    //! \brief destructor
    ~MeshDiscretization();

   protected:
    // friend functions and operators
    friend bool getInformation<MeshDiscretization>(
        Context&, std::ostream&, const MeshDiscretization&) noexcept;
    friend bool operator==(const MeshDiscretization&,
                           const MeshDiscretization&) noexcept;
    //! \return a mutable pointer to the underlying parallel mesh
    [[nodiscard]] std::shared_ptr<Mesh<true>> getMutableParallelMeshPointer()
        const noexcept;
    //! \return a mutable pointer to the underlying sequential mesh
    [[nodiscard]] std::shared_ptr<Mesh<false>> getMutableSequentialMeshPointer()
        const noexcept;
    //! \return a pointer to the underlying parallel mesh
    [[nodiscard]] std::shared_ptr<const Mesh<true>> getParallelMeshPointer()
        const noexcept;
    //! \return a pointer to the underlying sequential mesh
    [[nodiscard]] std::shared_ptr<const Mesh<false>> getSequentialMeshPointer()
        const noexcept;
    /*!
     * \return the location identifier in the main mesh
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] id: attribute in the mesh
     *
     * \see `MeshDiscretization::getLocationIdentifier` for details
     */
    [[nodiscard]] std::optional<LocationIdentifier>
    getParallelLocationIdentifier(Context&,
                                  const Mesh<true>&,
                                  const size_type) const noexcept;
    /*!
     * \return the location identifier in the main mesh
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] id: attribute in the mesh
     *
     * \see `MeshDiscretization::getLocationIdentifier` for details
     */
    [[nodiscard]] std::optional<LocationIdentifier>
    getSequentialLocationIdentifier(Context&,
                                    const Mesh<false>&,
                                    const size_type) const noexcept;
    /*!
     * \return a mutable pointer to the underlying parallel mesh
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] std::shared_ptr<Mesh<true>> getMutableParallelMeshPointer(
        Context&, const Mesh<true>&) const noexcept;
    /*!
     * \return a mutable pointer to the underlying sequential mesh
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    [[nodiscard]] std::shared_ptr<Mesh<false>> getMutableSequentialMeshPointer(
        Context&, const Mesh<false>&) const noexcept;
    /*!
     * \return a pointer to the underlying parallel mesh
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    [[nodiscard]] std::shared_ptr<const Mesh<true>> getParallelMeshPointer(
        Context&, const Mesh<true>&) const noexcept;
    /*!
     * \return a pointer to the underlying sequential mesh
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    [[nodiscard]] std::shared_ptr<const Mesh<false>> getSequentialMeshPointer(
        Context&, const Mesh<false>&) const noexcept;
    /*!
     * \return the parallel sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    std::shared_ptr<const SubMesh<true>> getParallelSubMeshPointer(
        Context&, const Parameter&, const Location) const noexcept;
    /*!
     * \return the parallel sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    std::shared_ptr<SubMesh<true>> getParallelMutableSubMeshPointer(
        Context&, const Parameter&, const Location) const noexcept;
    /*!
     * \return the sequential sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    std::shared_ptr<const SubMesh<false>> getSequentialSubMeshPointer(
        Context&, const Parameter&, const Location) const noexcept;
    /*!
     * \return the sequential sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    std::shared_ptr<SubMesh<false>> getSequentialMutableSubMeshPointer(
        Context&, const Parameter&, const Location) const noexcept;

#ifdef MGIS_HAVE_TFEL
    /*!
     * \return the registred points in 2D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, Point<2>, std::less<>>>
    getPoints2D(Context&) const noexcept;
    /*!
     * \return the registred points in 3D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, Point<3>, std::less<>>>
    getPoints3D(Context&) const noexcept;
    /*!
     * \return the registred points in 2D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, std::vector<Point<2>>, std::less<>>>
    getPointsSets2D(Context&) const noexcept;
    /*!
     * \return the registred points in 3D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, std::vector<Point<3>>, std::less<>>>
    getPointsSets3D(Context&) const noexcept;
    /*!
     * \return the point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     */
    [[nodiscard]] std::optional<Point<2>> getPoint2D(
        Context&, std::string_view) const noexcept;
    /*!
     * \return the point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     */
    [[nodiscard]] std::optional<Point<3>> getPoint3D(
        Context&, std::string_view) const noexcept;
    /*!
     * \return the points set with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     */
    [[nodiscard]] OptionalReference<const std::vector<Point<2>>> getPointsSet2D(
        Context&, std::string_view) const noexcept;
    /*!
     * \return the point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     */
    [[nodiscard]] OptionalReference<const std::vector<Point<3>>> getPointsSet3D(
        Context&, std::string_view) const noexcept;
#endif /* MGIS_HAVE_TFEL */
    //! \brief internal structure to implement the pimpl idiom
    struct Implementation;
    //! \brief pointer to the internal implementation
    std::shared_ptr<Implementation> pimpl;
  };  // end of MeshDiscretization

  using MaterialIdentifier = MeshDiscretization::MaterialIdentifier;
  /*!
   * \brief a simple alias for the identifier of a single
   * boundary (associated with a unique attribute)
   */
  using BoundaryIdentifier = MeshDiscretization::BoundaryIdentifier;
  /*!
   * \brief a simple alias for the identifier of a single material or a
   * single boundary (associated with a unique attribute)
   */
  using LocationIdentifier = MeshDiscretization::LocationIdentifier;
  /*!
   * \return if the given location identifier is invalid
   * \param[in] l: location identifier
   */
  [[nodiscard]] constexpr bool isInvalid(const LocationIdentifier& l) noexcept {
    const auto mok = isValid(l.material_identifier);
    const auto bok = isValid(l.boundary_identifier);
    const auto b1 = (!mok) && (!bok);  // none is valid
    const auto b2 = mok && bok;        // both are valid
    return b1 || b2;
  }  // end of is Invalid
  /*!
   * \brief compare two mesh discretisations to see if they points to the same
   * underlying implementation
   *
   * \param[in] lhs: left hand side
   * \param[in] rhs: right hand side
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool operator==(
      const MeshDiscretization&, const MeshDiscretization&) noexcept;
  /*!
   * \brief compare two mesh discretisations
   *
   * \param[in] lhs: left hand side
   * \param[in] rhs: right hand side
   *
   * \note material names and boundary names may different in both
   * discretisations.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool operator!=(
      const MeshDiscretization&, const MeshDiscretization&) noexcept;
  /*!
   * \return the space dimension
   * \param[in] m: mesh discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getSpaceDimension(
      const MeshDiscretization&) noexcept;
  /*!
   * \return the list of materials attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getMaterialsAttributes(const MeshDiscretization&) noexcept;
  /*!
   * \return the list of boundaries attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getBoundariesAttributes(const MeshDiscretization&) noexcept;

  /*!
   * \return the list of materials identifiers described by the given
   * parameter.
   * \param[in] throwing: throwing attributes
   * \param[in] m: mesh discretization
   * \param[in] p: parameter
   *
   * \note The parameter may hold:
   *
   * - an integer
   * - a string
   * - a vector of parameters which must be either strings and integers.
   *
   * Integers are directly intepreted as materials identifiers.
   *
   * Strings are intepreted as regular expressions which allows the selection
   * of materials by names.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<size_type> getMaterialsIdentifiers(
      attributes::Throwing, const MeshDiscretization&, const Parameter&);
  /*!
   * \return the list of boundaries identifiers described by the given
   * parameter.
   * \param[in] throwing: throwing attributes
   * \param[in] m: mesh discretization
   * \param[in] p: parameter
   *
   * \note The parameter may hold:
   *
   * - an integer
   * - a string
   * - a vector of parameters which must be either strings and integers.
   *
   * Integers are directly intepreted as boundaries identifiers.
   *
   * Strings are intepreted as regular expressions which allows the selection
   * of boundaries by names.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<size_type>
  getBoundariesIdentifiers(attributes::Throwing,
                           const MeshDiscretization&,
                           const Parameter&);

  /*!
   * \brief display information about a mesh discretization
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] m: mesh discretization
   */
  template <>
  MFEM_MGIS_EXPORT bool getInformation<MeshDiscretization>(
      Context&, std::ostream&, const MeshDiscretization&) noexcept;

#ifdef MFEM_USE_MPI

  /*!
   * \return the MPI communicator associated with the mesh discretization
   * \param[in] m: mesh discretization
   *
   * \note If a sequential computation is described, `MPI_COMM_WORLD` is
   * returned.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] MPI_Comm getMPICommunicator(
      const MeshDiscretization&) noexcept;

#endif /* MFEM_USE_MPI */

  /*!
   * \return if the current process is the main one (the process of rank 0)
   * \param[in] m: mesh discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isMainProcess(
      const MeshDiscretization&) noexcept;

  /*!
   * \return if the given location identifier is consistent with the mesh
   * discretization
   *
   * \param[in,out] ctx: execution context
   * \param[in] m: mesh discretization
   * \param[in] l: location identifier
   *
   * This check fails if:
   *
   * - the identifier is invalid
   * - if the material identifier (if valid) is not mesh attribute
   * - if the boundary identifier (if valid) is not boundary mesh attribute
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool check(Context&,
                                            const MeshDiscretization&,
                                            const LocationIdentifier&) noexcept;

#ifdef MGIS_HAVE_TFEL

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<Point<N>> makePoint(
          Context&, const MeshDiscretization&, const Parameter&) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsSet(
          Context&, const MeshDiscretization&, const Parameter&) noexcept;

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsOnCurve(
          Context&, const MeshDiscretization&, const Parameters&) noexcept;

  // partial specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<2>> makePoint<2>(
      Context&, const MeshDiscretization&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<3>> makePoint<3>(
      Context&, const MeshDiscretization&, const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsSet<2>(Context&,
                   const MeshDiscretization&,
                   const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsSet<3>(Context&,
                   const MeshDiscretization&,
                   const Parameter&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsOnCurve<2>(Context&,
                       const MeshDiscretization&,
                       const Parameters&) noexcept;
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsOnCurve<3>(Context&,
                       const MeshDiscretization&,
                       const Parameters&) noexcept;

#endif /* MGIS_HAVE_TFEL */

}  // end of namespace mfem_mgis

namespace mgis::internal {
  /*!
   * \brief partial specialization to integrate the LocationIdentifier class
   * in MGIS's error handling scheme
   */
  template <>
  struct InvalidValueTraits<::mfem_mgis::LocationIdentifier> {
    //! \brief tag indicating that this class is properly specialized
    static constexpr bool isSpecialized = true;
    //! \brief return an invalid location identifier
    static constexpr auto getValue() noexcept {
      return ::mfem_mgis::LocationIdentifier{};
    }
  };

}  // end of namespace mgis::internal

#include "MFEMMGIS/MeshDiscretization.ixx"

#endif /* LIB_MFEM_MGIS_MESHDISCRETIZATION_HXX */
