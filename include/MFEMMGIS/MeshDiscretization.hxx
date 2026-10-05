/*!
 * \file   include/MFEMMGIS/MeshDiscretization.hxx
 * \brief  This file declares the `MeshDiscretization` class
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
       * \param[in] ids: list of attributes
       */
      AttributesList(const std::vector<size_type>& ids) : attributes(ids) {
        std::sort(this->attributes.begin(), this->attributes.end());
      }  // end of AttributesList
      /*!
       * \brief comparison operator
       * \return if this list is lexicographically lower than rhs
       * \param[in] rhs: right hand side
       */
      [[nodiscard]] bool operator<(const AttributesList& rhs) const noexcept {
        return this->attributes < rhs.attributes;
      }  // end of operator<
      //! \return the list of attributes
      [[nodiscard]] const std::vector<size_type>& getAttributes()
          const noexcept {
        return this->attributes;
      }

     private:
      //! \brief list of attributes
      std::vector<size_type> attributes;
    };
    //! \brief location on which submeshes can be defined
    enum struct Location {
      ON_MATERIALS,  //!< on materials
      ON_BOUNDARIES  //!< on boundaries
    };
    /*!
     * \brief a helper structure used to distinguish the identifiers of a
     * single material (associated with a unique attribute) and identifiers of
     * a single boundary (associated with a unique boundary attribute)
     */
    template <MeshDiscretization::Location>
    struct RawLocationIdentifier {
      //! \brief attribute associated with the material or the boundary
      size_type id;
      /*!
       * \brief comparison operator
       * \return the ordering of the two identifiers
       */
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
      /*!
       * \brief comparison operator
       * \return the ordering of the two identifiers
       */
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
    //! \brief string associated to the `GeneralVerbosityLevel` parameter
    static const char* const GeneralVerbosityLevel;
    //! \return the list of valid parameters
    [[nodiscard]] static std::vector<std::string> getParametersList() noexcept;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - `Parallel` (boolean): if true, a parallel computation is to be
     *    performed. This value is assumed to be false by default.
     * - `MeshFileName` (string): mesh file.
     * - `MeshReadMode` (string): how to read the mesh. Supported values are
     *   "FromScratch" and "Restart". The default value is "FromScratch".
     * - `NumberOfUniformRefinements` (int): number of uniform refinements
     *   applied to the mesh
     * - `Materials` (map): mapping between material names and identifiers
     * - `Boundaries` (map): mapping between boundary names and identifiers
     * - `Points` (map): points to be added to the mesh
     * - `PointsSets` (map): points sets to be added to the mesh
     * - `GeneralVerbosityLevel` (int): with large positive numbers, expect more
     * verbosity
     */
    MeshDiscretization(mgis::Context& ctx, const Parameters& params);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     */
    MeshDiscretization(std::shared_ptr<Mesh<true>> m);
    /*!
     * \brief constructor
     * \param[in] m: mesh
     */
    MeshDiscretization(std::shared_ptr<Mesh<false>> m);
    //! \brief move constructor
    MeshDiscretization(MeshDiscretization&&) noexcept;
    //! \brief copy constructor
    MeshDiscretization(const MeshDiscretization&) noexcept;
    /*!
     * \brief check if a mesh is managed by this mesh discretization
     * \return if the given mesh is managed by this mesh discretization
     * \param[in] m: parallel mesh
     */
    bool manages(const Mesh<true>& m) const noexcept;
    /*!
     * \brief check if a mesh is managed by this mesh discretization
     * \return if the given mesh is managed by this mesh discretization
     * \param[in] m: sequential mesh
     */
    bool manages(const Mesh<false>& m) const noexcept;
    /*!
     * \brief check if a mesh is defined on materials of the main mesh
     * \return if the given mesh is defined on (a subset of) the
     * materials of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this method fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnMaterials(
        Context& ctx, const Mesh<true>& m) const noexcept;
    /*!
     * \brief check if a mesh is defined on materials of the main mesh
     * \return if the given mesh is defined on (a subset of) the
     * materials of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this method fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnMaterials(
        Context& ctx, const Mesh<false>& m) const noexcept;
    /*!
     * \brief check if a mesh is defined on boundaries of the main mesh
     * \return if the given mesh is defined on (a subset of) the
     * boundaries of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this method fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnBoundaries(
        Context& ctx, const Mesh<true>& m) const noexcept;
    /*!
     * \brief check if a mesh is defined on boundaries of the main mesh
     * \return if the given mesh is defined on (a subset of) the
     * boundaries of the main mesh.
     * \param[in, out]  ctx: execution context
     * \param[in]  m: mesh
     *
     * \note this method fails if the given mesh is not managed
     */
    [[nodiscard]] std::optional<bool> isDefinedOnBoundaries(
        Context& ctx, const Mesh<false>& m) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return a pointer to the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    template <bool parallel>
    std::shared_ptr<SubMesh<parallel>> getMutableSubMeshPointer(
        Context& ctx, const Parameter& p, const Location l) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return a pointer to the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    template <bool parallel>
    std::shared_ptr<SubMesh<parallel>> getSubMeshPointer(
        Context& ctx, const Parameter& p, const Location l) noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return a pointer to the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    template <bool parallel>
    std::shared_ptr<const SubMesh<parallel>> getSubMeshPointer(
        Context& ctx, const Parameter& p, const Location l) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    template <bool parallel>
    OptionalReference<SubMesh<parallel>> getMutableSubMeshReference(
        Context& ctx, const Parameter& p, const Location l) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    template <bool parallel>
    OptionalReference<SubMesh<parallel>> getSubMesh(Context& ctx,
                                                    const Parameter& p,
                                                    const Location l) noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    template <bool parallel>
    OptionalReference<const SubMesh<parallel>> getSubMesh(
        Context& ctx, const Parameter& p, const Location l) const noexcept;
    /*!
     * \brief get the shared pointer associated with a managed mesh
     * \return the shared pointer associated with the given mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    template <bool parallel>
    std::shared_ptr<Mesh<parallel>> getMutableMeshPointer(
        Context& ctx, const Mesh<parallel>& m) const noexcept;
    /*!
     * \brief get the shared pointer associated with a managed mesh
     * \return the shared pointer associated with the given mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    template <bool parallel>
    std::shared_ptr<Mesh<parallel>> getMeshPointer(
        Context& ctx, const Mesh<parallel>& m) noexcept;
    /*!
     * \brief get the shared pointer associated with a managed mesh
     * \return the shared pointer associated with the given mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    template <bool parallel>
    std::shared_ptr<const Mesh<parallel>> getMeshPointer(
        Context& ctx, const Mesh<parallel>& m) const noexcept;
#ifdef MFEM_USE_MPI
    /*!
     * \brief convert an attribute of a managed mesh to a location identifier
     * \return the location identifier in the main mesh
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] id: attribute in the mesh
     *
     * \note The rationale behind this method is that submesh may be created on
     * boundaries. In this case, the boundary attributes used to create
     * the submesh become standard attributes of the submesh. This may lead
     * to ambiguity when it comes to determine where a partial quadrature space
     * is defined for instance. The returned location identifier does not have
     * such ambiguity.
     */
    [[nodiscard]] std::optional<LocationIdentifier> getLocationIdentifier(
        Context& ctx, const Mesh<true>& m, const size_type id) const noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief convert an attribute of a managed mesh to a location identifier
     * \return the location identifier in the main mesh
     *
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     * \param[in] id: attribute in the mesh
     *
     * \note The rationale behind this method is that submesh may be created on
     * boundaries. In this case, the boundary attributes used to create
     * the submesh become standard attributes of the submesh. This may lead
     * to ambiguity when it comes to determine where a partial quadrature space
     * is defined for instance. The returned location identifier does not have
     * such ambiguity.
     */
    [[nodiscard]] std::optional<LocationIdentifier> getLocationIdentifier(
        Context& ctx, const Mesh<false>& m, const size_type id) const noexcept;
    /*!
     * \brief set material names
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     */
    [[nodiscard]] bool setMaterialsNames(
        Context& ctx, const std::map<size_type, std::string>& ids) noexcept;
    /*!
     * \brief set boundary names
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     */
    [[nodiscard]] bool setBoundariesNames(
        Context& ctx, const std::map<size_type, std::string>& ids) noexcept;
    /*!
     * \brief get the name of a location
     * \return the name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \param[in, out] ctx: execution context
     * \param[in] id: location identifier
     * \note the method only fails if the identifier is invalid or not defined
     * in the mesh
     */
    [[nodiscard]] std::optional<std::string> getLocationName(
        Context& ctx, const LocationIdentifier& id) const noexcept;
    /*!
     * \brief get the name of a material
     * \return the material name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \param[in, out] ctx: execution context
     * \param[in] id: material identifier
     * \note the method only fails if the material identifier is not defined in
     * the mesh
     */
    [[nodiscard]] std::optional<std::string> getMaterialName(
        Context& ctx, const size_type id) const noexcept;
    /*!
     * \brief get the name of a boundary
     * \return the boundary name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \param[in, out] ctx: execution context
     * \param[in] id: boundary identifier
     * \note the method only fails if the boundary identifier is not defined in
     * the mesh
     */
    [[nodiscard]] std::optional<std::string> getBoundaryName(
        Context& ctx, const size_type id) const noexcept;
    /*!
     * \brief get the identifier of a material
     * \return the material identifier described by the given parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: identifier or name of the material
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] std::optional<size_type> getMaterialIdentifier(
        Context& ctx, const Parameter& p) const noexcept;
    /*!
     * \brief get the identifier of a boundary
     * \return the boundary identifier described by the given parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: identifier or name of the boundary
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] std::optional<size_type> getBoundaryIdentifier(
        Context& ctx, const Parameter& p) const noexcept;
    /*!
     * \brief get the identifiers of a set of materials
     * \return the list of materials identifiers described by the given
     * parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings or integers.
     *
     * Integers are directly interpreted as materials identifiers.
     *
     * Strings are interpreted as regular expressions which allow the selection
     * of materials by names.
     */
    [[nodiscard]] std::optional<std::vector<size_type>> getMaterialsIdentifiers(
        Context& ctx, const Parameter& p) const noexcept;
    /*!
     * \brief get the identifiers of a set of boundaries
     * \return the list of boundaries identifiers described by the given
     * parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter
     *
     * \note The parameter may hold:
     *
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings or integers.
     *
     * Integers are directly interpreted as boundaries identifiers.
     *
     * Strings are interpreted as regular expressions which allow the selection
     * of boundaries by names.
     */
    [[nodiscard]] std::optional<std::vector<size_type>>
    getBoundariesIdentifiers(Context& ctx, const Parameter& p) const noexcept;
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
     * identifiers)
     */
    [[nodiscard]] std::map<size_type, std::string> getMaterialsNames()
        const noexcept;
    /*!
     * \return the names of the boundaries (and their mapping with their
     * identifiers)
     */
    [[nodiscard]] std::map<size_type, std::string> getBoundariesNames()
        const noexcept;
#ifdef MGIS_HAVE_TFEL
    /*!
     * \brief add a point
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \param[in] pt: coordinates of the point
     */
    [[nodiscard]] bool addPoint(Context& ctx,
                                std::string_view n,
                                const Point<2>& pt) noexcept;
    /*!
     * \brief add a point
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \param[in] pt: coordinates of the point
     */
    [[nodiscard]] bool addPoint(Context& ctx,
                                std::string_view n,
                                const Point<3>& pt) noexcept;
    /*!
     * \brief add a points set
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     * \param[in] pts: list of points
     */
    [[nodiscard]] bool addPointsSet(Context& ctx,
                                    std::string_view n,
                                    const std::vector<Point<2>>& pts) noexcept;
    /*!
     * \brief add a points set
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     * \param[in] pts: list of points
     */
    [[nodiscard]] bool addPointsSet(Context& ctx,
                                    std::string_view n,
                                    const std::vector<Point<3>>& pts) noexcept;
    /*!
     * \brief get a point by its name
     * \return the point with the given name
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] std::optional<Point<N>> getPoint(
            Context& ctx, std::string_view n) const noexcept;
    /*!
     * \brief get the registered points
     * \return the registered points
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] OptionalReference<
            const std::map<std::string, Point<N>, std::less<>>>  //
        getPoints(Context& ctx) const noexcept;
    /*!
     * \brief get the registered points sets
     * \return the registered points sets
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] OptionalReference<
            const std::map<std::string, std::vector<Point<N>>, std::less<>>>  //
        getPointsSets(Context& ctx) const noexcept;
    /*!
     * \brief get a points set by its name
     * \return the points set with the given name
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     */
    template <size_type N>
    requires((N == 2) || (N == 3))                                    //
        [[nodiscard]] OptionalReference<const std::vector<Point<N>>>  //
        getPointsSet(Context& ctx, std::string_view n) const noexcept;
#endif /* MGIS_HAVE_TFEL */

    //! \brief destructor
    ~MeshDiscretization();

   protected:
    // friend functions and operators
    /*!
     * \brief display information about a mesh discretization
     * \param[in, out] ctx: execution context
     * \param[out] os: output stream
     * \param[in] m: mesh discretization
     * \return true on success
     */
    friend bool getInformation<MeshDiscretization>(
        Context& ctx, std::ostream& os, const MeshDiscretization& m) noexcept;
    friend bool operator==(const MeshDiscretization& lhs,
                           const MeshDiscretization& rhs) noexcept;
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
     * \brief get the shared pointer associated with a managed mesh
     * \return a mutable pointer to the given parallel mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    [[nodiscard]] std::shared_ptr<Mesh<true>> getMutableParallelMeshPointer(
        Context& ctx, const Mesh<true>& m) const noexcept;
    /*!
     * \brief get the shared pointer associated with a managed mesh
     * \return a mutable pointer to the given sequential mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    [[nodiscard]] std::shared_ptr<Mesh<false>> getMutableSequentialMeshPointer(
        Context& ctx, const Mesh<false>& m) const noexcept;
    /*!
     * \brief get the shared pointer associated with a managed mesh
     * \return a pointer to the given parallel mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    [[nodiscard]] std::shared_ptr<const Mesh<true>> getParallelMeshPointer(
        Context& ctx, const Mesh<true>& m) const noexcept;
    /*!
     * \brief get the shared pointer associated with a managed mesh
     * \return a pointer to the given sequential mesh, if managed by
     * this mesh discretization
     * \param[in, out] ctx: execution context
     * \param[in] m: mesh
     */
    [[nodiscard]] std::shared_ptr<const Mesh<false>> getSequentialMeshPointer(
        Context& ctx, const Mesh<false>& m) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return the parallel sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    std::shared_ptr<const SubMesh<true>> getParallelSubMeshPointer(
        Context& ctx, const Parameter& p, const Location l) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return the parallel sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    std::shared_ptr<SubMesh<true>> getParallelMutableSubMeshPointer(
        Context& ctx, const Parameter& p, const Location l) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return the sequential sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    std::shared_ptr<const SubMesh<false>> getSequentialSubMeshPointer(
        Context& ctx, const Parameter& p, const Location l) const noexcept;
    /*!
     * \brief get or create the sub mesh associated with the given ids
     * \return the sequential sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     * \param[in] l: location on which the submesh is built (materials or
     * boundaries)
     *
     * \note the parameter may contain an integer, a string, a vector of
     * parameters which are either strings or integers.
     */
    std::shared_ptr<SubMesh<false>> getSequentialMutableSubMeshPointer(
        Context& ctx, const Parameter& p, const Location l) const noexcept;

#ifdef MGIS_HAVE_TFEL
    /*!
     * \brief get the registered points in 2D
     * \return the registered points in 2D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, Point<2>, std::less<>>>
    getPoints2D(Context& ctx) const noexcept;
    /*!
     * \brief get the registered points in 3D
     * \return the registered points in 3D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, Point<3>, std::less<>>>
    getPoints3D(Context& ctx) const noexcept;
    /*!
     * \brief get the registered points sets in 2D
     * \return the registered points sets in 2D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, std::vector<Point<2>>, std::less<>>>
    getPointsSets2D(Context& ctx) const noexcept;
    /*!
     * \brief get the registered points sets in 3D
     * \return the registered points sets in 3D
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, std::vector<Point<3>>, std::less<>>>
    getPointsSets3D(Context& ctx) const noexcept;
    /*!
     * \brief get a point by its name
     * \return the point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     */
    [[nodiscard]] std::optional<Point<2>> getPoint2D(
        Context& ctx, std::string_view n) const noexcept;
    /*!
     * \brief get a point by its name
     * \return the point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     */
    [[nodiscard]] std::optional<Point<3>> getPoint3D(
        Context& ctx, std::string_view n) const noexcept;
    /*!
     * \brief get a points set by its name
     * \return the points set with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     */
    [[nodiscard]] OptionalReference<const std::vector<Point<2>>> getPointsSet2D(
        Context& ctx, std::string_view n) const noexcept;
    /*!
     * \brief get a points set by its name
     * \return the points set with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     */
    [[nodiscard]] OptionalReference<const std::vector<Point<3>>> getPointsSet3D(
        Context& ctx, std::string_view n) const noexcept;
#endif /* MGIS_HAVE_TFEL */
    //! \brief internal structure to implement the pimpl idiom
    struct Implementation;
    //! \brief pointer to the internal implementation
    std::shared_ptr<Implementation> pimpl;
  };  // end of MeshDiscretization

  /*!
   * \brief a simple alias for the identifier of a single
   * material (associated with a unique attribute)
   */
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
   * \brief check if a location identifier is invalid
   * \return if the given location identifier is invalid
   * \param[in] l: location identifier
   */
  [[nodiscard]] constexpr bool isInvalid(const LocationIdentifier& l) noexcept {
    const auto mok = isValid(l.material_identifier);
    const auto bok = isValid(l.boundary_identifier);
    const auto b1 = (!mok) && (!bok);  // none is valid
    const auto b2 = mok && bok;        // both are valid
    return b1 || b2;
  }  // end of isInvalid
  /*!
   * \return a description of the location
   * \param[in] l: location identifier
   *
   * Examples of the return values are:
   *
   * - material '1'
   * - bounadary '1'
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::string getLocationDescription(
      const LocationIdentifier&);
  /*!
   * \brief compare two mesh discretisations to see if they point to the same
   * underlying implementation
   *
   * \return true if both discretisations share the same implementation
   * \param[in] lhs: left hand side
   * \param[in] rhs: right hand side
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool operator==(
      const MeshDiscretization& lhs, const MeshDiscretization& rhs) noexcept;
  /*!
   * \brief compare two mesh discretisations
   *
   * \return true if the discretisations do not share the same implementation
   * \param[in] lhs: left hand side
   * \param[in] rhs: right hand side
   *
   * \note material names and boundary names may differ in both
   * discretisations.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool operator!=(
      const MeshDiscretization& lhs, const MeshDiscretization& rhs) noexcept;
  /*!
   * \brief get the space dimension of the mesh
   * \return the space dimension
   * \param[in] m: mesh discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getSpaceDimension(
      const MeshDiscretization& m) noexcept;
  /*!
   * \brief get the attributes of the materials
   * \return the list of materials attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getMaterialsAttributes(const MeshDiscretization& m) noexcept;
  /*!
   * \brief get the attributes of the boundaries
   * \return the list of boundaries attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getBoundariesAttributes(const MeshDiscretization& m) noexcept;

  /*!
   * \brief get the identifiers of a set of materials
   * \return the list of materials identifiers described by the given
   * parameter.
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] m: mesh discretization
   * \param[in] p: parameter
   *
   * \note The parameter may hold:
   *
   * - an integer
   * - a string
   * - a vector of parameters which must be either strings or integers.
   *
   * Integers are directly interpreted as materials identifiers.
   *
   * Strings are interpreted as regular expressions which allow the selection
   * of materials by names.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<size_type> getMaterialsIdentifiers(
      attributes::Throwing throwing,
      const MeshDiscretization& m,
      const Parameter& p);
  /*!
   * \brief get the identifiers of a set of boundaries
   * \return the list of boundaries identifiers described by the given
   * parameter.
   * \param[in] throwing: dummy attribute to indicate that this function may
   * throw an exception
   * \param[in] m: mesh discretization
   * \param[in] p: parameter
   *
   * \note The parameter may hold:
   *
   * - an integer
   * - a string
   * - a vector of parameters which must be either strings or integers.
   *
   * Integers are directly interpreted as boundaries identifiers.
   *
   * Strings are interpreted as regular expressions which allow the selection
   * of boundaries by names.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::vector<size_type>
  getBoundariesIdentifiers(attributes::Throwing throwing,
                           const MeshDiscretization& m,
                           const Parameter& p);

  /*!
   * \brief display information about a mesh discretization
   *
   * \param[in, out] ctx: execution context
   * \param[out] os: output stream
   * \param[in] m: mesh discretization
   * \return true on success
   */
  template <>
  MFEM_MGIS_EXPORT bool getInformation<MeshDiscretization>(
      Context& ctx, std::ostream& os, const MeshDiscretization& m) noexcept;

#ifdef MFEM_USE_MPI

  /*!
   * \brief get the MPI communicator of a mesh discretization
   * \return the MPI communicator associated with the mesh discretization
   * \param[in] m: mesh discretization
   *
   * \note If a sequential computation is described, `MPI_COMM_WORLD` is
   * returned.
   */
  MFEM_MGIS_EXPORT [[nodiscard]] MPI_Comm getMPICommunicator(
      const MeshDiscretization& m) noexcept;

#endif /* MFEM_USE_MPI */

  /*!
   * \brief check if the current process is the main one
   * \return if the current process is the main one (the process of rank 0)
   * \param[in] m: mesh discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isMainProcess(
      const MeshDiscretization& m) noexcept;

  /*!
   * \brief check a location identifier against a mesh discretization
   * \return if the given location identifier is consistent with the mesh
   * discretization
   *
   * \param[in, out] ctx: execution context
   * \param[in] m: mesh discretization
   * \param[in] l: location identifier
   *
   * This check fails if:
   *
   * - the identifier is invalid
   * - the material identifier (if valid) is not a mesh attribute
   * - the boundary identifier (if valid) is not a boundary mesh attribute
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool check(
      Context& ctx,
      const MeshDiscretization& m,
      const LocationIdentifier& l) noexcept;

#ifdef MGIS_HAVE_TFEL

  /*!
   * \brief create a point from a parameter
   * \return the point defined by the given parameter, empty on failure.
   * The parameter may be the name of a point of the mesh.
   * \param[in, out] ctx: execution context
   * \param[in] m: mesh discretization
   * \param[in] p: name of a point or coordinates of the point
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<Point<N>> makePoint(
          Context& ctx,
          const MeshDiscretization& m,
          const Parameter& p) noexcept;

  /*!
   * \brief create a points set from a parameter
   * \return the points set defined by the given parameter, empty on failure.
   * The parameter may be the name of a points set of the mesh.
   * The points may be given by their names.
   * \param[in, out] ctx: execution context
   * \param[in] m: mesh discretization
   * \param[in] p: name of a points set, list of points or parameters defining
   * a curve
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsSet(
          Context& ctx,
          const MeshDiscretization& m,
          const Parameter& p) noexcept;

  /*!
   * \brief discretize a curve
   * \return the points of the curve defined by the given parameters, empty on
   * failure.
   * The points defining the curve may be given by their names.
   * \param[in, out] ctx: execution context
   * \param[in] m: mesh discretization
   * \param[in] p: parameters defining the curve
   */
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      [[nodiscard]] std::optional<std::vector<Point<N>>> makePointsOnCurve(
          Context& ctx,
          const MeshDiscretization& m,
          const Parameters& p) noexcept;

  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<2>> makePoint<2>(
      Context&, const MeshDiscretization&, const Parameter&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<Point<3>> makePoint<3>(
      Context&, const MeshDiscretization&, const Parameter&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsSet<2>(Context&,
                   const MeshDiscretization&,
                   const Parameter&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsSet<3>(Context&,
                   const MeshDiscretization&,
                   const Parameter&) noexcept;
  //! \brief 2D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<2>>>
  makePointsOnCurve<2>(Context&,
                       const MeshDiscretization&,
                       const Parameters&) noexcept;
  //! \brief 3D specialisation
  template <>
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::vector<Point<3>>>
  makePointsOnCurve<3>(Context&,
                       const MeshDiscretization&,
                       const Parameters&) noexcept;

#endif /* MGIS_HAVE_TFEL */

}  // end of namespace mfem_mgis

namespace mgis::internal {
  /*!
   * \brief specialization to integrate the LocationIdentifier class
   * in MGIS's error handling scheme
   */
  template <>
  struct InvalidValueTraits<::mfem_mgis::LocationIdentifier> {
    //! \brief tag indicating that this class is properly specialized
    static constexpr bool isSpecialized = true;
    //! \return an invalid location identifier
    static constexpr auto getValue() noexcept {
      return ::mfem_mgis::LocationIdentifier{};
    }
  };

}  // end of namespace mgis::internal

#include "MFEMMGIS/MeshDiscretization.ixx"

#endif /* LIB_MFEM_MGIS_MESHDISCRETIZATION_HXX */
