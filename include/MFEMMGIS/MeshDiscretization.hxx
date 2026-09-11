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

  // forward declaration
  struct Parameter;
  struct Parameters;

  //! \brief a simple class used to handle the life time of the mesh
  struct MFEM_MGIS_EXPORT [[nodiscard]] MeshDiscretization {
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
     * \brief return the material name associated with the given identifier
     * \param[in, out] ctx: execution context
     * \param[in] id: material identifier
     * \return the material name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \note the method only fails if the material identifier is not defined in
     * the mesh
     */
    [[nodiscard]] std::optional<std::string> getMaterialName(
        Context&, const size_type) const noexcept;
    /*!
     * \brief return the boundary name associated with the given identifier
     * \param[in, out] ctx: execution context
     * \param[in] id: boundary identifier
     * \return the boundary name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
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
     * \brief return the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     * \return the mesh
     */
    template <bool parallel>
    [[nodiscard]] Mesh<parallel>& getMesh() noexcept;
    /*!
     * \brief return the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     * \return the mesh
     */
    template <bool parallel>
    [[nodiscard]] const Mesh<parallel>& getMesh() const noexcept;
    /*!
     * \brief return a pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     * \return a pointer to the mesh
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<Mesh<parallel>> getMeshPointer() noexcept;
    /*!
     * \brief return a pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     * \return a pointer to the mesh
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<const Mesh<parallel>> getMeshPointer()
        const noexcept;
    /*!
     * \brief return a mutable pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     * \return a mutable pointer to the mesh
     */
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<Mesh<parallel>> getMutableMeshPointer()
        const noexcept;
    //! \return if this object is built to run parallel computations
    [[nodiscard]] bool describesAParallelComputation() const noexcept;
    /*!
     * \brief return the names of the materials (and their mapping with their
     * identifiers
     */
    [[nodiscard]] std::map<size_type, std::string> getMaterialsNames()
        const noexcept;
    /*!
     * \brief return the names of the boundaries (and their mapping with their
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
     * \brief return the point with the given name
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \return the point with the given name
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] std::optional<Point<N>> getPoint(
            Context&, std::string_view) const noexcept;
    /*!
     * \brief return the registered points
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \return the registered points
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] OptionalReference<
            const std::map<std::string, Point<N>, std::less<>>>  //
        getPoints(Context&) const noexcept;
    /*!
     * \brief return the registered points sets
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \return the registered points sets
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] OptionalReference<
            const std::map<std::string, std::vector<Point<N>>, std::less<>>>  //
        getPointsSets(Context&) const noexcept;
    /*!
     * \brief return the registered set of points
     * \tparam N: space dimension (2 or 3)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     * \return the points set with the given name
     */
    template <size_type N>
    requires((N == 2) || (N == 3))                                    //
        [[nodiscard]] OptionalReference<const std::vector<Point<N>>>  //
        getPointsSet(Context&, std::string_view) const noexcept;
#endif /* MGIS_HAVE_TFEL */

    /*!
     * \brief return the sub mesh associated with the given ids
     * \tparam parallel: whether to get the parallel sub mesh or not
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    template <bool parallel>
    std::shared_ptr<SubMesh<parallel>> getSubMesh(
        Context&, const Parameter&) const noexcept;

    //! \brief destructor
    ~MeshDiscretization();

   protected:
    /*!
     * \brief return a mutable pointer to the underlying parallel mesh
     * \return a mutable pointer to the parallel mesh
     */
    [[nodiscard]] std::shared_ptr<Mesh<true>> getMutableParallelMeshPointer()
        const noexcept;
    /*!
     * \brief return a mutable pointer to the underlying sequential mesh
     * \return a mutable pointer to the sequential mesh
     */
    [[nodiscard]] std::shared_ptr<Mesh<false>> getMutableSequentialMeshPointer()
        const noexcept;
    /*!
     * \brief return a pointer to the underlying parallel mesh
     * \return a pointer to the parallel mesh
     */
    [[nodiscard]] std::shared_ptr<const Mesh<true>> getParallelMeshPointer()
        const noexcept;
    /*!
     * \brief return a pointer to the underlying sequential mesh
     * \return a pointer to the sequential mesh
     */
    [[nodiscard]] std::shared_ptr<const Mesh<false>> getSequentialMeshPointer()
        const noexcept;
    /*!
     * \brief return the parallel sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    std::shared_ptr<SubMesh<true>> getParallelSubMesh(
        Context&, const Parameter&) const noexcept;
    /*!
     * \brief return the sequential sub mesh associated with the given ids
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter containing the list of ids
     *
     * \note the parameter may contain a integer, a string, a vector of
     * parameters which are either string or integers.
     */
    std::shared_ptr<SubMesh<false>> getSequentialSubMesh(
        Context&, const Parameter&) const noexcept;

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

  /*!
   * \brief compare two mesh discretisations to see if they point to the same
   * underlying meshes
   *
   * \param[in] lhs: left hand side
   * \param[in] rhs: right hand side
   *
   * \note material names and boundary names may different in both
   * discretisations.
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
   * \brief return the space dimension
   * \param[in] m: mesh discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] size_type getSpaceDimension(
      const MeshDiscretization&) noexcept;
  /*!
   * \brief return the list of materials attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getMaterialsAttributes(const MeshDiscretization&) noexcept;
  /*!
   * \brief return the list of boundaries attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getBoundariesAttributes(const MeshDiscretization&) noexcept;

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
  /*!
   * \return if the current process is the main one (the process of rank 0)
   * \param[in] m: mesh discretization
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool isMainProcess(
      const MeshDiscretization&) noexcept;

#endif /* MFEM_USE_MPI */

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

#include "MFEMMGIS/MeshDiscretization.ixx"

#endif /* LIB_MFEM_MGIS_MESHDISCRETIZATION_HXX */
