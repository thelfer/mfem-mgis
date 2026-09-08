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
    //! \brief string associated to the `NumberOfUniformRefinements` parameter
    static const char* const NumberOfUniformRefinements;
    //! \brief string associated to the `VerbosityLevel` parameter
    static const char* const GeneralVerbosityLevel;
    //!
    [[noreturn]] static void reportInvalidParallelMesh();
    //!
    [[noreturn]] static void reportInvalidSequentialMesh();
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
     * \return the material name associated with the given identifier, if it is
     * defined. If not defined, an empty string is returned
     *
     * \param[in, out] ctx: execution context
     * \param[in] id: material identifier
     *
     * \note the method only fails if the material identifier is not defined in
     * the mesh
     */
    [[nodiscard]] std::optional<std::string> getMaterialName(
        Context&, const size_type) const noexcept;
    /*!
     * \return the boundary name associated with the given identifier, if it is
     * defined. If not defined, an empty string is returned
     *
     * \param[in, out] ctx: execution context
     * \param[in] id: boundary identifier
     *
     * \note the method only fails is the boundary identifier is defined in the
     * mesh
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
     * \return the material identifier by the given parameter.
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
    //! \return the mesh
    template <bool parallel>
    [[nodiscard]] Mesh<parallel>& getMesh();
    //! \return the mesh
    template <bool parallel>
    [[nodiscard]] const Mesh<parallel>& getMesh() const;
    //! \return the mesh
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<Mesh<parallel>> getMeshPointer();
    //! \return the mesh
    template <bool parallel>
    [[nodiscard]] std::shared_ptr<const Mesh<parallel>> getMeshPointer() const;
    //! \return if this object is built to run parallel computations
    [[nodiscard]] bool describesAParallelComputation() const;
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
     * \return the point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] std::optional<Point<N>> getPoint(
            Context&, std::string_view) const noexcept;
    /*!
     * \return the registred points
     * \param[in, out] ctx: execution context
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] std::optional<std::map<std::string, Point<N>>> getPoints(
            Context&) const noexcept;
    /*!
     * \return the registred set of points
     * \param[in, out] ctx: execution context
     */
    template <size_type N>
    requires((N == 2) || (N == 3))  //
        [[nodiscard]] std::optional<std::vector<Point<N>>> getPointsSet(
            Context&, std::string_view) const noexcept;
#endif /* MGIS_HAVE_TFEL */

    //! \brief destructor
    ~MeshDiscretization();

   protected:
    /*!
     * \brief set material names
     * \param[in] ids: mapping between mesh identifiers and names
     */
    [[deprecated]] void setMaterialsNames(
        const std::map<size_type, std::string>&);
    /*!
     * \brief set material names
     * \param[in] ids: mapping between mesh identifiers and names
     */
    [[deprecated]] void setBoundariesNames(
        const std::map<size_type, std::string>&);
    //
    /*!
     * \return the material identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    [[deprecated, nodiscard]] size_type getMaterialIdentifier(
        const Parameter&) const;
    /*!
     * \return the material identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    [[deprecated, nodiscard]] size_type getBoundaryIdentifier(
        const Parameter&) const;
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
    [[deprecated, nodiscard]] std::vector<size_type> getMaterialsIdentifiers(
        const Parameter&) const;
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
    [[deprecated, nodiscard]] std::vector<size_type> getBoundariesIdentifiers(
        const Parameter&) const;
    /*!
     * \brief set names of materials
     * \param[in] 1: dummy parameter indicated that this function may throw
     * \param[in] ids: mapping between mesh identifiers and names
     */
    void setMaterialsNames(attributes::Throwing,
                           const std::map<size_type, std::string>&);
    /*!
     * \brief set names of boundaries
     * \param[in] 1: dummy parameter indicated that this function may throw
     * \param[in] ids: mapping between mesh identifiers and names
     */
    void setBoundariesNames(attributes::Throwing,
                            const std::map<size_type, std::string>&);
#ifdef MGIS_HAVE_TFEL
    /*!
     * \return the registred points
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] std::optional<std::map<std::string, Point<2>>> getPoints2D(
        Context&) const noexcept;
    /*!
     * \return the registred points
     * \param[in, out] ctx: execution context
     */
    [[nodiscard]] std::optional<std::map<std::string, Point<3>>> getPoints3D(
        Context&) const noexcept;
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
    [[nodiscard]] std::optional<std::vector<Point<2>>> getPointsSet2D(
        Context&, std::string_view) const noexcept;
    /*!
     * \return the point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     */
    [[nodiscard]] std::optional<std::vector<Point<3>>> getPointsSet3D(
        Context&, std::string_view) const noexcept;
#endif /* MGIS_HAVE_TFEL */

#ifdef MFEM_USE_MPI
    //! \brief parallel mesh
    std::shared_ptr<Mesh<true>> parallel_mesh;
#endif /* MFEM_USE_MPI */
    //! \brief sequential mesh
    std::shared_ptr<Mesh<false>> sequential_mesh;
    //! \brief mapping between materials identifiers and names
    std::map<size_type, std::string> materials_names;
    //! \brief mapping between materials boundaries and names
    std::map<size_type, std::string> boundaries_names;
#ifdef MGIS_HAVE_TFEL
    //! \brief points declared by the user, only valid for a 2D mesh
    std::map<std::string, Point<2>, std::less<>> points2D;
    //! \brief points declared by the user, only valid for a 3D mesh
    std::map<std::string, Point<3>, std::less<>> points3D;
    //! \brief points set declared by the user, only valid for a 2D mesh
    std::map<std::string, std::vector<Point<2>>, std::less<>> pointsSet2D;
    //! \brief points set declared by the user, only valid for a 3D mesh
    std::map<std::string, std::vector<Point<3>>, std::less<>> pointsSet3D;
#endif /* MGIS_HAVE_TFEL */
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
      const MeshDiscretization&);
  /*!
   * \brief return the list of materials attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getMaterialsAttributes(const MeshDiscretization&);
  /*!
   * \brief return the list of boundaries attributes
   * \param[in] m: mesh discretisation
   */
  MFEM_MGIS_EXPORT [[nodiscard]] const mfem::Array<size_type>&
  getBoundariesAttributes(const MeshDiscretization&);

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
