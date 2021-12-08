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
#include <memory>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameter;
  struct Parameters;

  /*!
   * \brief a simple class used to:
   * - handle the life time of the mesh and the finite element collection.
   * - create and handle a finite element space
   */
  struct MFEM_MGIS_EXPORT FiniteElementDiscretization {
    //! \brief string associated to the `Parallel` parameter
    static const char* const Parallel;
    //! \brief string associated to the `MeshFileName` parameter
    static const char* const MeshFileName;
    //! \brief string associated to the `FiniteElementFamily` parameter
    static const char* const FiniteElementFamily;
    //! \brief string associated to the `FiniteElementOrder` parameter
    static const char* const FiniteElementOrder;
    //! \brief string associated to the `UnknownsSize` parameter
    static const char* const UnknownsSize;
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
     * - `Parallel` (boolean): if true, a parallel computation is to be be
     *    performed. This value if assumed to be false by default.
     * - `MeshFileName` (string): mesh file.
     * - `FiniteElementFamily` (string): name of the finite element family to be
     *   used. Supported families are:
     *   - `H1`:
     *   The default value if `H1`.
     * - `FiniteElementOrder` (int): order of the polynomial approximation.
     * - `UnknownsSize` (int): number of components of the unknows
     * - `NumberOfUniformRefinements` (int): number of uniform refinements
     *   applied to the mesh
     * - `GeneralVerbosityLevel` (int): with large positive numbers, expect more
     * verbosity
     */
    FiniteElementDiscretization(const Parameters&);
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
    /*!
     * \brief set material names
     * \param[in] ids: mapping between mesh identifiers and names
     */
    void setMaterialsNames(const std::map<size_type, std::string>&);
    /*!
     * \brief set material names
     * \param[in] ids: mapping between mesh identifiers and names
     */
    void setBoundariesNames(const std::map<size_type, std::string>&);
    /*!
     * \return the material identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    size_type getMaterialIdentifier(const Parameter&) const;
    /*!
     * \return the material identifier by the given parameter.
     * \note The parameter may hold an integer or a string.
     */
    size_type getBoundaryIdentifier(const Parameter&) const;
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
    std::vector<size_type> getMaterialsIdentifiers(const Parameter&) const;
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
    std::vector<size_type> getBoundariesIdentifiers(const Parameter&) const;
    /*!
     * \brief set material names
     * \param[in] m: material names
     */
    //! \return the mesh
    template <bool parallel>
    Mesh<parallel>& getMesh();
    //! \return the mesh
    template <bool parallel>
    const Mesh<parallel>& getMesh() const;
    //! \return the finite element space
    template <bool parallel>
    FiniteElementSpace<parallel>& getFiniteElementSpace();
    //! \return the finite element space
    template <bool parallel>
    const FiniteElementSpace<parallel>& getFiniteElementSpace() const;
    //! \return the finite element collection
    const FiniteElementCollection& getFiniteElementCollection() const;
    //! \return if this object is built to run parallel computations
    bool describesAParallelComputation() const;
    //! \brief destructor
    ~FiniteElementDiscretization();

   private:
    //! \brief mesh
#ifdef MFEM_USE_MPI
    std::shared_ptr<Mesh<true>> parallel_mesh;
#endif /* MFEM_USE_MPI */
    std::shared_ptr<Mesh<false>> sequential_mesh;
    //! \brief finite element collection
    std::shared_ptr<const FiniteElementCollection> fec;
    //! \brief finite element space
#ifdef MFEM_USE_MPI
    std::unique_ptr<FiniteElementSpace<true>> parallel_fe_space;
#endif /* MFEM_USE_MPI */
    std::unique_ptr<FiniteElementSpace<false>> sequential_fe_space;
    //! \brief mapping between materials identifiers and names
    std::map<size_type, std::string> materials_names;
    //! \brief mapping between materials boundaries and names
    std::map<size_type, std::string> boundaries_names;
  };  // end of FiniteElementDiscretization

  /*!
   * \brief return the space dimension
   * \param[in] fed: finite element discretization
   */
  MFEM_MGIS_EXPORT size_type getSpaceDimension(const FiniteElementDiscretization&);

  /*!
   * \brief return the total number of unknowns
   * \param[in] fed: finite element discretization
   */
  MFEM_MGIS_EXPORT size_type getTrueVSize(const FiniteElementDiscretization&);

  /*!
   * \brief return the list of materials attributes
   * \param[in] fed: finite element discretisation
   */
  MFEM_MGIS_EXPORT const mfem::Array<size_type>& getMaterialsAttributes(
      const FiniteElementDiscretization&);
  /*!
   * \brief return the list of boundaries attributes
   * \param[in] fed: finite element discretisation
   */
  MFEM_MGIS_EXPORT const mfem::Array<size_type>& getBoundariesAttributes(
      const FiniteElementDiscretization&);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/FiniteElementDiscretization.ixx"

#endif /* LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX */
