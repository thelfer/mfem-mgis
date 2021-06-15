/*!
 * \file   include/MFEMMGIS/FiniteElementDiscretization.hxx
 * \brief
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#ifndef LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX
#define LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX

#include <memory>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
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
     * - `GeneralVerbosityLevel` (int): with large positive numbers, expect more verbosity
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
  };  // end of FiniteElementDiscretization

  /*!
   * \brief return the total number of unknowns
   * \param[in] fed: finite element discretization
   */
  MFEM_MGIS_EXPORT size_type getTrueVSize(const FiniteElementDiscretization&);


  /*!
   * \brief Extract the file extension
   * \param[in] s: string corresponding to a file name
   */
  static std::string getFileExt(const std::string& s); 

  
  /*!
   * \brief load a mesh (sequential)
   * \param[in] s: string corresponding to a file name
   *
   * \note MED format is handled in addition to standard MFEM 
   *       input formats.
   */
  MFEM_MGIS_EXPORT std::shared_ptr<Mesh<false>> loadMeshSequential(
       const std::string& mesh_name,
       int generate_edges = 0,
       int refine = 1,
       bool fix_orientation = true);

  
}  // end of namespace mfem_mgis

#include "MFEMMGIS/FiniteElementDiscretization.ixx"

#endif /* LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_HXX */
