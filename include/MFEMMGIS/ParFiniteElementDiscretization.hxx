/*!
 * \file   include/MFEMMGIS/FiniteElementDiscretization.hxx
 * \brief
 * \author Guillaume Latu
 * \date 02/02/2021
 */

#ifndef LIB_MFEM_MGIS_PARFINITEELEMENTDISCRETIZATION_HXX
#define LIB_MFEM_MGIS_PARFINITEELEMENTDISCRETIZATION_HXX

#ifdef MFEM_USE_MPI

#include <memory>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

  /*!
   * \brief a simple class used to:
   * - handle the life time of the mesh and the finite element collection.
   * - create and handle a parallel finite element space
   */
  struct MFEM_MGIS_EXPORT ParFiniteElementDiscretization {
    /*!
     * \brief constructor
     * \param[in] m: parallel mesh
     * \param[in] c: collection
     * \param[in] d: size of the unknowns
     *
     * \note this methods creates the finite element space.
     */
    ParFiniteElementDiscretization(
        std::shared_ptr<mfem::ParMesh>,
        std::shared_ptr<const mfem::FiniteElementCollection>,
        const size_type);
    /*!
     * \brief constructor
     * \param[in] m: parallel mesh
     * \param[in] c: collection
     * \param[in] s: parallel finite element space
     */
    ParFiniteElementDiscretization(
        std::shared_ptr<mfem::ParMesh>,
        std::shared_ptr<const mfem::FiniteElementCollection>,
        std::unique_ptr<mfem::ParFiniteElementSpace>);
    //! \return the mesh
    mfem::ParMesh& getMesh();
    //! \return the mesh
    const mfem::ParMesh& getMesh() const;
    //! \return the finite element space
    mfem::ParFiniteElementSpace& getFiniteElementSpace();
    //! \return the finite element space
    const mfem::ParFiniteElementSpace& getFiniteElementSpace() const;
    //! \return the finite element collection
    const mfem::FiniteElementCollection& getFiniteElementCollection() const;
    //! \return a Finite Element Discretization
    std::shared_ptr<FiniteElementDiscretization> getFiniteElementDiscretization();
    //! \brief destructor
    ~ParFiniteElementDiscretization();

   private:
    //! \brief mesh
    std::shared_ptr<mfem::ParMesh> pmesh;
    //! \brief finite element collection
    std::shared_ptr<const mfem::FiniteElementCollection> fec;
    //! \brief finite element space
    std::unique_ptr<mfem::ParFiniteElementSpace> pfe_space;
    //! \brief finite element discretization
    std::shared_ptr<FiniteElementDiscretization> fed;
  };  // end of ParFiniteElementDiscretization

}  // end of namespace mfem_mgis

#endif /* MFEM_USE_MPI */

#endif /* LIB_MFEM_MGIS_PARFINITEELEMENTDISCRETIZATION_HXX */
