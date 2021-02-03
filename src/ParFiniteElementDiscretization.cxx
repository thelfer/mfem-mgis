/*!
 * \file   src/ParFiniteElementDiscretization.cxx
 * \brief
 * \author Guillaume Latu
 * \date 02/02/2021
 */

#include <utility>
#include <mfem/fem/pfespace.hpp>
#include "MGIS/Raise.hxx"

#ifdef MFEM_USE_MPI
#include "MFEMMGIS/ParFiniteElementDiscretization.hxx"

namespace mfem_mgis {

  ParFiniteElementDiscretization::ParFiniteElementDiscretization(
      std::shared_ptr<mfem::ParMesh> m,
      std::shared_ptr<const mfem::FiniteElementCollection> c,
      const size_type d)
      : pmesh(std::move(m)), fec(std::move(c)) {
    this->pfe_space = std::make_unique<mfem::ParFiniteElementSpace>(
        this->pmesh.get(), this->fec.get(), d);
    mfem::Mesh *parentmesh=this->pmesh.get();
    this->fed = std::make_unique<FiniteElementDiscretization>(
      std::make_shared<mfem::Mesh>(*parentmesh), 
     this->fec, d);       
  }  // end of ParFiniteElementDiscretization

  ParFiniteElementDiscretization::ParFiniteElementDiscretization(
      std::shared_ptr<mfem::ParMesh> m,
      std::shared_ptr<const mfem::FiniteElementCollection> c,
      std::unique_ptr<mfem::ParFiniteElementSpace> s)
      : pmesh(std::move(m)), fec(std::move(c)), pfe_space(std::move(s)) {
    this->fed = std::make_unique<FiniteElementDiscretization>(
      std::make_shared<mfem::Mesh>(*this->pmesh.get()), 
      this->fec,
      std::make_unique<mfem::FiniteElementSpace>(*this->pfe_space.get()));       
    if (this->pmesh.get() != this->pfe_space->GetMesh()) {
      mgis::raise(
          "ParFiniteElementDiscretization::ParFiniteElementDiscretization: "
          "mesh pointer don't match the mesh on which the finite element space "
          "is built");
    }
  }  // end of ParFiniteElementDiscretization

  mfem::ParMesh& ParFiniteElementDiscretization::getMesh() {
    return *(this->pmesh);
  }  // end of getMesh

  const mfem::ParMesh& ParFiniteElementDiscretization::getMesh() const {
    return *(this->pmesh);
  }  // end of getMesh

  mfem::ParFiniteElementSpace&
  ParFiniteElementDiscretization::getFiniteElementSpace() {
    return *(this->pfe_space);
  }  // end of getFiniteElementSpace

  const mfem::ParFiniteElementSpace&
  ParFiniteElementDiscretization::getFiniteElementSpace() const {
    return *(this->pfe_space);
  }  // end of getFiniteElementSpace

  const mfem::FiniteElementCollection&
  ParFiniteElementDiscretization::getFiniteElementCollection() const {
    return *(this->fec);
  }  // end of getFiniteElementCollection

  std::shared_ptr<FiniteElementDiscretization> 
  ParFiniteElementDiscretization::getFiniteElementDiscretization() {
    return (this->fed);
  }  // end of getFiniteElementDiscretization

  ParFiniteElementDiscretization::~ParFiniteElementDiscretization() = default;

} // end of namespace mfem_mgis

#endif /* MFEM_USE_MPI */
