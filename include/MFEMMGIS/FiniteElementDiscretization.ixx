/*!
 * \file   include/MFEMMGIS/FiniteElementDiscretization.ixx
 * \brief
 * \author Thomas Helfer
 * \date   13/02/2021
 */

#ifndef LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_IXX
#define LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_IXX

namespace mfem_mgis {

  template <bool parallel>
  inline Mesh<parallel>& FiniteElementDiscretization::getMesh() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        FiniteElementDiscretization::reportInvalidParallelMesh();
      }
      return *(this->parallel_mesh);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        FiniteElementDiscretization::reportInvalidSequentialMesh();
      }
      return *(this->sequential_mesh);
    }
  }  // end of getMesh

  template <bool parallel>
  inline const Mesh<parallel>& FiniteElementDiscretization::getMesh() const {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        FiniteElementDiscretization::reportInvalidParallelMesh();
      }
      return *(this->parallel_mesh);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        FiniteElementDiscretization::reportInvalidSequentialMesh();
      }
      return *(this->sequential_mesh);
    }
  }  // end of getMesh

  template <bool parallel>
  std::shared_ptr<Mesh<parallel>>
  FiniteElementDiscretization::getMeshPointer() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        FiniteElementDiscretization::reportInvalidParallelMesh();
      }
      return this->parallel_mesh;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        FiniteElementDiscretization::reportInvalidSequentialMesh();
      }
      return this->sequential_mesh;
    }
  }  // end of getMeshPointer

  template <bool parallel>
  std::shared_ptr<const Mesh<parallel>>
  FiniteElementDiscretization::getMeshPointer() const {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        FiniteElementDiscretization::reportInvalidParallelMesh();
      }
      return this->parallel_mesh;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        FiniteElementDiscretization::reportInvalidSequentialMesh();
      }
      return this->sequential_mesh;
    }
  }  // end of getMeshPointer

  template <bool parallel>
  inline FiniteElementSpace<parallel>&
  FiniteElementDiscretization::getFiniteElementSpace() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_fe_space.get()) {
        FiniteElementDiscretization::reportInvalidParallelFiniteElementSpace();
      }
      return *(this->parallel_fe_space);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_fe_space.get()) {
        FiniteElementDiscretization::
            reportInvalidSequentialFiniteElementSpace();
      }
      return *(this->sequential_fe_space);
    }
  }  // end of getFiniteElementSpace

  template <bool parallel>
  const FiniteElementSpace<parallel>&
  FiniteElementDiscretization::getFiniteElementSpace() const {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_fe_space.get()) {
        FiniteElementDiscretization::reportInvalidParallelFiniteElementSpace();
      }
      return *(this->parallel_fe_space);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_fe_space.get()) {
        FiniteElementDiscretization::
            reportInvalidSequentialFiniteElementSpace();
      }
      return *(this->sequential_fe_space);
    }

  }  // end of getFiniteElementSpace

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_FINITEELEMENTDISCRETIZATION_IXX */
