/*!
 * \file   include/MFEMMGIS/MeshDiscretization.ixx
 * \brief
 * \author Thomas Helfer
 * \date   06/03/2026
 */

#ifndef LIB_MFEM_MGIS_MESHDISCRETIZATION_IXX
#define LIB_MFEM_MGIS_MESHDISCRETIZATION_IXX

namespace mfem_mgis {

  template <bool parallel>
  inline Mesh<parallel>& MeshDiscretization::getMesh() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        MeshDiscretization::reportInvalidParallelMesh();
      }
      return *(this->parallel_mesh);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        MeshDiscretization::reportInvalidSequentialMesh();
      }
      return *(this->sequential_mesh);
    }
  }  // end of getMesh

  template <bool parallel>
  inline const Mesh<parallel>& MeshDiscretization::getMesh() const {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        MeshDiscretization::reportInvalidParallelMesh();
      }
      return *(this->parallel_mesh);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        MeshDiscretization::reportInvalidSequentialMesh();
      }
      return *(this->sequential_mesh);
    }
  }  // end of getMesh

  template <bool parallel>
  std::shared_ptr<Mesh<parallel>> MeshDiscretization::getMeshPointer() {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        MeshDiscretization::reportInvalidParallelMesh();
      }
      return this->parallel_mesh;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        MeshDiscretization::reportInvalidSequentialMesh();
      }
      return this->sequential_mesh;
    }
  }  // end of getMeshPointer

  template <bool parallel>
  std::shared_ptr<const Mesh<parallel>> MeshDiscretization::getMeshPointer()
      const {
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      if (!this->parallel_mesh.get()) {
        MeshDiscretization::reportInvalidParallelMesh();
      }
      return this->parallel_mesh;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (!this->sequential_mesh.get()) {
        MeshDiscretization::reportInvalidSequentialMesh();
      }
      return this->sequential_mesh;
    }
  }  // end of getMeshPointer

#ifdef MGIS_HAVE_TFEL
  template <size_type N>
  requires((N == 2) || (N == 3))  //
      std::optional<Point<N>> MeshDiscretization::getPoint(Context& ctx,
                                                           std::string_view n)
  const noexcept {
    if constexpr (N == 2) {
      return this->getPoint2D(ctx, n);
    } else {
      return this->getPoint3D(ctx, n);
    }
  }  // end of getPoint

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      std::optional<std::vector<Point<N>>> MeshDiscretization::getPointsSet(
          Context& ctx, std::string_view n)
  const noexcept {
    if constexpr (N == 2) {
      return this->getPointsSet2D(ctx,n);
    } else {
      return this->getPointsSet3D(ctx,n);
    }
  }  // end of getPointsSet

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      std::optional<std::map<std::string,Point<N>>> MeshDiscretization::getPoints(Context& ctx)
  const noexcept {
    if constexpr (N == 2) {
      return this->getPoints2D(ctx);
    } else {
      return this->getPoints3D(ctx);
    }
  }  // end of getPoints

#endif /* MGIS_HAVE_TFEL */


}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MESHDISCRETIZATION_IXX */
