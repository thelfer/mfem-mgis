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
  inline Mesh<parallel>& MeshDiscretization::getMesh() noexcept {
    return *(this->template getMeshPointer<parallel>());
  }  // end of getMesh

  template <bool parallel>
  inline const Mesh<parallel>& MeshDiscretization::getMesh() const noexcept {
    return *(this->template getMeshPointer<parallel>());
  }  // end of getMesh

  template <bool parallel>
  std::shared_ptr<Mesh<parallel>> MeshDiscretization::getMutableMeshPointer()
      const noexcept {
    if constexpr (parallel) {
      return this->getMutableParallelMeshPointer();
    } else {
      return this->getMutableSequentialMeshPointer();
    }
  }  // end of getMeshPointer

  template <bool parallel>
  std::shared_ptr<Mesh<parallel>>
  MeshDiscretization::getMeshPointer() noexcept {
    return this->template getMutableMeshPointer<parallel>();
  }

  template <bool parallel>
  std::shared_ptr<const Mesh<parallel>> MeshDiscretization::getMeshPointer()
      const noexcept {
    if constexpr (parallel) {
      return this->getParallelMeshPointer();
    } else {
      return this->getSequentialMeshPointer();
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
  requires((N == 2) || (N == 3))                      //
      OptionalReference<const std::vector<Point<N>>>  //
      MeshDiscretization::getPointsSet(Context& ctx, std::string_view n)
  const noexcept {
    if constexpr (N == 2) {
      return this->getPointsSet2D(ctx, n);
    } else {
      return this->getPointsSet3D(ctx, n);
    }
  }  // end of getPointsSet

  template <size_type N>
  requires((N == 2) || (N == 3))                                             //
      OptionalReference<const std::map<std::string, Point<N>, std::less<>>>  //
      MeshDiscretization::getPoints(Context& ctx)
  const noexcept {
    if constexpr (N == 2) {
      return this->getPoints2D(ctx);
    } else {
      return this->getPoints3D(ctx);
    }
  }  // end of getPoints

  template <size_type N>
  requires((N == 2) || (N == 3))  //
      OptionalReference<
          const std::map<std::string, std::vector<Point<N>>, std::less<>>>  //
      MeshDiscretization::getPointsSets(Context& ctx)
  const noexcept {
    if constexpr (N == 2) {
      return this->getPointsSets2D(ctx);
    } else {
      return this->getPointsSets3D(ctx);
    }
  }  // end of getPointsSets

#endif /* MGIS_HAVE_TFEL */

  template <bool parallel>
  std::shared_ptr<SubMesh<parallel>> MeshDiscretization::getSubMesh(
      Context& ctx, const Parameter& p) const noexcept {
    if constexpr (parallel) {
      return this->getParallelSubMesh(ctx, p);
    } else {
      return this->getSequentialSubMesh(ctx, p);
    }
  }  // end of getSubMesh

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MESHDISCRETIZATION_IXX */
