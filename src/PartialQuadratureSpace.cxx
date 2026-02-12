/*!
 * \file   src/PartialQuadratureSpace.cxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#include <cmath>
#include <iterator>
#include <algorithm>
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MGIS/Raise.hxx"
#ifdef MFEM_USE_MPI
#include "MFEMMGIS/MPI.hxx"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"

namespace mfem_mgis {

  template <bool parallel>
  static size_type buildPartialQuadratureSpaceOffsets(
      std::unordered_map<size_type, size_type>& offsets,
      std::unordered_map<size_type, size_type>& number_of_quadrature_points,
      const FiniteElementSpace<parallel>& fespace,
      const size_type m,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>&
          integration_rule_selector) {
    auto ng = size_type{};
    for (size_type i = 0; i != fespace.GetNE(); ++i) {
      if (fespace.GetAttribute(i) != m) {
        continue;
      }
      const auto& fe = *(fespace.GetFE(i));
      const auto& tr = *(fespace.GetElementTransformation(i));
      offsets[i] = ng;
      const auto& ir = integration_rule_selector(fe, tr);
      const auto lng = ir.GetNPoints();
      ng += lng;
      number_of_quadrature_points[i] = lng;
    }
    return ng;
  }  // end of buildPartialQuadratureSpaceOffsets

  void PartialQuadratureSpace::treatInvalidElementIndex(const size_type id,
                                                        const size_type i) {
    mgis::raise(
        "PartialQuadratureSpace::getOffset: "
        "invalid element index '" +
        std::to_string(i) + "' for material '" + std::to_string(id) + "'");
  }  // end of treatInvalidElementIndex

  PartialQuadratureSpace::PartialQuadratureSpace(
      const FiniteElementDiscretization& fed,
      const size_type m,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>& irs)
      : fe_discretization(fed), integration_rule_selector(irs), id(m) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& fespace =
          this->fe_discretization.getFiniteElementSpace<true>();
      this->ng = buildPartialQuadratureSpaceOffsets<true>(
          this->offsets, this->number_of_quadrature_points, fespace, this->id,
          this->integration_rule_selector);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      const auto& fespace =
          this->fe_discretization.getFiniteElementSpace<false>();
      this->ng = buildPartialQuadratureSpaceOffsets<false>(
          this->offsets, this->number_of_quadrature_points, fespace, this->id,
          this->integration_rule_selector);
    }
  }  // end of PartialQuadratureSpace

  const mfem::IntegrationRule& PartialQuadratureSpace::getIntegrationRule(
      const mfem::FiniteElement& e,
      const mfem::ElementTransformation& tr) const {
    return this->integration_rule_selector(e, tr);
  }

  const FiniteElementDiscretization&
  PartialQuadratureSpace::getFiniteElementDiscretization() const {
    return this->fe_discretization;
  }  // end of getFiniteElementDiscretization

  std::optional<size_type> PartialQuadratureSpace::getNumberOfQuadraturePoints(
      Context& ctx, const size_type e) const noexcept {
    const auto p = this->number_of_quadrature_points.find(e);
    if (p == this->number_of_quadrature_points.end()) {
      return ctx.registerErrorMessage("invalid element index '" +
                                      std::to_string(e) + "'");
    }
    return p->second;
  }  // end of getNumberOfQuadraturePoints

  PartialQuadratureSpace::~PartialQuadratureSpace() = default;

  template <bool parallel>
  [[nodiscard]] static std::optional<std::map<mfem::Geometry::Type, size_type>>
  getNumberOfQuadraturePointsByGeometricElementType(
      Context& ctx, const PartialQuadratureSpace& s) noexcept {
    const auto& fed = s.getFiniteElementDiscretization();
    const auto& mesh = fed.getMesh<parallel>();
    auto qmapping = std::map<mfem::Geometry::Type, size_type>{};
    for (const auto [e, o] : s.getOffsets()) {
      const auto gtype = mesh.GetElementGeometry(e);
      if (qmapping.contains(gtype)) {
        continue;
      }
      const auto on = s.getNumberOfQuadraturePoints(ctx, e);
      if (isInvalid(on)) {
        return {};
      }
      qmapping.insert({gtype, *on});
    }
    return qmapping;
  }  // end of getNumberOfQuadraturePointsByGeometricElementType

  std::optional<PartialQuadratureSpaceInformation> getInformation(
      Context& ctx, const PartialQuadratureSpace& s) noexcept {
    const auto& fed = s.getFiniteElementDiscretization();
    auto info = PartialQuadratureSpaceInformation{};
    info.identifier = s.getId();
    const auto oname = fed.getMaterialName(ctx, s.getId());
    if (isInvalid(oname)) {
      return {};
    }
    info.name = *oname;
    info.number_of_finite_elements = getNumberOfCells(s);
    info.number_of_quadrature_points = getNumberOfElements(s);
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto oqmapping =
          getNumberOfQuadraturePointsByGeometricElementType<true>(ctx, s);
      if (isInvalid(oqmapping)) {
        return {};
      }
      info.number_of_quadrature_points_by_geometric_type = *oqmapping;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      const auto oqmapping =
          getNumberOfQuadraturePointsByGeometricElementType<false>(ctx, s);
      if (isInvalid(oqmapping)) {
        return {};
      }
      info.number_of_quadrature_points_by_geometric_type = *oqmapping;
    }
    return info;
  }  // end of getInformation

  [[nodiscard]] static std::optional<std::string_view> to_string(
      Context& ctx, const mfem::Geometry::Type t) noexcept {
    if (t == mfem::Geometry::INVALID) {
      return ctx.registerErrorMessage("invalid geometric type");
    } else if (t == mfem::Geometry::POINT) {
      return "POINT";
    } else if (t == mfem::Geometry::SEGMENT) {
      return "SEGMENT";
    } else if (t == mfem::Geometry::TRIANGLE) {
      return "TRIANGLE";
    } else if (t == mfem::Geometry::SQUARE) {
      return "QUADRANGLE (SQUARE)";
    } else if (t == mfem::Geometry::TETRAHEDRON) {
      return "TETRAHEDRON";
    } else if (t == mfem::Geometry::CUBE) {
      return "HEXAHEDRON (CUBE)";
    } else if (t == mfem::Geometry::PRISM) {
      return "PRISM";
    } else if (t == mfem::Geometry::PYRAMID) {
      return "PYRAMID";
    }
    return ctx.registerErrorMessage("unsupported geometric type");
  }  // end of to_string

  bool info(Context& ctx,
            std::ostream& os,
            const PartialQuadratureSpace& s) noexcept {
    auto oinfo = getInformation(ctx, s);
    if (isInvalid(oinfo)) {
      return false;
    }
    return info(ctx, os, *oinfo);
  }  // end of info

  bool info(Context& ctx,
            std::ostream& os,
            const PartialQuadratureSpaceInformation& info) noexcept {
    auto success = true;
    os << "- material:";
    if (!info.name.empty()) {
      os << " '" << info.name << "' ";
    }
    os << " " << std::to_string(info.identifier) << '\n';
    os << "- number of elements: " << info.number_of_finite_elements << '\n';
    os << "- number of quadrature points: " << info.number_of_quadrature_points
       << '\n';
    if (!info.number_of_quadrature_points_by_geometric_type.empty()) {
      os << "- number of quadrature points per element type:\n";
      for (const auto& [g, n] :
           info.number_of_quadrature_points_by_geometric_type) {
        const auto ogn = to_string(ctx, g);
        const auto ok = isValid(ogn);
        if (ok) {
          os << "  - " << *ogn << ": " << n << '\n';
        }
        success = success && ok;
      }
    }
    return success;
  }  // end of info

  std::optional<PartialQuadratureSpaceInformation> synchronize(
      Context& ctx, const PartialQuadratureSpaceInformation& info) noexcept {
#ifdef MFEM_USE_MPI
    auto r = PartialQuadratureSpaceInformation{};
    int nprocesses;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
    std::vector<size_type> ids(nprocesses);
    MPI_Allgather(&(info.identifier), 1, mpi_type<size_type>, ids.data(), 1,
                  mpi_type<size_type>, MPI_COMM_WORLD);
    if (std::adjacent_find(ids.begin(), ids.end(), std::not_equal_to<>()) !=
        ids.end()) {
      return ctx.registerErrorMessage(
          "the given information do not refer to the same material on all "
          "processes");
    }
    //           n, as_mpi_datatype<T>::value(), MPI_COMM_WORLD);
    //     int MPI_Allreduce(const void* sendbuf, void* recvbuf, int count,
    //                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
    // get all ids
    //     using T = typename C::value_type;
    //     Index n = getSize(buf);
    //     Vector<T> output(MPI::nProc() * n);
    //
    //     Int ierr;
    //     if constexpr (std::is_fundamental_v<T>) {
    //       ierr = MPI_Allgather(
    //           buf.data(), n, as_mpi_datatype<T>::value(), output.data(),
    //           n, as_mpi_datatype<T>::value(), MPI_COMM_WORLD);
    //     }
    //     else {
    //       ierr = MPI_Allgather(buf.data(), n * sizeof(T), MPI_BYTE,
    //       output.data(), n * sizeof(T), MPI_BYTE, MPI_COMM_WORLD);
    //     }
    //     MPI::check(ierr);
    //
    //     return output;
    return r;
#else  /* MFEM_USE_MPI */
    return info;
#endif /* MFEM_USE_MPI */
  }    // end of synchronize

}  // end of namespace mfem_mgis
