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
#include "MFEMMGIS/MPI.hxx"
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
      : fe_discretization(fed),
#ifdef MFEM_USE_MPI
        parallel_fespace(fed.describesAParallelComputation()
                             ? &(fed.getFiniteElementSpace<true>())
                             : nullptr),
        sequential_fespace(!fed.describesAParallelComputation()
                               ? &(fed.getFiniteElementSpace<false>())
                               : nullptr),
#else  /* MFEM_USE_MPI */
        sequential_fespace(fed.getFiniteElementSpace<false>()),
#endif /* MFEM_USE_MPI */
        integration_rule_selector(irs),
        id(m) {
    this->initialize(throwing);
  }  // end of PartialQuadratureSpace

#ifdef MFEM_USE_MPI
  [[nodiscard]] static const FiniteElementSpace<true>*
  initializeParallelFiniteElementSpace(attributes::Throwing,
                                       const FiniteElementDiscretization& fed,
                                       const LocationIdentifier l) {
    if (isInvalid(l)) {
      raise("invalid location identifier");
    }
    if (!fed.describesAParallelComputation()) {
      return nullptr;
    }
    if (isValid(l.material_identifier)) {
      return &(fed.getFiniteElementSpace<true>());
    }
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    auto fes_manager = fed.getFiniteElementSpacesManager();
    auto ofes =
        fes_manager.getFiniteElementSpace<true>(
            ctx,
            FiniteElementSpacesManager::GetFiniteElementSpaceOnSubMeshArguments{
                .location = MeshDiscretization::Location::ON_BOUNDARIES,
                .identifiers = Parameter{l.boundary_identifier->id},
                .number_of_components =
                    fed.getFiniteElementSpace<true>().GetVDim()}) |
        or_raise;
    return ofes.get();
  }    // end of initializeParallelFiniteElementSpace
#endif /* MFEM_USE_MPI */

  [[nodiscard]] static const FiniteElementSpace<false>*
  initializeSequentialFiniteElementSpace(attributes::Throwing,
                                         const FiniteElementDiscretization& fed,
                                         const LocationIdentifier l) {
    if (isInvalid(l)) {
      raise("invalid location identifier");
    }
    if (fed.describesAParallelComputation()) {
      return nullptr;
    }
    if (isValid(l.material_identifier)) {
      return &(fed.getFiniteElementSpace<false>());
    }
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    auto fes_manager = fed.getFiniteElementSpacesManager();
    auto ofes =
        fes_manager.getFiniteElementSpace<false>(
            ctx,
            FiniteElementSpacesManager::GetFiniteElementSpaceOnSubMeshArguments{
                .location = MeshDiscretization::Location::ON_BOUNDARIES,
                .identifiers = l.boundary_identifier->id,
                .number_of_components =
                    fed.getFiniteElementSpace<false>().GetVDim()}) |
        or_raise;
    return ofes.get();
  }  // end of initializeSequentialFiniteElementSpace

  [[nodiscard]] static size_type getIdentifier(attributes::Throwing,
                                               const LocationIdentifier& l) {
    if (isInvalid(l)) {
      raise("invalid location identifier");
    }
    if (isValid(l.material_identifier)) {
      return l.material_identifier.value().id;
    }
    return l.boundary_identifier.value().id;
  }  // end of getIdentifier

  PartialQuadratureSpace::PartialQuadratureSpace(
      const FiniteElementDiscretization& fed,
      const LocationIdentifier& l,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>& irs)
      : fe_discretization(fed),
#ifdef MFEM_USE_MPI
        parallel_fespace(
            initializeParallelFiniteElementSpace(throwing, fed, l)),
        sequential_fespace(
            initializeSequentialFiniteElementSpace(throwing, fed, l)),
#else  /* MFEM_USE_MPI */
        sequential_fespace(
            initializeSequentialFiniteElementSpace(throwing, fed, l)),
#endif /* MFEM_USE_MPI */
        integration_rule_selector(irs),
        id(getIdentifier(throwing, l)) {
    this->initialize(throwing);
  }  // end of PartialQuadratureSpace

#ifdef MFEM_USE_MPI
  PartialQuadratureSpace::PartialQuadratureSpace(
      const FiniteElementDiscretization& fed,
      const FiniteElementSpace<true>& fespace,
      const size_type l,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>& irs)
      : fe_discretization(fed),
        parallel_fespace(&fespace),
        integration_rule_selector(irs),
        id(l) {
    this->initialize(throwing);
  }    // end of PartialQuadratureSpace
#endif /* MFEM_USE_MPI */

  PartialQuadratureSpace::PartialQuadratureSpace(
      const FiniteElementDiscretization& fed,
      const FiniteElementSpace<false>& fespace,
      const size_type l,
      const std::function<const mfem::IntegrationRule&(
          const mfem::FiniteElement&, const mfem::ElementTransformation&)>& irs)
      : fe_discretization(fed),
        sequential_fespace(&fespace),
        integration_rule_selector(irs),
        id(l) {
    this->initialize(throwing);
  }  // end of PartialQuadratureSpace

  void PartialQuadratureSpace::initialize(attributes::Throwing) {
    if (!(this->integration_rule_selector)) {
      raise("invalid quadrature rule selector");
    }
    auto fespaces_manager =
        this->fe_discretization.getFiniteElementSpacesManager();
#ifdef MFEM_USE_MPI
    if (this->parallel_fespace != nullptr) {
      if (!fespaces_manager.manages(*(this->parallel_fespace))) {
        raise(
            "the given finite element space is not managed "
            "by the given finite element discretization");
      }
      const auto& attributes = this->parallel_fespace->GetParMesh()->attributes;
      if (attributes.Find(this->id) == -1) {
        raise("invalid location identifier (" + std::to_string(this->id) + ")");
      }
      this->ng = buildPartialQuadratureSpaceOffsets<true>(
          this->offsets, this->number_of_quadrature_points,
          *(this->parallel_fespace), this->id, this->integration_rule_selector);
      return;
    }
#endif /* MFEM_USE_MPI */
    if (!fespaces_manager.manages(*(this->sequential_fespace))) {
      raise(
          "the given finite element space is not managed by "
          "the given finite element discretization");
    }
    const auto& attributes = this->sequential_fespace->GetMesh()->attributes;
    if (attributes.Find(this->id) == -1) {
      raise("invalid location identifier (" + std::to_string(this->id) + ")");
    }
    this->ng = buildPartialQuadratureSpaceOffsets<false>(
        this->offsets, this->number_of_quadrature_points,
        *(this->sequential_fespace), this->id, this->integration_rule_selector);
  }  // end of initialize

  const mfem::IntegrationRule& PartialQuadratureSpace::getIntegrationRule(
      const mfem::FiniteElement& e,
      const mfem::ElementTransformation& tr) const {
    return this->integration_rule_selector(e, tr);
  }

  std::string PartialQuadratureSpace::getLocationName() const noexcept {
    auto ctx = Context{};
    const auto ol = [&ctx, this] {
#ifdef MFEM_USE_MPI
      if (this->parallel_fespace != nullptr) {
        return this->fe_discretization.getLocationIdentifier(
            ctx, *(this->parallel_fespace->GetParMesh()), this->getId());
      }
#endif /* MFEM_USE_MPI */
      ctx.assertOrTerminate(this->sequential_fespace != nullptr,
                            "internal error");
      return this->fe_discretization.getLocationIdentifier(
          ctx, *(this->sequential_fespace->GetMesh()), this->getId());
    }();
    ctx.assertOrTerminate(isValid(ol), "internal error");
    if (isValid(ol->material_identifier)) {
      const auto oname = this->fe_discretization.getMaterialName(
          ctx, ol->material_identifier->id);
      if (isValid(oname)) {
        if (!oname->empty()) {
          return *oname;
        }
      }
      return "material (" + std::to_string(this->getId()) + ")";
    }
    ctx.assertOrTerminate(isValid(ol->boundary_identifier), "internal error");
    const auto oname = this->fe_discretization.getBoundaryName(
        ctx, ol->boundary_identifier->id);
    if (isValid(oname)) {
      if (!oname->empty()) {
        return *oname;
      }
    }
    return "boundary (" + std::to_string(this->getId()) + ")";
  }  // end of getLocationName

  const MeshDiscretization& PartialQuadratureSpace::getMeshDiscretization()
      const noexcept {
    return static_cast<const MeshDiscretization&>(this->fe_discretization);
  }  // end of getMeshDiscretization

  const FiniteElementDiscretization&
  PartialQuadratureSpace::getFiniteElementDiscretization() const noexcept {
    return this->fe_discretization;
  }  // end of getFiniteElementDiscretization

  OptionalReference<const Mesh<true>> PartialQuadratureSpace::getParallelMesh(
      Context& ctx) const noexcept {
#ifdef MFEM_USE_MPI
    auto ofes = this->getParallelFiniteElementSpace(ctx);
    if (isInvalid(ofes)) {
      return {};
    }
    return {ofes->GetParMesh()};
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of getParallelMesh

  OptionalReference<const Mesh<false>>
  PartialQuadratureSpace::getSequentialMesh(Context& ctx) const noexcept {
    auto ofes = this->getSequentialFiniteElementSpace(ctx);
    if (isInvalid(ofes)) {
      return {};
    }
    return {ofes->GetMesh()};
  }  // end of getSequentialMesh

  OptionalReference<const FiniteElementSpace<true>>
  PartialQuadratureSpace::getParallelFiniteElementSpace(
      Context& ctx) const noexcept {
#ifdef MFEM_USE_MPI
    if (this->parallel_fespace == nullptr) {
      return ctx.registerErrorMessage(
          "the partial quadrature space is not built on a parallel mesh");
    }
    return {this->parallel_fespace};
#else  /* MFEM_USE_MPI */
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of getParallelFiniteElementSpace

  OptionalReference<const FiniteElementSpace<false>>
  PartialQuadratureSpace::getSequentialFiniteElementSpace(
      Context& ctx) const noexcept {
    if (this->sequential_fespace == nullptr) {
      return ctx.registerErrorMessage(
          "the partial quadrature space is not built on a sequential mesh");
    }
    return {this->sequential_fespace};
  }  // end of getSequentialFiniteElementSpace

  bool PartialQuadratureSpace::isDefinedOnAMaterial() const {
    auto ctx = Context{};
#ifdef MFEM_USE_MPI
    if (this->parallel_fespace != nullptr) {
      const auto ook = this->fe_discretization.isDefinedOnMaterials(
          ctx, *(this->parallel_fespace->GetParMesh()));
      ctx.assertOrTerminate(isValid(ook), "internal error");
      return *ook;
    }
#endif /* MFEM_USE_MPI */
    ctx.assertOrTerminate(this->sequential_fespace != nullptr,
                          "internal error");
    const auto ook = this->fe_discretization.isDefinedOnMaterials(
        ctx, *(this->sequential_fespace->GetMesh()));
    ctx.assertOrTerminate(isValid(ook), "internal error");
    return *ook;
  }  // end of isDefinedOnAMaterial

  bool PartialQuadratureSpace::isDefinedOnABoundary() const {
    auto ctx = Context{};
#ifdef MFEM_USE_MPI
    if (this->parallel_fespace != nullptr) {
      const auto ook = this->fe_discretization.isDefinedOnBoundaries(
          ctx, *(this->parallel_fespace->GetParMesh()));
      ctx.assertOrTerminate(isValid(ook), "internal error");
      return *ook;
    }
#endif /* MFEM_USE_MPI */
    ctx.assertOrTerminate(this->sequential_fespace != nullptr,
                          "internal error");
    const auto ook = this->fe_discretization.isDefinedOnBoundaries(
        ctx, *(this->sequential_fespace->GetMesh()));
    ctx.assertOrTerminate(isValid(ook), "internal error");
    return *ook;
  }  // end of isDefinedOnABoundary

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
  [[nodiscard]] static std::map<mfem::Geometry::Type, size_type>
  getNumberOfElementsByGeometricElementType(
      const PartialQuadratureSpace& s) noexcept {
    auto ctx = Context{};
    auto or_die = ctx.getFatalFailureHandler();
    const auto& mesh = s.getMesh<parallel>(ctx) | or_die;
    auto emapping = std::map<mfem::Geometry::Type, size_type>{};
    for (const auto [e, o] : s.getOffsets()) {
      const auto gtype = mesh.GetElementGeometry(e);
      ++(emapping[gtype]);
    }
    return emapping;
  }  // end of getNumberOfElementsByGeometricElementType

  template <bool parallel>
  [[nodiscard]] static std::optional<std::map<mfem::Geometry::Type, size_type>>
  getNumberOfQuadraturePointsByGeometricElementType(
      Context& ctx, const PartialQuadratureSpace& s) noexcept {
    const auto& omesh = s.getMesh<parallel>(ctx);
    if (isInvalid(omesh)) {
      return {};
    }
    auto qmapping = std::map<mfem::Geometry::Type, size_type>{};
    for (const auto [e, o] : s.getOffsets()) {
      const auto gtype = omesh->GetElementGeometry(e);
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

  std::optional<PartialQuadratureSpaceInformation> getLocalInformation(
      Context& ctx, const PartialQuadratureSpace& s) noexcept {
    const auto& fed = s.getFiniteElementDiscretization();
    auto info = PartialQuadratureSpaceInformation{};
    info.identifier = s.getId();
    info.name = s.getLocationName();
    info.number_of_cells = getNumberOfCells(s);
    info.number_of_quadrature_points = getNumberOfElements(s);
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      info.number_of_cells_by_geometric_type =
          getNumberOfElementsByGeometricElementType<true>(s);
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
      info.number_of_cells_by_geometric_type =
          getNumberOfElementsByGeometricElementType<false>(s);
      const auto oqmapping =
          getNumberOfQuadraturePointsByGeometricElementType<false>(ctx, s);
      if (isInvalid(oqmapping)) {
        return {};
      }
      info.number_of_quadrature_points_by_geometric_type = *oqmapping;
    }
#ifdef MFEM_USE_MPI
    info.communicator = getMPICommunicator(fed);
#endif /* MFEM_USE_MPI */
    return info;
  }  // end of getLocalInformation

  std::optional<PartialQuadratureSpaceInformation> getInformation(
      Context& ctx, const PartialQuadratureSpace& s) noexcept {
    const auto linfo = getLocalInformation(ctx, s);
    const auto& fed = s.getFiniteElementDiscretization();
    if (!isValidOnAllProcesses(fed, linfo)) {
      return {};
    }
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return synchronize(ctx, *linfo);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return linfo;
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
      return "QUADRANGLE (mfem::Geometry::SQUARE)";
    } else if (t == mfem::Geometry::TETRAHEDRON) {
      return "TETRAHEDRON";
    } else if (t == mfem::Geometry::CUBE) {
      return "HEXAHEDRON (mfem::Geometry::CUBE)";
    } else if (t == mfem::Geometry::PRISM) {
      return "PRISM";
    } else if (t == mfem::Geometry::PYRAMID) {
      return "PYRAMID";
    }
    return ctx.registerErrorMessage("unsupported geometric type");
  }  // end of to_string

  template <>
  bool getInformation<PartialQuadratureSpace>(
      Context& ctx,
      std::ostream& os,
      const PartialQuadratureSpace& s) noexcept {
    auto oinfo = getInformation(ctx, s);
    if (isInvalid(oinfo)) {
      return false;
    }
    return getInformation(ctx, os, *oinfo);
  }  // end of getInformation

  template <>
  bool getInformation<PartialQuadratureSpaceInformation>(
      Context& ctx,
      std::ostream& os,
      const PartialQuadratureSpaceInformation& info) noexcept {
    auto success = true;
    os << "- material or boundary identifier:";
    if (!info.name.empty()) {
      os << " '" << info.name << "' ";
    }
    os << " " << std::to_string(info.identifier) << '\n';
    os << "- number of elements: " << info.number_of_cells << '\n';
    os << "- number of quadrature points: " << info.number_of_quadrature_points
       << '\n';
    if (!info.number_of_cells_by_geometric_type.empty()) {
      os << "- number of elements per geometric type:\n";
      for (const auto& [g, n] : info.number_of_cells_by_geometric_type) {
        const auto ogn = to_string(ctx, g);
        const auto ok = isValid(ogn);
        if (ok) {
          os << "  - " << *ogn << ": " << n << '\n';
        }
        success = success && ok;
      }
    }
    if (!info.number_of_quadrature_points_by_geometric_type.empty()) {
      os << "- number of quadrature points per geometric type:\n";
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
  }  // end of getInformation

  std::optional<PartialQuadratureSpaceInformation> synchronize(
      Context& ctx, const PartialQuadratureSpaceInformation& info) noexcept {
#ifdef MFEM_USE_MPI
    auto r = PartialQuadratureSpaceInformation{};
    // paranoïc check
    int nprocesses;
    MPI_Comm_size(info.communicator, &nprocesses);
    std::vector<size_type> ids(nprocesses);
    MPI_Allgather(&(info.identifier), 1, mpi_type<size_type>, ids.data(), 1,
                  mpi_type<size_type>, info.communicator);
    if (std::adjacent_find(ids.begin(), ids.end(), std::not_equal_to<>()) !=
        ids.end()) {
      return ctx.registerErrorMessage(
          "the given information do not refer to the same material on all "
          "processes");
    }
    r.identifier = info.identifier;
    // assuming the name is the same on all process as the identifier is
    r.name = info.name;
    //
    r.number_of_cells = info.number_of_cells;
    MPI_Allreduce(MPI_IN_PLACE, &(r.number_of_cells), 1, mpi_type<size_type>,
                  MPI_SUM, info.communicator);
    r.number_of_quadrature_points = info.number_of_quadrature_points;
    MPI_Allreduce(MPI_IN_PLACE, &(r.number_of_quadrature_points), 1,
                  mpi_type<size_type>, MPI_SUM, info.communicator);
    //
    auto nelts = std::array<size_type, mfem::Geometry::NumGeom>{};
    for (const auto [g, n] : info.number_of_cells_by_geometric_type) {
      nelts[static_cast<std::size_t>(g)] = n;
    }
    MPI_Allreduce(MPI_IN_PLACE, nelts.data(), mfem::Geometry::NumGeom,
                  mpi_type<size_type>, MPI_SUM, info.communicator);
    for (size_type i = 0; const auto& n : nelts) {
      if (n != 0) {
        r.number_of_cells_by_geometric_type.insert(
            {static_cast<mfem::Geometry::Type>(i), n});
      }
      ++i;
    }
    //
    auto nqpoints = std::array<size_type, mfem::Geometry::NumGeom>{};
    for (const auto [g, n] :
         info.number_of_quadrature_points_by_geometric_type) {
      nqpoints[static_cast<std::size_t>(g)] = n;
    }
    for (size_type g = 0; auto& n : nqpoints) {
      std::vector<size_type> all_nqpoints(nprocesses);
      MPI_Allgather(&n, 1, mpi_type<size_type>, all_nqpoints.data(), 1,
                    mpi_type<size_type>, info.communicator);
      const auto pmax =
          std::max_element(all_nqpoints.begin(), all_nqpoints.end());
      if (pmax == all_nqpoints.end()) {  // avoid warning
        abort("internal error");
      }
      // paranoïac check
      for (const auto& n2 : all_nqpoints) {
        if ((n2 != 0) && (n2 != *pmax)) {
          const auto ogn = to_string(ctx, static_cast<mfem::Geometry::Type>(g));
          if (!isValid(ogn)) {
            return {};
          }
          return ctx.registerErrorMessage(
              "the number of quadrature points for geometric type '" +
              std::string{*ogn} + "' is not the same on all processes");
        }
      }
      n = *pmax;
      ++g;
    }
    for (size_type i = 0; const auto& n : nqpoints) {
      if (n != 0) {
        r.number_of_quadrature_points_by_geometric_type.insert(
            {static_cast<mfem::Geometry::Type>(i), n});
      }
      ++i;
    }
    return r;
#else  /* MFEM_USE_MPI */
    return info;
#endif /* MFEM_USE_MPI */
  }    // end of synchronize

  bool areEquivalent(const PartialQuadratureSpace& s1,
                     const PartialQuadratureSpace& s2) noexcept {
    if (&s1 == &s2) {
      return true;
    }
    const auto& fed1 = s1.getFiniteElementDiscretization();
    const auto& fed2 = s2.getFiniteElementDiscretization();
    if (fed1.describesAParallelComputation() !=
        fed2.describesAParallelComputation()) {
      return false;
    }
    if (fed1.describesAParallelComputation()) {
      auto fes_manager = fed1.getFiniteElementSpacesManager();
#ifdef MFEM_USE_MPI
      if (!fes_manager.manages(fed2.getFiniteElementSpace<true>())) {
        return false;
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      auto fes_manager = fed1.getFiniteElementSpacesManager();
      if (!fes_manager.manages(fed2.getFiniteElementSpace<false>())) {
        return false;
      }
    }
    if (s1.getId() != s2.getId()) {
      return false;
    }
    const auto success = [&s1, &s2] {
      if (getSpaceSize(s1) != getSpaceSize(s2)) {
        return false;
      }
      const auto& offsets1 = s1.getOffsets();
      const auto& offsets2 = s2.getOffsets();
      if (offsets1.size() != offsets2.size()) {
        return false;
      }
      auto p1 = offsets1.begin();
      auto p2 = offsets2.begin();
      for (; p1 != offsets1.end(); ++p1, ++p2) {
        if ((p1->first != p2->first) || (p1->second != p2->second)) {
          return false;
        }
      }
      return true;
    }();
    return isTrueOnAllProcesses(fed1, success);
  }  // end of areEquivalent

}  // end of namespace mfem_mgis
