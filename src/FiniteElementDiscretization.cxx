/*!
 * \file   src/FiniteElementDiscretization.cxx
 * \brief
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#include <regex>
#include <cctype>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <mfem/mesh/mesh.hpp>
#include <mfem/fem/fespace.hpp>
#ifdef MFEM_USE_MPI
#include <mfem/mesh/pmesh.hpp>
#include <mfem/fem/pfespace.hpp>
#endif
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

  const char* const FiniteElementDiscretization::FiniteElementFamily =
      "FiniteElementFamily";
  const char* const FiniteElementDiscretization::FiniteElementOrder =
      "FiniteElementOrder";
  const char* const FiniteElementDiscretization::UnknownsSize = "UnknownsSize";

  void FiniteElementDiscretization::reportInvalidParallelFiniteElementSpace() {
    raise(
        "FiniteElementDiscretization::reportInvalidParallelFiniteElementSpace: "
        "no parallel finite element space defined");
  }  // end of reportInvalidParallelFiniteElementSpace

  void
  FiniteElementDiscretization::reportInvalidSequentialFiniteElementSpace() {
    raise(
        "FiniteElementDiscretization::"
        "reportInvalidSequentialFiniteElementSpace: "
        "no sequential finite element space defined");
  }  // end of reportInvalidSequentialFiniteElementSpace

  [[nodiscard]] static std::vector<std::string>
  getFiniteElementDiscretizationParametersList() {
    return {FiniteElementDiscretization::FiniteElementFamily,
            FiniteElementDiscretization::FiniteElementOrder,
            FiniteElementDiscretization::UnknownsSize};
  }  // end of getFiniteElementDiscretizationParametersList

  std::vector<std::string> FiniteElementDiscretization::getParametersList() {
    auto d = MeshDiscretization::getParametersList();
    const auto names = getFiniteElementDiscretizationParametersList();
    d.insert(d.end(), names.begin(), names.end());
    return d;
  }  // end of getParametersList

  template <bool parallel>
  std::pair<std::shared_ptr<const FiniteElementCollection>,
            std::unique_ptr<FiniteElementSpace<parallel>>>
  buildFiniteElementCollectionAndSpace(MeshDiscretization& m,
                                       const Parameters& params) {
    checkParameters(throwing, params,
                    getFiniteElementDiscretizationParametersList());
    const auto& fe_family = get_if<std::string>(
        throwing, params, FiniteElementDiscretization::FiniteElementFamily,
        "H1");
    const auto fe_order = get_if<int>(
        throwing, params, FiniteElementDiscretization::FiniteElementOrder, 1);
    const auto u_size =
        get<int>(throwing, params, FiniteElementDiscretization::UnknownsSize);
    // building the finite element collection
    if (fe_family != "H1") {
      raise(
          "FiniteElementDiscretization::FiniteElementDiscretization: "
          "unsupported finite element family '" +
          fe_family + "'");
    }
    auto fec =
        std::make_shared<mfem::H1_FECollection>(fe_order, getSpaceDimension(m));
    // building the finite element space
    if constexpr (parallel) {
#ifdef MFEM_USE_MPI
      return {fec, std::make_unique<FiniteElementSpace<true>>(
                       m.getMeshPointer<true>().get(), fec.get(), u_size)};
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      return {fec, std::make_unique<FiniteElementSpace<false>>(
                       m.getMeshPointer<false>().get(), fec.get(), u_size)};
    }
  }  // end of buildFiniteElementCollectionAndSpace

  FiniteElementDiscretization::FiniteElementDiscretization(
      const MeshDiscretization& m, const Parameters& params)
      : MeshDiscretization(m) {
    CatchTimeSection("FED::Constructor");
    checkParameters(throwing, params,
                    FiniteElementDiscretization::getParametersList());
    if (this->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      std::tie(this->fec, this->parallel_fe_space) =
          buildFiniteElementCollectionAndSpace<true>(*this, params);
#else
      reportUnsupportedParallelComputations();
#endif
    } else {
      std::tie(this->fec, this->sequential_fe_space) =
          buildFiniteElementCollectionAndSpace<false>(*this, params);
    }
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<true>> m, const Parameters& params)
      : MeshDiscretization(m) {
    CatchTimeSection("FED::Constructor");
    checkParameters(throwing, params,
                    FiniteElementDiscretization::getParametersList());
#ifdef MFEM_USE_MPI
    std::tie(this->fec, this->parallel_fe_space) =
        buildFiniteElementCollectionAndSpace<true>(*this, params);
#else
    reportUnsupportedParallelComputations();
#endif
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<false>> m, const Parameters& params)
      : MeshDiscretization(m) {
    std::tie(this->fec, this->sequential_fe_space) =
        buildFiniteElementCollectionAndSpace<false>(*this, params);
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      const Parameters& params)
      : MeshDiscretization(extract(
            throwing, params, MeshDiscretization::getParametersList())) {
    CatchTimeSection("FED::Constructor");
    if (this->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      std::tie(this->fec, this->parallel_fe_space) =
          buildFiniteElementCollectionAndSpace<true>(
              *this, remove(params, MeshDiscretization::getParametersList()));
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      std::tie(this->fec, this->sequential_fe_space) =
          buildFiniteElementCollectionAndSpace<false>(
              *this, remove(params, MeshDiscretization::getParametersList()));
    }
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<true>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      const size_type d)
      : MeshDiscretization(std::move(m)), fec(std::move(c)) {
    if (fec.get() == nullptr) {
      raise("invalid finite element collection");
    }
#ifdef MFEM_USE_MPI
    this->parallel_fe_space = std::make_unique<FiniteElementSpace<true>>(
        this->parallel_mesh.get(), this->fec.get(), d);
#else  /* MFEM_USE_MPI */
    static_cast<void>(d);
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<true>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      std::unique_ptr<FiniteElementSpace<true>> s)
      : MeshDiscretization(std::move(m)),
        fec(std::move(c))
#ifdef MFEM_USE_MPI
        ,
        parallel_fe_space(std::move(s))
#endif /* MFEM_USE_MPI */
  {
#ifdef MFEM_USE_MPI
    if (this->fec.get() == nullptr) {
      raise("invalid finite element collection");
    }
    if (this->parallel_fe_space.get() == nullptr) {
      raise("invalid finite element space");
    }
    if (this->parallel_mesh.get() != this->parallel_fe_space->GetMesh()) {
      raise(
          "FiniteElementDiscretization::FiniteElementDiscretization: "
          "mesh pointer don't match the mesh on which the finite element space "
          "is built");
    }
#else  /* MFEM_USE_MPI */
    static_cast<void>(s);
    reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }    // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<false>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      const size_type d)
      : MeshDiscretization(m), fec(std::move(c)) {
    if (this->fec.get() == nullptr) {
      raise("invalid finite element collection");
    }
    this->sequential_fe_space = std::make_unique<FiniteElementSpace<false>>(
        this->sequential_mesh.get(), this->fec.get(), d);
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<false>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      std::unique_ptr<FiniteElementSpace<false>> s)
      : MeshDiscretization(std::move(m)),
        fec(std::move(c)),
        sequential_fe_space(std::move(s)) {
    if (this->fec.get() == nullptr) {
      raise("invalid finite element collection");
    }
    if (this->sequential_fe_space.get() == nullptr) {
      raise("invalid finite element space");
    }
    if (this->sequential_mesh.get() != this->sequential_fe_space->GetMesh()) {
      raise(
          "FiniteElementDiscretization::FiniteElementDiscretization: "
          "mesh pointer don't match the mesh on which the finite element space "
          "is built");
    }
  }  // end of FiniteElementDiscretization

  const FiniteElementCollection&
  FiniteElementDiscretization::getFiniteElementCollection() const noexcept {
    return *(this->fec);
  }  // end of getFiniteElementCollection

  std::shared_ptr<const FiniteElementCollection>
  FiniteElementDiscretization::getFiniteElementCollectionPointer()
      const noexcept {
    return this->fec;
  }  // end of getFiniteElementCollection

  FiniteElementDiscretization::~FiniteElementDiscretization() = default;

  size_type getTrueVSize(const FiniteElementDiscretization& fed) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return fed.getFiniteElementSpace<true>().GetTrueVSize();
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return fed.getFiniteElementSpace<false>().GetTrueVSize();
  }  // end of getTrueVSize

  template <>
  bool getInformation<FiniteElementDiscretization>(
      Context& ctx,
      std::ostream& os,
      const FiniteElementDiscretization& fed) noexcept {
    if (!getInformation(ctx, os, static_cast<const MeshDiscretization&>(fed))) {
      return false;
    }
    os << "\n\n# Finite element space\n\n"
       << "- true vector size: " << getTrueVSize(fed) << '\n';
    return true;
  }  // end of info

}  // end of namespace mfem_mgis
