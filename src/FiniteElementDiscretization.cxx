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
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

  const char* const FiniteElementDiscretization::UnknownsSize = "UnknownsSize";

  //! list of valid parametres when the mesh is alredy built
  [[nodiscard]] static std::vector<std::string>
  getFiniteElementDiscretizationParametersList() {
    auto d =
        FiniteElementSpacesManager::getFiniteElementCollectionParametersList();
    d.push_back(FiniteElementDiscretization::UnknownsSize);
    return d;
  }  // end of getFiniteElementDiscretizationParametersList

  std::vector<std::string> FiniteElementDiscretization::getParametersList() {
    auto d = MeshDiscretization::getParametersList();
    const auto names = getFiniteElementDiscretizationParametersList();
    d.insert(d.end(), names.begin(), names.end());
    return d;
  }  // end of getParametersList

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

  FiniteElementDiscretization::FiniteElementDiscretization(
      Context& ctx,
      const FiniteElementSpacesManager& m,
      const Parameters& params)
      : MeshDiscretization(m.getMeshDiscretization()), fespaces_manager(m) {
    CatchTimeSection(ctx, "FED::Constructor");
    checkParameters(
        throwing, params,
        std::vector<std::string>{FiniteElementDiscretization::UnknownsSize});
    auto or_raise = ctx.getThrowingFailureHandler();
    const auto usize =
        get<int>(throwing, params, FiniteElementDiscretization::UnknownsSize);
    if (this->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      this->parallel_fe_space =
          this->fespaces_manager.getFiniteElementSpace<true>(ctx, usize) |
          or_raise;
#else
      reportUnsupportedParallelComputations();
#endif
    } else {
      this->sequential_fe_space =
          this->fespaces_manager.getFiniteElementSpace<false>(ctx, usize) |
          or_raise;
    }
  }

  FiniteElementDiscretization::FiniteElementDiscretization(
      Context& ctx, const Parameters& params)
      : FiniteElementDiscretization(
            ctx,
            MeshDiscretization(
                ctx,
                extract(
                    throwing, params, MeshDiscretization::getParametersList())),
            extract(throwing,
                    params,
                    getFiniteElementDiscretizationParametersList())) {
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      Context& ctx, const MeshDiscretization& m, const Parameters& params)
      : FiniteElementDiscretization(
            ctx,
            FiniteElementSpacesManager(
                ctx,
                m,
                extract(throwing,
                        params,
                        FiniteElementSpacesManager::
                            getFiniteElementCollectionParametersList())),
            remove(params,
                   FiniteElementSpacesManager::
                       getFiniteElementCollectionParametersList())) {
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      Context& ctx, std::shared_ptr<Mesh<true>> m, const Parameters& params)
      : FiniteElementDiscretization(
            ctx,
            MeshDiscretization(m),
            extract(throwing,
                    params,
                    getFiniteElementDiscretizationParametersList())) {
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      Context& ctx, std::shared_ptr<Mesh<false>> m, const Parameters& params)
      : FiniteElementDiscretization(
            ctx,
            MeshDiscretization(m),
            extract(throwing,
                    params,
                    getFiniteElementDiscretizationParametersList())) {
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      Context& ctx,
      std::shared_ptr<Mesh<true>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      const size_type d)
      : FiniteElementDiscretization(
            ctx,
            FiniteElementSpacesManager(ctx, m, c),
            Parameters{{FiniteElementDiscretization::UnknownsSize, d}}) {}

  FiniteElementDiscretization::FiniteElementDiscretization(
      Context& ctx,
      std::shared_ptr<Mesh<false>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      const size_type d)
      : FiniteElementDiscretization(
            ctx,
            FiniteElementSpacesManager(ctx, m, c),
            Parameters{{FiniteElementDiscretization::UnknownsSize, d}}) {}

  const FiniteElementCollection&
  FiniteElementDiscretization::getFiniteElementCollection() const noexcept {
    return this->fespaces_manager.getFiniteElementCollection();
  }  // end of getFiniteElementCollection

  std::shared_ptr<const FiniteElementCollection>
  FiniteElementDiscretization::getFiniteElementCollectionPointer()
      const noexcept {
    return this->fespaces_manager.getFiniteElementCollectionPointer();
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