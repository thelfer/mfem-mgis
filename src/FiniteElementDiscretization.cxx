/*!
 * \file   src/FiniteElementDiscretization.cxx
 * \brief
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#include <utility>
#include <fstream>
#include <mfem/mesh/mesh.hpp>
#include <mfem/fem/fespace.hpp>
#ifdef MFEM_USE_MPI
#include <mfem/mesh/pmesh.hpp>
#include <mfem/fem/pfespace.hpp>
#endif
#include "MGIS/Raise.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

  const char* const FiniteElementDiscretization::Parallel = "Parallel";
  const char* const FiniteElementDiscretization::MeshFileName = "MeshFileName";
  const char* const FiniteElementDiscretization::FiniteElementFamily =
      "FiniteElementFamily";
  const char* const FiniteElementDiscretization::FiniteElementOrder =
      "FiniteElementOrder";
  const char* const FiniteElementDiscretization::UnknownsSize = "UnknownsSize";
  const char* const FiniteElementDiscretization::NumberOfUniformRefinements =
      "NumberOfUniformRefinements";
  const char* const FiniteElementDiscretization::GeneralVerbosityLevel =
      "GeneralVerbosityLevel";

  void FiniteElementDiscretization::reportInvalidParallelMesh() {
    raise(
        "FiniteElementDiscretization::reportInvalidParallelMesh: "
        "no parallel mesh defined");
  }  // end of reportInvalidParallelMesh

  void FiniteElementDiscretization::reportInvalidSequentialMesh() {
    raise(
        "FiniteElementDiscretization::reportInvalidSequentialMesh: "
        "no sequential mesh defined");
  }  // end of reportInvalidSequentialMesh

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

  std::vector<std::string> FiniteElementDiscretization::getParametersList() {
    return {FiniteElementDiscretization::Parallel,
            FiniteElementDiscretization::MeshFileName,
            FiniteElementDiscretization::FiniteElementFamily,
            FiniteElementDiscretization::FiniteElementOrder,
            FiniteElementDiscretization::UnknownsSize,
            FiniteElementDiscretization::NumberOfUniformRefinements,
            FiniteElementDiscretization::GeneralVerbosityLevel};
  }  // end of getParametersList

  FiniteElementDiscretization::FiniteElementDiscretization(
      const Parameters& params) {
    checkParameters(params, FiniteElementDiscretization::getParametersList());
    const auto parallel =
        get_if<bool>(params, FiniteElementDiscretization::Parallel, false);
    const auto& mesh_file =
        get<std::string>(params, FiniteElementDiscretization::MeshFileName);
    const auto& fe_family = get_if<std::string>(
        params, FiniteElementDiscretization::FiniteElementFamily, "H1");
    const auto fe_order =
        get_if<int>(params, FiniteElementDiscretization::FiniteElementOrder, 1);
    const auto u_size =
        get<int>(params, FiniteElementDiscretization::UnknownsSize);
    const auto nrefinement = get_if<int>(
        params, FiniteElementDiscretization::NumberOfUniformRefinements, 0);
    if (parallel) {
#ifdef MFEM_USE_MPI
      auto smesh = loadMeshSequential(mesh_file, 0, 1, true);
      this->parallel_mesh =
          std::make_shared<Mesh<true>>(MPI_COMM_WORLD, *smesh);
      for (size_type i = 0; i < nrefinement; ++i) {
        this->parallel_mesh->UniformRefinement();
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->sequential_mesh =
        loadMeshSequential(mesh_file, 0, 1, true);
      for (size_type i = 0; i < nrefinement; ++i) {
        this->sequential_mesh->UniformRefinement();
      }
    }
    // building the finite element collection
    if (fe_family == "H1") {
      if (parallel) {
#ifdef MFEM_USE_MPI
        this->fec = std::make_shared<mfem::H1_FECollection>(
            fe_order, this->parallel_mesh->Dimension());
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      } else {
        this->fec = std::make_shared<mfem::H1_FECollection>(
            fe_order, this->sequential_mesh->Dimension());
      }
    } else {
      raise(
          "FiniteElementDiscretization::FiniteElementDiscretization: "
          "unsupported finite element family '" +
          fe_family + "'");
    }
    // building the finite element space
    if (parallel) {
#ifdef MFEM_USE_MPI
      this->parallel_fe_space = std::make_unique<FiniteElementSpace<true>>(
          this->parallel_mesh.get(), this->fec.get(), u_size);
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      this->sequential_fe_space = std::make_unique<FiniteElementSpace<false>>(
          this->sequential_mesh.get(), this->fec.get(), u_size);
    }
  }  // end of FiniteElementDiscretization

#ifdef MFEM_USE_MPI

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<true>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      const size_type d)
      : parallel_mesh(std::move(m)), fec(std::move(c)) {
    this->parallel_fe_space = std::make_unique<FiniteElementSpace<true>>(
        this->parallel_mesh.get(), this->fec.get(), d);
  }  // end of FiniteElementDiscretization

#else /* MFEM_USE_MPI */

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<true>>,
      std::shared_ptr<const FiniteElementCollection>,
      const size_type) {
    reportUnsupportedParallelComputations();
  }  // end of FiniteElementDiscretization

#endif /* MFEM_USE_MPI */

#ifdef MFEM_USE_MPI

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<true>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      std::unique_ptr<FiniteElementSpace<true>> s)
      : parallel_mesh(std::move(m)),
        fec(std::move(c)),
        parallel_fe_space(std::move(s)) {
    if (this->parallel_mesh.get() != this->parallel_fe_space->GetMesh()) {
      raise(
          "FiniteElementDiscretization::FiniteElementDiscretization: "
          "mesh pointer don't match the mesh on which the finite element space "
          "is built");
    }
  }  // end of FiniteElementDiscretization

#else /* MFEM_USE_MPI */

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<true>>,
      std::shared_ptr<const FiniteElementCollection>,
      std::unique_ptr<FiniteElementSpace<true>>) {
    reportUnsupportedParallelComputations();
  }  // end of FiniteElementDiscretization

#endif /* MFEM_USE_MPI */

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<false>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      const size_type d)
      : sequential_mesh(std::move(m)), fec(std::move(c)) {
    this->sequential_fe_space = std::make_unique<FiniteElementSpace<false>>(
        this->sequential_mesh.get(), this->fec.get(), d);
  }  // end of FiniteElementDiscretization

  FiniteElementDiscretization::FiniteElementDiscretization(
      std::shared_ptr<Mesh<false>> m,
      std::shared_ptr<const FiniteElementCollection> c,
      std::unique_ptr<FiniteElementSpace<false>> s)
      : sequential_mesh(std::move(m)),
        fec(std::move(c)),
        sequential_fe_space(std::move(s)) {
    if (this->sequential_mesh.get() != this->sequential_fe_space->GetMesh()) {
      raise(
          "FiniteElementDiscretization::FiniteElementDiscretization: "
          "mesh pointer don't match the mesh on which the finite element space "
          "is built");
    }
  }  // end of FiniteElementDiscretization

  bool FiniteElementDiscretization::describesAParallelComputation() const {
#ifdef MFEM_USE_MPI
    return this->parallel_mesh.get() != nullptr;
#else  /* MFEM_USE_MPI */
    return false;
#endif /* MFEM_USE_MPI */
  }    // end of describesAParallelComputation

  const FiniteElementCollection&
  FiniteElementDiscretization::getFiniteElementCollection() const {
    return *(this->fec);
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

  std::shared_ptr<Mesh<false>>
  loadMeshSequential(
                     const std::string& mesh_name,
                     int generate_edges,
                     int refine,
                     bool fix_orientation) {
#ifdef MFEM_USE_MED
    auto extension = getFileExt(mesh_name);
    if (extension == "med") {
      auto medmesh = std::make_shared<Mesh<false>>();
      std::string per_name = mesh_name;
      per_name.replace (per_name.length()-4,4,".per");
      std::ifstream per_file(per_name.c_str());
      if (per_file.good()) {
        medmesh->ImportMED(mesh_name, 0, per_name);
      } else {
        medmesh->ImportMED(mesh_name, 0, "");
      } 
      //medmesh->CheckElementOrientation(fix_orientation);
      //medmesh->CheckBdrElementOrientation(fix_orientation);
      return medmesh;
    }
#endif  /* MFEM_USE_MED */
    auto smesh = std::make_shared<Mesh<false>>(mesh_name.c_str(), generate_edges, refine);
    return smesh;
    
  }  // end of loadMeshSequential
  
}  // end of namespace mfem_mgis
