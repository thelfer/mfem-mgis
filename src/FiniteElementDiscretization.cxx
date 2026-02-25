/*!
 * \file   src/FiniteElementDiscretization.cxx
 * \brief
 * \author Thomas Helfer
 * \date 16/12/2020
 */

#include <iostream>
#include <cctype>
#include <utility>
#include <regex>
#include <fstream>
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

#ifdef MFEM_USE_MED

  /*!
   * \brief Extract the file extension
   * \param[in] s: string corresponding to a file name
   */
  static std::string getFileExt(const std::string& s) {
    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos) {
      return (s.substr(i + 1, s.length() - i));
    }

    return ("");
  }  // end of getFileExt

#endif /* MFEM_USE_MED */

  /*!
   * \brief load a mesh (sequential)
   * \param[in] s: string corresponding to a file name
   *
   * \note MED format is handled in addition to standard MFEM
   *       input formats.
   */
  static std::shared_ptr<Mesh<false>> loadMeshSequential(
      const std::string& mesh_name,
      int generate_edges = 0,
      int refine = 1,
      bool /* fix_orientation */ = true) {
    CatchTimeSection("FED::LoadMesh");
#ifdef MFEM_USE_MED
    const auto extension = getFileExt(mesh_name);
    if (extension == "med") {
      auto medmesh = std::make_shared<Mesh<false>>();
      std::string per_name = mesh_name;
      per_name.replace(per_name.length() - 4, 4, ".per");
      std::ifstream per_file(per_name.c_str());
      if (per_file.good()) {
        medmesh->ImportMED(mesh_name, 0, per_name);
      } else {
        medmesh->ImportMED(mesh_name, 0, "");
      }
      // medmesh->CheckElementOrientation(fix_orientation);
      // medmesh->CheckBdrElementOrientation(fix_orientation);
      return medmesh;
    }
#endif /* MFEM_USE_MED */
    auto smesh = std::make_shared<Mesh<false>>(mesh_name.c_str(),
                                               generate_edges, refine);
    return smesh;
  }  // end of loadMeshSequential

  static std::shared_ptr<Mesh<true>> loadMeshParallel(
      const std::string& mesh_name) {
    CatchTimeSection("FED::LoadMeshInParallel");
#ifdef MFEM_USE_MED
    const auto extension = getFileExt(mesh_name);
    if (extension == "med") {
      std::cout << "Aborting. The option '--restart' is not handled with a "
                   "'.med' file format"
                << std::endl;
      std::abort();
    }
#endif /* MFEM_USE_MED */
#ifdef MFEM_USE_MPI
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    std::string fname(mfem::MakeParFilename(mesh_name, myid));
    std::ifstream ifs(fname);
    MFEM_VERIFY(ifs.good(), "Checkpoint file " << fname << " not found.");
    auto smesh = std::make_shared<Mesh<true>>(Mesh<true>(MPI_COMM_WORLD, ifs));
    //    auto smesh = std::make_shared<Mesh<true>>(mesh_name.c_str(),
    //                                               generate_edges, refine);
    return smesh;
#else
    std::cout
        << "Aborting. The option Restart is not compatible with sequential run."
        << "Enable to load mesh file '" << mesh_name << "'\n"
        << std::endl;
    std::abort();
#endif
  }  // end of loadMeshParallel

  const char* const FiniteElementDiscretization::Parallel = "Parallel";
  const char* const FiniteElementDiscretization::MeshFileName = "MeshFileName";
  const char* const FiniteElementDiscretization::MeshReadMode = "MeshReadMode";
  const char* const FiniteElementDiscretization::FiniteElementFamily =
      "FiniteElementFamily";
  const char* const FiniteElementDiscretization::FiniteElementOrder =
      "FiniteElementOrder";
  const char* const FiniteElementDiscretization::Materials = "Materials";
  const char* const FiniteElementDiscretization::Boundaries = "Boundaries";
  const char* const FiniteElementDiscretization::UnknownsSize = "UnknownsSize";
  const char* const FiniteElementDiscretization::NumberOfUniformRefinements =
      "NumberOfUniformRefinements";
  const char* const FiniteElementDiscretization::GeneralVerbosityLevel =
      "GeneralVerbosityLevel";

  const mfem::Array<size_type>& getMaterialsAttributes(
      const FiniteElementDiscretization& fed) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& mesh = fed.getMesh<true>();
      return mesh.attributes;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    const auto& mesh = fed.getMesh<false>();
    return mesh.attributes;
  }  // end of getMaterialsAttributes

  const mfem::Array<size_type>& getBoundariesAttributes(
      const FiniteElementDiscretization& fed) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& mesh = fed.getMesh<true>();
      return mesh.bdr_attributes;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    const auto& mesh = fed.getMesh<false>();
    return mesh.bdr_attributes;
  }  // end of getBoundariesAttributes

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
            FiniteElementDiscretization::MeshReadMode,
            FiniteElementDiscretization::FiniteElementFamily,
            FiniteElementDiscretization::FiniteElementOrder,
            FiniteElementDiscretization::UnknownsSize,
            FiniteElementDiscretization::NumberOfUniformRefinements,
            FiniteElementDiscretization::Materials,
            FiniteElementDiscretization::Boundaries,
            FiniteElementDiscretization::GeneralVerbosityLevel};
  }  // end of getParametersList

  FiniteElementDiscretization::FiniteElementDiscretization(
      const Parameters& params) {
    CatchTimeSection("FED::Constructor");
    auto extractMap = [](const Parameters& parameters) {
      auto m = std::map<size_type, std::string>{};
      for (const auto& p : parameters) {
        m[get<int>(p.second)] = p.first;
      }
      return m;
    };
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
    const auto mesh_mode = get_if<std::string>(
        params, FiniteElementDiscretization::MeshReadMode, "FromScratch");
    if (parallel) {
#ifdef MFEM_USE_MPI

      size_type ref_level = 0;
      if (mesh_mode == "FromScratch") {
        auto smesh = loadMeshSequential(mesh_file, 0, 1, true);
        // Perform a uniform refinement on the sequential mesh if it doesn't
        // have enough elements. Assume that each subdomain should have at leat
        // 8 elements Not superior to nrefinement
        if (nrefinement > 0) {
          double numberOfProcs = double(mfem::Mpi::WorldSize());
          while ((double(smesh->GetNE()) / numberOfProcs) < 8 &&
                 ref_level < nrefinement) {
            smesh->UniformRefinement();
            ref_level++;
          }
        }
        this->parallel_mesh =
            std::make_shared<Mesh<true>>(MPI_COMM_WORLD, *smesh);
      } else if (mesh_mode == "Restart") {
        this->parallel_mesh = loadMeshParallel(mesh_file);
      } else {
        raise("Wrong MeshReadMode value");
      }
      for (size_type i = ref_level; i < nrefinement; ++i) {
        CatchNestedTimeSection("FED::Run_ParUniformRefinement");
        this->parallel_mesh->UniformRefinement();
      }
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    } else {
      if (mesh_mode == "Restart") {
        raise(
            "Aborting. The option 'Restart' is not handled while running a "
            "sequential program");
      }
      this->sequential_mesh = loadMeshSequential(mesh_file, 0, 1, true);
      for (size_type i = 0; i < nrefinement; ++i) {
        CatchNestedTimeSection("FED::Run_SeqUniformRefinement");
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
    // declaring materials and boundaries
    if (contains(params, FiniteElementDiscretization::Materials)) {
      this->setMaterialsNames(extractMap(
          get<Parameters>(params, FiniteElementDiscretization::Materials)));
    }
    if (contains(params, FiniteElementDiscretization::Boundaries)) {
      this->setBoundariesNames(extractMap(
          get<Parameters>(params, FiniteElementDiscretization::Boundaries)));
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
  FiniteElementDiscretization::getFiniteElementCollection() const noexcept {
    return *(this->fec);
  }  // end of getFiniteElementCollection

  std::shared_ptr<const FiniteElementCollection>
  FiniteElementDiscretization::getFiniteElementCollectionPointer()
      const noexcept {
    return this->fec;
  }  // end of getFiniteElementCollection

  static bool isValidMeshObjectName(const std::string& n) {
    if (n.empty()) {
      return false;
    }
    auto p = n.begin();
    if (std::isdigit(*p)) {
      return false;
    }
    for (; p != n.end(); ++p) {
      if ((!std::isalpha(*p)) && (!(std::isdigit(*p))) && (*p != '_')) {
        return false;
      }
      if (std::isspace(*p)) {
        return false;
      }
    }
    return true;
  }  // end of isValidMeshObjectName

  static void setMeshObjectNames(std::map<size_type, std::string>& ids,
                                 const std::map<size_type, std::string>& nids,
                                 const mfem::Array<size_type>& attributes,
                                 const std::string& m,
                                 const std::string& t) {
    // checks that the given identifiers are ok
    for (const auto& [a, n] : nids) {
      if (attributes.Find(a) == -1) {
        raise(m + ": no " + t + " associated with attribute '" +
              std::to_string(a) + "'");
      }
      if (!isValidMeshObjectName(n)) {
        raise(m + ": " + n + " is not a valid " + t + " identifier");
      }
      const auto p = ids.find(a);
      if (p != ids.end()) {
        raise(m + ": a name has already been associated with " + t + " " +
              std::to_string(a) + " ('" + p->second + "') ");
      }
    }
    // declaring attributes
    ids.insert(nids.begin(), nids.end());
  }  // end of setMeshObjectNames

  void FiniteElementDiscretization::setMaterialsNames(
      const std::map<size_type, std::string>& ids) {
    setMeshObjectNames(this->materials_names, ids,
                       getMaterialsAttributes(*this), "setMaterialsNames",
                       "material");
  }  // end of setMaterialsNames

  void FiniteElementDiscretization::setBoundariesNames(
      const std::map<size_type, std::string>& ids) {
    setMeshObjectNames(this->boundaries_names, ids,
                       getBoundariesAttributes(*this), "setBoundariesNames",
                       "boundary");
  }  // end of setBoundariesNames

  std::optional<std::string> FiniteElementDiscretization::getMaterialName(
      Context& ctx, const size_type id) const noexcept {
    const auto& mids = getMaterialsAttributes(*this);
    if (mids.Find(id) == -1) {
      return ctx.registerErrorMessage("no material id '" + std::to_string(id) +
                                      "' defined in the mesh");
    }
    const auto p = this->materials_names.find(id);
    if (p == this->materials_names.end()) {
      return std::string{};
    }
    return p->second;
  }  // end of getMaterialName

  std::optional<std::string> FiniteElementDiscretization::getBoundaryName(
      Context& ctx, const size_type id) const noexcept {
    const auto& bids = getBoundariesAttributes(*this);
    if (bids.Find(id) == -1) {
      return ctx.registerErrorMessage("no boundary id '" + std::to_string(id) +
                                      "' defined in the mesh");
    }
    const auto p = this->boundaries_names.find(id);
    if (p == this->boundaries_names.end()) {
      return std::string{};
    }
    return p->second;
  }  // end of getBoundaryName

  static std::vector<size_type> selectMeshObjectsIdentifiers(
      const mfem::Array<size_type>& attributes,
      const size_type id,
      const std::string& t,
      const std::string& m) {
    if (attributes.Find(id) == -1) {
      raise(m + ": no " + t + " associated with attribute '" +
            std::to_string(id) + "'");
    }
    return {id};
  }  // end of selectMeshObjectsIdentifiers

  std::vector<size_type> selectMeshObjectsIdentifiers(
      const std::map<size_type, std::string>& names,
      const std::string& id,
      const std::string& t,
      const std::string& m) {
    auto r = std::vector<size_type>{};
    try {
      std::regex e(id);
      for (const auto& [a, n] : names) {
        if (std::regex_match(n, e)) {
          r.push_back(a);
        }
      }
      if (r.empty()) {
        raise(m + ": no " + t + " matching regular expression '" + id + "'");
      }
    } catch (std::exception& e) {
      raise(m + ": invalid regular expression '" + id + "'");
    }
    return r;
  }  // end of selectMeshObjectsIdentifiers

  std::vector<size_type> selectMeshObjectsIdentifiers(
      const mfem::Array<size_type>& attributes,
      const std::map<size_type, std::string>& names,
      const std::vector<Parameter>& ids,
      const std::string& t,
      const std::string& m) {
    if (ids.empty()) {
      raise(m + ": empty list of identifiers");
    }
    auto r = std::vector<size_type>{};
    auto append = [&r, &m](const auto& nids) {
      for (const auto& id : nids) {
        if (std::find(std::cbegin(r), std::cend(r), id) != std::cend(r)) {
          raise(m + ": identifier '" + std::to_string(id) +
                "' multiply selected");
        }
        r.push_back(id);
      }
    };
    for (const auto& id : ids) {
      if (is<size_type>(id)) {
        const auto i = get<size_type>(id);
        append(selectMeshObjectsIdentifiers(attributes, i, t, m));
      } else if (is<std::string>(id)) {
        const auto& n = get<std::string>(id);
        append(selectMeshObjectsIdentifiers(names, n, t, m));
      } else {
        raise(m + ": invalid parameter");
      }
    }
    return r;
  }  // end of selectMeshObjectsIdentifiers

  std::vector<size_type> selectMeshObjectsIdentifiers(
      const mfem::Array<size_type>& attributes,
      const std::map<size_type, std::string>& names,
      const Parameter& p,
      const std::string& t,
      const std::string& m) {
    if (is<size_type>(p)) {
      const auto id = get<size_type>(p);
      return selectMeshObjectsIdentifiers(attributes, id, t, m);
    } else if (is<std::string>(p)) {
      const auto& id = get<std::string>(p);
      return selectMeshObjectsIdentifiers(names, id, t, m);
    }
    if (!is<std::vector<Parameter>>(p)) {
      raise(m + ": invalid parameter type");
    }
    const auto& ids = get<std::vector<Parameter>>(p);
    return selectMeshObjectsIdentifiers(attributes, names, ids, t, m);
  }  // end of selectMeshObjectsIdentifiers

  std::vector<size_type> FiniteElementDiscretization::getMaterialsIdentifiers(
      const Parameter& p) const {
    return selectMeshObjectsIdentifiers(getMaterialsAttributes(*this),
                                        this->materials_names, p, "material",
                                        "getMaterialsIdentifiers");
  }  // end of getMaterialsIdentifiers

  std::vector<size_type> FiniteElementDiscretization::getBoundariesIdentifiers(
      const Parameter& p) const {
    return selectMeshObjectsIdentifiers(getBoundariesAttributes(*this),
                                        this->boundaries_names, p, "boundary",
                                        "getBoundariesIdentifiers");
  }  // end of getBoundariesIdentifiers

  size_type FiniteElementDiscretization::getMaterialIdentifier(
      const Parameter& p) const {
    if (is<size_type>(p)) {
      const auto id = get<size_type>(p);
      const auto ids = getMaterialsAttributes(*this);
      if (ids.Find(id) == -1) {
        raise(
            "getMaterialIdentifier: "
            "no material id for identifier '" +
            std::to_string(id) + "'");
      }
      return id;
    }
    if (!is<std::string>(p)) {
      raise("getMaterialIdentifier: invalid parameter type");
    }
    const auto& n = get<std::string>(p);
    for (const auto& [id, name] : this->materials_names) {
      if (name == n) {
        return id;
      }
    }
    raise("getMaterialIdentifier: no material named '" + n + "'");
  }  // end of getMaterialIdentifier

  size_type FiniteElementDiscretization::getBoundaryIdentifier(
      const Parameter& p) const {
    if (is<size_type>(p)) {
      const auto id = get<size_type>(p);
      const auto ids = getBoundariesAttributes(*this);
      if (ids.Find(id) == -1) {
        raise(
            "getBoundaryIdentifier: "
            "no boundary id associated with identifier '" +
            std::to_string(id) + "'");
      }
      return id;
    }
    if (!is<std::string>(p)) {
      raise("getBoundaryIdentifier: invalid parameter type");
    }
    const auto& n = get<std::string>(p);
    for (const auto& [id, name] : this->boundaries_names) {
      if (name == n) {
        return id;
      }
    }
    raise("getMaterialIdentifier: no boundary named '" + n + "'");
  }  // end of getBoundaryIdentifier

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

  size_type getSpaceDimension(const FiniteElementDiscretization& fed) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return fed.getMesh<true>().SpaceDimension();
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return fed.getMesh<false>().SpaceDimension();
  }  // end of getSpaceDimension

}  // end of namespace mfem_mgis
