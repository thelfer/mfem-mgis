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

  //! \brief remove extra spaces on the rigth
  [[nodiscard]] static std::string trim_right(const std::string& s) noexcept {
    auto r = std::string{s};
    r.erase(std::find_if(
                r.rbegin(), r.rend(),
                [](std::string::value_type ch) { return !std::isspace(ch); })
                .base(),
            r.end());
    return r;
  }

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
        << std::endl;
    std::abort();
    return std::shared_ptr<Mesh<true>>(nullptr);
#endif
  }  // end of loadMeshParallel

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

  [[nodiscard]] static size_type count(
      const std::map<size_type, std::string>& names,
      const std::string& name) noexcept {
    auto c = size_type{};
    for (const auto& [id, n] : names) {
      static_cast<void>(id);
      if (n == name) {
        ++c;
      }
    }
    return c;
  }  // end of count

  [[nodiscard]] static std::optional<size_type> key(
      const std::map<size_type, std::string>& names,
      const std::string& name) noexcept {
    for (const auto& [id, n] : names) {
      if (n == name) {
        return id;
      }
    }
    return {};
  }  // end of key

  static void setMeshObjectNames(attributes::Throwing,
                                 std::map<size_type, std::string>& ids,
                                 const std::map<size_type, std::string>& nids,
                                 const mfem::Array<size_type>& attributes,
                                 const std::string& m,
                                 const std::string& t) {
    // checks that the given identifiers are ok
    for (const auto& [a, n] : nids) {
      if (count(nids, n) != 1) {
        raise(m + ": name " + n + " multiply defined");
      }
      if (attributes.Find(a) == -1) {
        raise(m + ": no " + t + " associated with attribute '" +
              std::to_string(a) + "'");
      }
      if (!isValidMeshObjectName(n)) {
        raise(m + ": " + n + " is not a valid " + t + " identifier");
      }
      auto oa = key(ids, n);
      if (oa.has_value()) {
        if (*oa != a) {
          raise(m + ": name " + n + " is already associated to another " + t);
        }
      }
      const auto p = ids.find(a);
      if (p != ids.end()) {
        if (p->second != n) {
          warning(getDefaultLogStream(), m, ": overwritting ", t, " name '",
                  p->second, "' by '", n, "'");
        }
      }
    }
    for (const auto& [a, n] : nids) {
      const auto p = ids.find(a);
      if (p != ids.end()) {
        ids.erase(p);
      }
    }
    // declaring attributes
    ids.insert(nids.begin(), nids.end());
  }  // end of setMeshObjectNames

  template <bool parallel>
  static void updateNamesFromAttributesSets(
      attributes::Throwing,
      std::map<size_type, std::string>& materials_names,
      std::map<size_type, std::string>& boundaries_names,
      const Mesh<parallel>& mesh) {
    // checks
    for (const auto& [id, n] : materials_names) {
      if (count(materials_names, n) != 1) {
        raise("material name '" + n + "' multiply defined");
      }
      if (count(boundaries_names, n) != 0) {
        raise("material name '" + n + "' also defined as a boundary name");
      }
    }
    for (const auto& [id, n] : boundaries_names) {
      // the name can't also be a material name, we checked that in the previous
      // loop
      if (count(boundaries_names, n) != 1) {
        raise("boundary name '" + n + "' multiply defined");
      }
    }
    //
    const auto& attr_sets = mesh.attribute_sets;
    const auto& bdr_attr_sets = mesh.bdr_attribute_sets;
    const auto mnames = attr_sets.GetAttributeSetNames();
    const auto bnames = bdr_attr_sets.GetAttributeSetNames();
    for (const auto& an : attr_sets.GetAttributeSetNames()) {
      if (!attr_sets.AttributeSetExists(an)) {
        // This seems very unlikely
        continue;
      }
      if (bdr_attr_sets.AttributeSetExists(an)) {
        warning(getDefaultLogStream(), "ignoring attribute set '", an,
                "' whose name is also associated to a boundary attribute "
                "set");
        continue;
      }
      const auto& mids = attr_sets.GetAttributeSet(an);
      if (mids.Size() != 1) {
        warning(getDefaultLogStream(), "ignoring attribute set '", an,
                "' which is associated to multiple materials");
        continue;
      }
      // attributes may have extra spaces on the right.
      // At this stage, I don't know if it comes from MED convertion of GMSH,
      // but it is better to get rid of them
      const auto n = trim_right(an);
      if (count(boundaries_names, n) != 0) {
        warning(getDefaultLogStream(), "ignoring attribute set '", n,
                "' which is associated by the user to a boundary");
        continue;
      }
      if (count(materials_names, n) != 0) {
        warning(getDefaultLogStream(), "ignoring attribute set '", n,
                "' which is already associated by the user to a material");
        continue;
      }
      if (materials_names.find(mids[0]) != materials_names.end()) {
        warning(getDefaultLogStream(), "ignoring attribute set '", n,
                "' for material (", mids[0],
                ") which is already names by the user to a material");
        continue;
      }
      if (!isValidMeshObjectName(n)) {
        warning(getDefaultLogStream(), "ignoring attribute set '", n,
                "' for material (", mids[0], ") as it is not a valid name");
        continue;
      }
      materials_names.insert({mids[0], n});
    }
    for (const auto& an : bdr_attr_sets.GetAttributeSetNames()) {
      if (!bdr_attr_sets.AttributeSetExists(an)) {
        // This seems very unlikely
        continue;
      }
      if (attr_sets.AttributeSetExists(an)) {
        warning(getDefaultLogStream(), "ignoring boundary attribute set '", an,
                "' whose name is also associated to a material attribute "
                "set");
        continue;
      }
      const auto& bids = bdr_attr_sets.GetAttributeSet(an);
      if (bids.Size() != 1) {
        warning(getDefaultLogStream(), "ignoring boundary attribute set '", an,
                "' which is associated to multiple materials");
        continue;
      }
      const auto n = trim_right(an);
      if (count(boundaries_names, n) != 0) {
        warning(getDefaultLogStream(), "ignoring boundary attribute set '", n,
                "' which is associated by the user to a boundary");
        continue;
      }
      if (count(materials_names, n) != 0) {
        warning(getDefaultLogStream(), "ignoring boundary attribute set '", n,
                "' which is already associated by the user to a material");
        continue;
      }
      if (boundaries_names.find(bids[0]) != boundaries_names.end()) {
        warning(getDefaultLogStream(), "ignoring boundary attribute set '", n,
                "' for material (", bids[0],
                ") which is already names by the user to a material");
        continue;
      }
      if (!isValidMeshObjectName(n)) {
        warning(getDefaultLogStream(), "ignoring attribute set '", n,
                "' for boundary (", bids[0], ") as it is not a valid name");
        continue;
      }
      boundaries_names.insert({bids[0], n});
    }
  }  // end of updateNamesFromAttributesSets

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
    auto mnames = [&params, extractMap]() -> std::map<size_type, std::string> {
      if (contains(params, FiniteElementDiscretization::Materials)) {
        return extractMap(
            get<Parameters>(params, FiniteElementDiscretization::Materials));
      }
      return {};
    }();
    auto bnames = [&params, extractMap]() -> std::map<size_type, std::string> {
      if (contains(params, FiniteElementDiscretization::Boundaries)) {
        return extractMap(
            get<Parameters>(params, FiniteElementDiscretization::Boundaries));
      }
      return {};
    }();
    if (parallel) {
      updateNamesFromAttributesSets<true>(throwing, mnames, bnames,
                                          *(this->parallel_mesh));
    } else {
      updateNamesFromAttributesSets<false>(throwing, mnames, bnames,
                                           *(this->sequential_mesh));
    }
    if (!mnames.empty()) {
      this->setMaterialsNames(throwing, mnames);
    }
    if (!bnames.empty()) {
      this->setBoundariesNames(throwing, bnames);
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

  bool FiniteElementDiscretization::setMaterialsNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    try {
      this->setMaterialsNames(throwing, ids);
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of setMaterialsNames

  bool FiniteElementDiscretization::setBoundariesNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    try {
      this->setBoundariesNames(throwing, ids);
    } catch (...) {
      return registerExceptionInErrorBacktrace(ctx);
    }
    return true;
  }  // end of setBoundariesNames

  void FiniteElementDiscretization::setMaterialsNames(
      const std::map<size_type, std::string>& ids) {
    this->setMaterialsNames(throwing, ids);
  }  // end of setMaterialsNames

  void FiniteElementDiscretization::setBoundariesNames(
      const std::map<size_type, std::string>& ids) {
    this->setBoundariesNames(throwing, ids);
  }  // end of setBoundariesNames

  void FiniteElementDiscretization::setMaterialsNames(
      attributes::Throwing, const std::map<size_type, std::string>& ids) {
    setMeshObjectNames(throwing, this->materials_names, ids,
                       getMaterialsAttributes(*this), "setMaterialsNames",
                       "material");
  }  // end of setMaterialsNames

  void FiniteElementDiscretization::setBoundariesNames(
      attributes::Throwing, const std::map<size_type, std::string>& ids) {
    setMeshObjectNames(throwing, this->boundaries_names, ids,
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

  [[nodiscard]] static std::vector<size_type> selectMeshObjectsIdentifiers(
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
    } catch (std::exception&) {
      raise(m + ": invalid regular expression '" + id + "'");
    }
    return r;
  }  // end of selectMeshObjectsIdentifiers

  [[nodiscard]] static std::vector<size_type> selectMeshObjectsIdentifiers(
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

  std::optional<std::vector<size_type>>
  FiniteElementDiscretization::getMaterialsIdentifiers(
      Context& ctx, const Parameter& p) const noexcept {
    try {
      return this->getMaterialsIdentifiers(p);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of getMaterialsIdentifiers

  std::optional<std::vector<size_type>>
  FiniteElementDiscretization::getBoundariesIdentifiers(
      Context& ctx, const Parameter& p) const noexcept {
    try {
      return this->getBoundariesIdentifiers(p);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of getBoundariesIdentifiers

  std::optional<size_type> FiniteElementDiscretization::getMaterialIdentifier(
      Context& ctx, const Parameter& p) const noexcept {
    try {
      return this->getMaterialIdentifier(p);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of getMaterialIdentifier

  std::optional<size_type> FiniteElementDiscretization::getBoundaryIdentifier(
      Context& ctx, const Parameter& p) const noexcept {
    try {
      return this->getBoundaryIdentifier(p);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
    }
    return {};
  }  // end of getBoundaryIdentifier

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

  std::map<size_type, std::string>
  FiniteElementDiscretization::getMaterialsNames() const noexcept {
    return this->materials_names;
  }  // end of getMaterialsNames

  std::map<size_type, std::string>
  FiniteElementDiscretization::getBoundariesNames() const noexcept {
    return this->boundaries_names;
  }  // end of getBoundariesNames

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

  template <>
  bool getInformation<FiniteElementDiscretization>(
      Context&,
      std::ostream& os,
      const FiniteElementDiscretization& fed) noexcept {
    const auto& mnames = fed.getMaterialsNames();
    os << "# Mesh\n\n"
       << "- space dimension: " << getSpaceDimension(fed);
    if (!mnames.empty()) {
      os << "\n\n## Materials\n";
      for (const auto& [id, n] : mnames) {
        os << "\n- '" << n << "' associated with identifier (" << id << ")";
      }
    }
    const auto& bnames = fed.getBoundariesNames();
    if (!bnames.empty()) {
      os << "\n\n## Boundaries\n";
      for (const auto& [id, n] : bnames) {
        os << "\n- '" << n << "' associated with identifier (" << id << ")";
      }
    }
    os << "\n\n# Finite element space\n\n"
       << "- true vector size: " << getTrueVSize(fed) << '\n';
    return true;
  }  // end of info

}  // end of namespace mfem_mgis
