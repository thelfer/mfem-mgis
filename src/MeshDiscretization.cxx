/*!
 * \file   src/MeshDiscretization.cxx
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
#include <mfem/general/error.hpp>
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
#include "MFEMMGIS/MeshDiscretization.hxx"

namespace mfem_mgis {

  //! \brief remove extra spaces on the right
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
   * input formats.
   */
  static std::shared_ptr<Mesh<false>> loadMeshSequential(
      mgis::Context& ctx,
      const std::string& mesh_name,
      int generate_edges = 0,
      int refine = 1,
      bool /* fix_orientation */ = true) {
    CatchTimeSection(ctx, "Mesh::LoadMesh");
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
      mgis::Context& ctx, const std::string& mesh_name) {
    CatchTimeSection(ctx, "Mesh::LoadMeshInParallel");
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

  [[nodiscard]] static bool setMeshObjectNames(
      Context& ctx,
      std::map<size_type, std::string>& ids,
      const std::map<size_type, std::string>& nids,
      const mfem::Array<size_type>& attributes,
      const std::string& m,
      const std::string& t) noexcept {
    // checks that the given identifiers are ok
    for (const auto& [a, n] : nids) {
      if (count(nids, n) != 1) {
        return ctx.registerErrorMessage(m + ": name " + n +
                                        " multiply defined");
      }
      if (attributes.Find(a) == -1) {
        return ctx.registerErrorMessage(m + ": no " + t +
                                        " associated with attribute '" +
                                        std::to_string(a) + "'");
      }
      if (!isValidMeshObjectName(n)) {
        return ctx.registerErrorMessage(m + ": " + n + " is not a valid " + t +
                                        " identifier");
      }
      auto oa = key(ids, n);
      if (oa.has_value()) {
        if (*oa != a) {
          return ctx.registerErrorMessage(
              m + ": name " + n + " is already associated to another " + t);
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
    return true;
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
        warning(getDefaultLogStream(), "ignoring material attribute set '", n,
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

  const char* const MeshDiscretization::Parallel = "Parallel";
  const char* const MeshDiscretization::MeshFileName = "MeshFileName";
  const char* const MeshDiscretization::MeshReadMode = "MeshReadMode";
  const char* const MeshDiscretization::Materials = "Materials";
  const char* const MeshDiscretization::Boundaries = "Boundaries";
  const char* const MeshDiscretization::Points = "Points";
  const char* const MeshDiscretization::PointsSets = "PointsSets";
  const char* const MeshDiscretization::NumberOfUniformRefinements =
      "NumberOfUniformRefinements";
  const char* const MeshDiscretization::GeneralVerbosityLevel =
      "GeneralVerbosityLevel";

  [[nodiscard]] static std::optional<std::vector<size_type>>
  selectMeshObjectsIdentifiers(Context& ctx,
                               const mfem::Array<size_type>& attributes,
                               const size_type id,
                               const std::string& t,
                               const std::string& m) noexcept {
    if (attributes.Find(id) == -1) {
      return ctx.registerErrorMessage(m + ": no " + t +
                                      " associated with attribute '" +
                                      std::to_string(id) + "'");
    }
    return std::vector<size_type>{id};
  }  // end of selectMeshObjectsIdentifiers

  [[nodiscard]] static std::optional<std::vector<size_type>>
  selectMeshObjectsIdentifiers(Context& ctx,
                               const std::map<size_type, std::string>& names,
                               const std::string& id,
                               const std::string& t,
                               const std::string& m) noexcept {
    auto r = std::vector<size_type>{};
    try {
      std::regex e(id);
      for (const auto& [a, n] : names) {
        if (std::regex_match(n, e)) {
          r.push_back(a);
        }
      }
      if (r.empty()) {
        return ctx.registerErrorMessage(
            m + ": no " + t + " matching regular expression '" + id + "'");
      }
    } catch (std::exception&) {
      return ctx.registerErrorMessage(m + ": invalid regular expression '" +
                                      id + "'");
    }
    return r;
  }  // end of selectMeshObjectsIdentifiers

  [[nodiscard]] static std::optional<std::vector<size_type>>
  selectMeshObjectsIdentifiers(Context& ctx,
                               const mfem::Array<size_type>& attributes,
                               const std::map<size_type, std::string>& names,
                               const std::vector<Parameter>& ids,
                               const std::string& t,
                               const std::string& m) noexcept {
    if (ids.empty()) {
      return ctx.registerErrorMessage(m + ": empty list of identifiers");
    }
    auto r = std::vector<size_type>{};
    auto append = [&ctx, &r, &m](const auto& nids) -> bool {
      for (const auto& id : nids) {
        if (std::find(std::cbegin(r), std::cend(r), id) != std::cend(r)) {
          return ctx.registerErrorMessage(m + ": identifier '" +
                                          std::to_string(id) +
                                          "' multiply selected");
        }
        r.push_back(id);
      }
      return true;
    };
    for (const auto& id : ids) {
      if (is<size_type>(id)) {
        const auto i = get<size_type>(throwing, id);
        const auto oids =
            selectMeshObjectsIdentifiers(ctx, attributes, i, t, m);
        if (isInvalid(oids)) {
          return {};
        }
        if (!append(*oids)) {
          return {};
        }
      } else if (is<std::string>(id)) {
        const auto& n = get<std::string>(throwing, id);
        const auto oids = selectMeshObjectsIdentifiers(ctx, names, n, t, m);
        if (isInvalid(oids)) {
          return {};
        }
        if (!append(*oids)) {
          return {};
        }
      } else {
        return ctx.registerErrorMessage(m + ": invalid parameter");
      }
    }
    return r;
  }  // end of selectMeshObjectsIdentifiers

  [[nodiscard]] static std::optional<std::vector<size_type>>
  selectMeshObjectsIdentifiers(Context& ctx,
                               const mfem::Array<size_type>& attributes,
                               const std::map<size_type, std::string>& names,
                               const Parameter& p,
                               const std::string& t,
                               const std::string& m) noexcept {
    if (is<size_type>(p)) {
      const auto id = get<size_type>(throwing, p);
      return selectMeshObjectsIdentifiers(ctx, attributes, id, t, m);
    } else if (is<std::string>(p)) {
      const auto& id = get<std::string>(throwing, p);
      return selectMeshObjectsIdentifiers(ctx, names, id, t, m);
    }
    if (!is<std::vector<Parameter>>(p)) {
      return ctx.registerErrorMessage(m + ": invalid parameter type");
    }
    const auto& ids = get<std::vector<Parameter>>(throwing, p);
    return selectMeshObjectsIdentifiers(ctx, attributes, names, ids, t, m);
  }  // end of selectMeshObjectsIdentifiers

  //! \brief Implementation class for MeshDiscretization using the PIMPL idiom
  struct MeshDiscretization::Implementation {
    /*!
     * \brief report an error when no parallel mesh is defined
     * \note This function never returns and aborts the computation
     */
    [[noreturn]] static void reportInvalidParallelMesh() noexcept {
      mfem::mfem_error(
          "MeshDiscretization::reportInvalidParallelMesh: "
          "no parallel mesh defined");
    }  // end of reportInvalidParallelMesh

    /*!
     * \brief report an error when no sequential mesh is defined
     * \note This function never returns and aborts the computation
     */
    [[noreturn]] static void reportInvalidSequentialMesh() noexcept {
      mfem::mfem_error(
          "MeshDiscretization::reportInvalidSequentialMesh: "
          "no sequential mesh defined");
    }  // end of reportInvalidSequentialMesh

    /*!
     * \brief return the space dimension of the mesh
     * \param[in] m: mesh discretization implementation
     * \return the space dimension
     */
    [[nodiscard]] static size_type getSpaceDimension(
        const Implementation& m) noexcept {
      if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
        return m.getMeshPointer<true>()->SpaceDimension();
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      }
      return m.getMeshPointer<false>()->SpaceDimension();
    }

    /*!
     * \brief return the materials attributes
     * \param[in] m: mesh discretization implementation
     * \return the materials attributes
     */
    static const mfem::Array<size_type>& getMaterialsAttributes(
        const Implementation& m) noexcept {
      if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
        const auto& mesh = *(m.getMeshPointer<true>());
        return mesh.attributes;
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      }
      const auto& mesh = *(m.getMeshPointer<false>());
      return mesh.attributes;
    }  // end of getMaterialsAttributes

    /*!
     * \brief return the boundaries attributes
     * \param[in] m: mesh discretization implementation
     * \return the boundaries attributes
     */
    static const mfem::Array<size_type>& getBoundariesAttributes(
        const Implementation& m) noexcept {
      if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
        const auto& mesh = *(m.getMeshPointer<true>());
        return mesh.bdr_attributes;
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      }
      const auto& mesh = *(m.getMeshPointer<false>());
      return mesh.bdr_attributes;
    }  // end of getBoundariesAttributes

#ifdef MGIS_HAVE_TFEL
    /*!
     * \brief add points from parameters
     * \param[in] throwing: throwing attributes
     * \param[in, out] m: mesh discretization implementation
     * \param[in] parameters: parameters containing point definitions
     */
    static void addPoints(attributes::Throwing,
                          MeshDiscretization::Implementation& m,
                          const Parameters& parameters) {
      auto ctx = Context{};
      auto or_raise = ctx.getThrowingFailureHandler();
      const auto d = getSpaceDimension(m);
      if ((d != 2) && (d != 3)) {
        raise("can only add points in 2D or 3D");
      }
      for (const auto& [n, p] : parameters) {
        if (d == 2) {
          const auto pt = makePoint<2>(ctx, p) | or_raise;
          m.addPoint(ctx, n, pt) | or_raise;
        } else {
          const auto pt = makePoint<3>(ctx, p) | or_raise;
          m.addPoint(ctx, n, pt) | or_raise;
        }
      }
    }  // end of addPoints

    /*!
     * \brief add points sets from parameters
     * \param[in] throwing: throwing attributes
     * \param[in, out] m: mesh discretization implementation
     * \param[in] parameters: parameters containing points set definitions
     */
    static void addPointsSets(attributes::Throwing,
                              MeshDiscretization::Implementation& m,
                              const Parameters& parameters) {
      auto ctx = Context{};
      auto or_raise = ctx.getThrowingFailureHandler();
      const auto d = getSpaceDimension(m);
      if ((d != 2) && (d != 3)) {
        raise("can only add points sets in 2D or 3D");
      }
      for (const auto& [n, p] : parameters) {
        if (d == 2) {
          const auto pts =
              makePointsSet<2>(ctx, m.pointsSets2D, m.points2D, p) | or_raise;
          m.addPointsSet(ctx, n, pts) | or_raise;
        } else {
          const auto pts =
              makePointsSet<3>(ctx, m.pointsSets3D, m.points3D, p) | or_raise;
          m.addPointsSet(ctx, n, pts) | or_raise;
        }
      }
    }  // end of addPointsSet
#endif /* MGIS_HAVE_TFEL */

    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context used for profiling
     * \param[in] params: parameters
     *
     * The following parameters are expected:
     *
     * - `Parallel` (boolean): if true, a parallel computation is to be
     * performed. This value is assumed to be false by default.
     * - `MeshFileName` (string): mesh file.
     * - `NumberOfUniformRefinements` (int): number of uniform refinements
     *   applied to the mesh
     * - `MeshReadMode` (string): how to read the mesh. Supported values are
     *   "FromScratch" and "Restart".
     * - `GeneralVerbosityLevel` (int): with large positive numbers, expect more
     *   verbosity
     * - `Materials` (map): mapping between material identifiers and names
     * - `Boundaries` (map): mapping between boundary identifiers and names
     * - `Points` (map): points to be added to the mesh
     * - `PointsSets` (map): points sets to be added to the mesh
     */
    Implementation(mgis::Context& ctx, const Parameters& params) {
      CatchTimeSection(ctx, "Mesh::Constructor");
      auto or_raise = ctx.getThrowingFailureHandler();
      auto extractMap = [](const Parameters& parameters) {
        auto m = std::map<size_type, std::string>{};
        for (const auto& p : parameters) {
          m[get<int>(throwing, p.second)] = p.first;
        }
        return m;
      };
      checkParameters(throwing, params,
                      MeshDiscretization::getParametersList());
      const auto parallel =
          get_if<bool>(throwing, params, MeshDiscretization::Parallel, false);
      const auto& mesh_file =
          get<std::string>(throwing, params, MeshDiscretization::MeshFileName);
      const auto nrefinement = get_if<int>(
          throwing, params, MeshDiscretization::NumberOfUniformRefinements, 0);
      const auto mesh_mode = get_if<std::string>(
          throwing, params, MeshDiscretization::MeshReadMode, "FromScratch");
      if (parallel) {
#ifdef MFEM_USE_MPI
        size_type ref_level = 0;
        if (mesh_mode == "FromScratch") {
          auto smesh = loadMeshSequential(ctx, mesh_file, 0, 1, true);
          // Perform a uniform refinement on the sequential mesh if it doesn't
          // have enough elements. Assume that each subdomain should have at
          // leat 8 elements Not superior to nrefinement
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
          this->parallel_mesh = loadMeshParallel(ctx, mesh_file);
        } else {
          raise("Wrong MeshReadMode value");
        }
        for (size_type i = ref_level; i < nrefinement; ++i) {
          CatchTimeSection(ctx, "Mesh::Run_ParUniformRefinement");
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
        this->sequential_mesh = loadMeshSequential(ctx, mesh_file, 0, 1, true);
        for (size_type i = 0; i < nrefinement; ++i) {
          CatchTimeSection(ctx, "Mesh::Run_SeqUniformRefinement");
          this->sequential_mesh->UniformRefinement();
        }
      }
      // building the finite element collection
      // declaring materials and boundaries
      auto mnames = [&params,
                     extractMap]() -> std::map<size_type, std::string> {
        if (contains(params, MeshDiscretization::Materials)) {
          return extractMap(
              get<Parameters>(throwing, params, MeshDiscretization::Materials));
        }
        return {};
      }();
      auto bnames = [&params,
                     extractMap]() -> std::map<size_type, std::string> {
        if (contains(params, MeshDiscretization::Boundaries)) {
          return extractMap(get<Parameters>(throwing, params,
                                            MeshDiscretization::Boundaries));
        }
        return {};
      }();
      if (parallel) {
#ifdef MFEM_USE_MPI
        updateNamesFromAttributesSets<true>(throwing, mnames, bnames,
                                            *(this->parallel_mesh));
#else
        reportUnsupportedParallelComputations();
#endif
      } else {
        updateNamesFromAttributesSets<false>(throwing, mnames, bnames,
                                             *(this->sequential_mesh));
      }
      if (!mnames.empty()) {
        this->setMaterialsNames(ctx, mnames) | or_raise;
      }
      if (!bnames.empty()) {
        this->setBoundariesNames(ctx, bnames) | or_raise;
      }
      //
      if (contains(params, MeshDiscretization::Points)) {
        addPoints(
            throwing, *this,
            get<Parameters>(throwing, params, MeshDiscretization::Points));
      }
      //
      if (contains(params, MeshDiscretization::PointsSets)) {
        addPointsSets(
            throwing, *this,
            get<Parameters>(throwing, params, MeshDiscretization::PointsSets));
      }
    }  // end of Implementation

#ifdef MFEM_USE_MPI

    /*!
     * \brief constructor
     * \param[in] m: parallel mesh
     */
    Implementation(std::shared_ptr<Mesh<true>> m)
        : parallel_mesh(std::move(m)) {
      if (this->parallel_mesh.get() == nullptr) {
        raise("invalid mesh");
      }
    }  // end of Implementation

#else /* MFEM_USE_MPI */

    /*!
     * \brief constructor (not supported in sequential mode)
     * \param[in] m: parallel mesh (not used)
     */
    Implementation(std::shared_ptr<Mesh<true>>) {
      reportUnsupportedParallelComputations();
    }  // end of Implementation

#endif /* MFEM_USE_MPI */

    /*!
     * \brief constructor
     * \param[in] m: sequential mesh
     */
    Implementation(std::shared_ptr<Mesh<false>> m)
        : sequential_mesh(std::move(m)) {
      if (this->sequential_mesh.get() == nullptr) {
        raise("invalid mesh");
      }
    }  // end of Implementation

    /*!
     * \brief return a mutable pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     * \return a mutable pointer to the mesh
     */
    template <bool parallel>
    std::shared_ptr<Mesh<parallel>> getMutableMeshPointer() const {
      if constexpr (parallel) {
#ifdef MFEM_USE_MPI
        if (!this->parallel_mesh.get()) {
          reportInvalidParallelMesh();
        }
        return this->parallel_mesh;
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      } else {
        if (!this->sequential_mesh.get()) {
          reportInvalidSequentialMesh();
        }
        return this->sequential_mesh;
      }
    }  // end of getMeshPointer

    /*!
     * \brief return a pointer to the mesh
     * \tparam parallel: whether to get the parallel mesh or not
     * \return a pointer to the mesh
     */
    template <bool parallel>
    std::shared_ptr<const Mesh<parallel>> getMeshPointer() const {
      if constexpr (parallel) {
#ifdef MFEM_USE_MPI
        if (!this->parallel_mesh.get()) {
          reportInvalidParallelMesh();
        }
        return this->parallel_mesh;
#else  /* MFEM_USE_MPI */
        reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
      } else {
        if (!this->sequential_mesh.get()) {
          reportInvalidSequentialMesh();
        }
        return this->sequential_mesh;
      }
    }  // end of getMeshPointer

    //   std::shared_ptr<SubMesh<true>> getParallelSubMesh(
    //       Context& ctx, const Parameter& p) const noexcept {
    //     return this->pimpl->getSubMesh<true>(ctx, p);
    //   }  // end of getParallelSubMesh
    //
    //   std::shared_ptr<SubMesh<false>> getSequentialSubMesh(
    //       Context& ctx, const Parameter& p) const noexcept {
    //     return this->pimpl->getSubMesh<false>(ctx, p);
    //   }  // end of getSequentialSubMesh

    /*!
     * \brief return if this object is built to run parallel computations
     * \return true if parallel computations are described
     */
    [[nodiscard]] bool describesAParallelComputation() const noexcept {
#ifdef MFEM_USE_MPI
      return this->parallel_mesh.get() != nullptr;
#else  /* MFEM_USE_MPI */
      return false;
#endif /* MFEM_USE_MPI */
    }  // end of describesAParallelComputation

    /*!
     * \brief set material names
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     * \return true on success, false on failure
     */
    [[nodiscard]] bool setMaterialsNames(
        Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
      return setMeshObjectNames(ctx, this->materials_names, ids,
                                getMaterialsAttributes(*this),
                                "setMaterialsNames", "material");
    }  // end of setMaterialsNames

    /*!
     * \brief set boundary names
     * \param[in, out] ctx: execution context
     * \param[in] ids: mapping between mesh identifiers and names
     * \return true on success, false on failure
     */
    [[nodiscard]] bool setBoundariesNames(
        Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
      return setMeshObjectNames(ctx, this->boundaries_names, ids,
                                getBoundariesAttributes(*this),
                                "setBoundariesNames", "boundary");
    }  // end of setBoundariesNames

    /*!
     * \brief return the material name associated with the given identifier
     * \param[in, out] ctx: execution context
     * \param[in] id: material identifier
     * \return the material name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \note the method only fails if the material identifier is not defined in
     * the mesh
     */
    [[nodiscard]] std::optional<std::string> getMaterialName(
        Context& ctx, const size_type id) const noexcept {
      const auto& mids = getMaterialsAttributes(*this);
      if (mids.Find(id) == -1) {
        return ctx.registerErrorMessage(
            "no material id '" + std::to_string(id) + "' defined in the mesh");
      }
      const auto p = this->materials_names.find(id);
      if (p == this->materials_names.end()) {
        return std::string{};
      }
      return p->second;
    }  // end of getMaterialName

    /*!
     * \brief return the boundary name associated with the given identifier
     * \param[in, out] ctx: execution context
     * \param[in] id: boundary identifier
     * \return the boundary name associated with the given identifier, if it is
     * defined. If the identifier exists but has no name, an empty string is
     * returned.
     * \note the method only fails if the boundary identifier is not defined in
     * the mesh
     */
    std::optional<std::string> getBoundaryName(
        Context& ctx, const size_type id) const noexcept {
      const auto& bids = getBoundariesAttributes(*this);
      if (bids.Find(id) == -1) {
        return ctx.registerErrorMessage(
            "no boundary id '" + std::to_string(id) + "' defined in the mesh");
      }
      const auto p = this->boundaries_names.find(id);
      if (p == this->boundaries_names.end()) {
        return std::string{};
      }
      return p->second;
    }  // end of getBoundaryName

    /*!
     * \brief return the list of materials identifiers described by the given
     * parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter
     * \return the list of materials identifiers
     * \note The parameter may hold:
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings and integers.
     * Integers are directly interpreted as materials identifiers.
     * Strings are interpreted as regular expressions which allows the selection
     * of materials by names.
     */
    [[nodiscard]] std::optional<std::vector<size_type>> getMaterialsIdentifiers(
        Context& ctx, const Parameter& p) const noexcept {
      return selectMeshObjectsIdentifiers(ctx, getMaterialsAttributes(*this),
                                          this->materials_names, p, "material",
                                          "getMaterialsIdentifiers");
    }  // end of getMaterialsIdentifiers

    /*!
     * \brief return the list of boundaries identifiers described by the given
     * parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter
     * \return the list of boundaries identifiers
     * \note The parameter may hold:
     * - an integer
     * - a string
     * - a vector of parameters which must be either strings and integers.
     * Integers are directly interpreted as boundaries identifiers.
     * Strings are interpreted as regular expressions which allows the selection
     * of boundaries by names.
     */
    [[nodiscard]] std::optional<std::vector<size_type>>
    getBoundariesIdentifiers(Context& ctx, const Parameter& p) const noexcept {
      return selectMeshObjectsIdentifiers(ctx, getBoundariesAttributes(*this),
                                          this->boundaries_names, p, "boundary",
                                          "getBoundariesIdentifiers");
    }  // end of getBoundariesIdentifiers

    /*!
     * \brief return the material identifier by the given parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter
     * \return the material identifier
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] std::optional<size_type> getMaterialIdentifier(
        Context& ctx, const Parameter& p) const noexcept {
      if (is<size_type>(p)) {
        const auto id = get<size_type>(throwing, p);
        const auto ids = getMaterialsAttributes(*this);
        if (ids.Find(id) == -1) {
          return ctx.registerErrorMessage(
              "getMaterialIdentifier: "
              "no material id for identifier '" +
              std::to_string(id) + "'");
        }
        return id;
      }
      if (!is<std::string>(p)) {
        return ctx.registerErrorMessage(
            "getMaterialIdentifier: invalid parameter type");
      }
      const auto& n = get<std::string>(throwing, p);
      for (const auto& [id, name] : this->materials_names) {
        if (name == n) {
          return id;
        }
      }
      return ctx.registerErrorMessage(
          "getMaterialIdentifier: no material named '" + n + "'");
    }  // end of getMaterialIdentifier

    /*!
     * \brief return the boundary identifier by the given parameter.
     * \param[in, out] ctx: execution context
     * \param[in] p: parameter
     * \return the boundary identifier
     * \note The parameter may hold an integer or a string.
     */
    [[nodiscard]] std::optional<size_type> getBoundaryIdentifier(
        Context& ctx, const Parameter& p) const noexcept {
      if (is<size_type>(p)) {
        const auto id = get<size_type>(throwing, p);
        const auto ids = getBoundariesAttributes(*this);
        if (ids.Find(id) == -1) {
          return ctx.registerErrorMessage(
              "getBoundaryIdentifier: "
              "no boundary id associated with identifier '" +
              std::to_string(id) + "'");
        }
        return id;
      }
      if (!is<std::string>(p)) {
        return ctx.registerErrorMessage(
            "getBoundaryIdentifier: invalid parameter type");
      }
      const auto& n = get<std::string>(throwing, p);
      for (const auto& [id, name] : this->boundaries_names) {
        if (name == n) {
          return id;
        }
      }
      return ctx.registerErrorMessage(
          "getBoundaryIdentifier: no boundary named '" + n + "'");
    }  // end of getBoundaryIdentifier

    /*!
     * \brief return the names of the materials (and their mapping with their
     * identifiers)
     * \return the mapping between material identifiers and names
     */
    [[nodiscard]] std::map<size_type, std::string> getMaterialsNames()
        const noexcept {
      return this->materials_names;
    }  // end of getMaterialsNames

    /*!
     * \brief return the names of the boundaries (and their mapping with their
     * identifiers)
     * \return the mapping between boundary identifiers and names
     */
    [[nodiscard]] std::map<size_type, std::string> getBoundariesNames()
        const noexcept {
      return this->boundaries_names;
    }  // end of getBoundariesNames

#ifdef MGIS_HAVE_TFEL

    /*!
     * \brief return the registered points in 2D
     * \param[in, out] ctx: execution context
     * \return the registered 2D points
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, Point<2>, std::less<>>>
    getPoints2D(Context& ctx) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 2) {
        return ctx.registerErrorMessage("can't return 2D points from a " +
                                        std::to_string(d) + "D mesh");
      }
      return {&(this->points2D)};
    }  // end of getPoints2D

    /*!
     * \brief return the registered points in 3D
     * \param[in, out] ctx: execution context
     * \return the registered 3D points
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, Point<3>, std::less<>>>
    getPoints3D(Context& ctx) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 3) {
        return ctx.registerErrorMessage("can't return 3D points from a " +
                                        std::to_string(d) + "D mesh");
      }
      return {&(this->points3D)};
    }  // end of getPoints3D

    /*!
     * \brief return the registered points sets in 2D
     * \param[in, out] ctx: execution context
     * \return the registered 2D points sets
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, std::vector<Point<2>>, std::less<>>>
    getPointsSets2D(Context& ctx) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 2) {
        return ctx.registerErrorMessage("can't return 2D points sets from a " +
                                        std::to_string(d) + "D mesh");
      }
      return {&(this->pointsSets2D)};
    }  // end of getPointsSets2D

    /*!
     * \brief return the registered points sets in 3D
     * \param[in, out] ctx: execution context
     * \return the registered 3D points sets
     */
    [[nodiscard]] OptionalReference<
        const std::map<std::string, std::vector<Point<3>>, std::less<>>>
    getPointsSets3D(Context& ctx) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 3) {
        return ctx.registerErrorMessage("can't return 3D points sets from a " +
                                        std::to_string(d) + "D mesh");
      }
      return {&(this->pointsSets3D)};
    }  // end of getPointsSets3D

    /*!
     * \brief add a 2D point
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \param[in] pt: coordinates of the point
     * \return true on success, false on failure
     */
    [[nodiscard]] bool addPoint(Context& ctx,
                                std::string_view n,
                                const Point<2>& pt) noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 2) {
        return ctx.registerErrorMessage("can't add a 2D point to a " +
                                        std::to_string(d) + "D mesh");
      }
      if (n.empty()) {
        return ctx.registerErrorMessage("empty name");
      }
      if (this->points2D.contains(n)) {
        return ctx.registerErrorMessage("a point named '" + std::string{n} +
                                        "' is already declared");
      }
      this->points2D.insert(std::pair<std::string, Point<2>>{n, pt});
      return true;
    }  // end of addPoint

    /*!
     * \brief add a 3D point
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \param[in] pt: coordinates of the point
     * \return true on success, false on failure
     */
    [[nodiscard]] bool addPoint(Context& ctx,
                                std::string_view n,
                                const Point<3>& pt) noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 3) {
        return ctx.registerErrorMessage("can't add a 3D point to a " +
                                        std::to_string(d) + "D mesh");
      }
      if (n.empty()) {
        return ctx.registerErrorMessage("empty name");
      }
      if (this->points3D.contains(n)) {
        return ctx.registerErrorMessage("a point named '" + std::string{n} +
                                        "' is already declared");
      }
      this->points3D.insert(std::pair<std::string, Point<3>>{n, pt});
      return true;
    }  // end of addPoint

    /*!
     * \brief return the 2D point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \return the point with the given name
     */
    [[nodiscard]] std::optional<Point<2>> getPoint2D(
        Context& ctx, std::string_view n) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 2) {
        return ctx.registerErrorMessage("can't return a 2D point from a " +
                                        std::to_string(d) + "D mesh");
      }
      const auto p = this->points2D.find(n);
      if (p == this->points2D.end()) {
        return ctx.registerErrorMessage("no point named '" + std::string{n} +
                                        "' declared");
      }
      return p->second;
    }  // end of getPoint2D

    /*!
     * \brief return the 3D point with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point
     * \return the point with the given name
     */
    [[nodiscard]] std::optional<Point<3>> getPoint3D(
        Context& ctx, std::string_view n) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 3) {
        return ctx.registerErrorMessage("can't return a 3D point from a " +
                                        std::to_string(d) + "D mesh");
      }
      const auto p = this->points3D.find(n);
      if (p == this->points3D.end()) {
        return ctx.registerErrorMessage("no point named '" + std::string{n} +
                                        "' declared");
      }
      return p->second;
    }  // end of getPoint3D

    /*!
     * \brief add a 2D points set
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point set
     * \param[in] pts: list of points
     * \return true on success, false on failure
     */
    [[nodiscard]] bool addPointsSet(Context& ctx,
                                    std::string_view n,
                                    const std::vector<Point<2>>& pts) noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 2) {
        return ctx.registerErrorMessage("can't add a 2D points set to a " +
                                        std::to_string(d) + "D mesh");
      }
      if (n.empty()) {
        return ctx.registerErrorMessage("empty name");
      }
      if (this->pointsSets2D.contains(n)) {
        return ctx.registerErrorMessage(
            "a points set named '" + std::string{n} + "' is already declared");
      }
      this->pointsSets2D.insert(
          std::pair<std::string, std::vector<Point<2>>>{n, pts});
      return true;
    }  // end of addPointsSet

    /*!
     * \brief add a 3D points set
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the point set
     * \param[in] pts: list of points
     * \return true on success, false on failure
     */
    [[nodiscard]] bool addPointsSet(Context& ctx,
                                    std::string_view n,
                                    const std::vector<Point<3>>& pts) noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 3) {
        return ctx.registerErrorMessage("can't add a 3D points set to a " +
                                        std::to_string(d) + "D mesh");
      }
      if (n.empty()) {
        return ctx.registerErrorMessage("empty name");
      }
      if (this->pointsSets3D.contains(n)) {
        return ctx.registerErrorMessage(
            "a points set named '" + std::string{n} + "' is already declared");
      }
      this->pointsSets3D.insert(
          std::pair<std::string, std::vector<Point<3>>>{n, pts});
      return true;
    }  // end of addPointsSet

    /*!
     * \brief return the 2D points set with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     * \return the points set with the given name
     */
    [[nodiscard]] OptionalReference<const std::vector<Point<2>>> getPointsSet2D(
        Context& ctx, std::string_view n) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 2) {
        return ctx.registerErrorMessage("can't return a 2D point from a " +
                                        std::to_string(d) + "D mesh");
      }
      const auto p = this->pointsSets2D.find(n);
      if (p == this->pointsSets2D.end()) {
        return ctx.registerErrorMessage("no points set named '" +
                                        std::string{n} + "' declared");
      }
      return {&(p->second)};
    }  // end of getPointsSet2D

    /*!
     * \brief return the 3D points set with the given name
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the points set
     * \return the points set with the given name
     */
    [[nodiscard]] OptionalReference<const std::vector<Point<3>>> getPointsSet3D(
        Context& ctx, std::string_view n) const noexcept {
      const auto d = getSpaceDimension(*this);
      if (d != 3) {
        return ctx.registerErrorMessage("can't return a 3D point from a " +
                                        std::to_string(d) + "D mesh");
      }
      const auto p = this->pointsSets3D.find(n);
      if (p == this->pointsSets3D.end()) {
        return ctx.registerErrorMessage("no points set named '" +
                                        std::string{n} + "' declared");
      }
      return {&(p->second)};
    }  // end of getPointsSet3D

#endif /* MGIS_HAVE_TFEL */

    ~Implementation() = default;

   private:
#ifdef MFEM_USE_MPI
    //! \brief parallel mesh
    std::shared_ptr<Mesh<true>> parallel_mesh;
#endif /* MFEM_USE_MPI */
    //! \brief sequential mesh
    std::shared_ptr<Mesh<false>> sequential_mesh;
    //! \brief mapping between materials identifiers and names
    std::map<size_type, std::string> materials_names;
    //! \brief mapping between boundaries identifiers and names
    std::map<size_type, std::string> boundaries_names;
#ifdef MGIS_HAVE_TFEL
    //! \brief points declared by the user, only valid for a 2D mesh
    std::map<std::string, Point<2>, std::less<>> points2D;
    //! \brief points declared by the user, only valid for a 3D mesh
    std::map<std::string, Point<3>, std::less<>> points3D;
    //! \brief points sets declared by the user, only valid for a 2D mesh
    std::map<std::string, std::vector<Point<2>>, std::less<>> pointsSets2D;
    //! \brief points sets declared by the user, only valid for a 3D mesh
    std::map<std::string, std::vector<Point<3>>, std::less<>> pointsSets3D;
#endif /* MGIS_HAVE_TFEL */
  };

  std::vector<std::string> MeshDiscretization::getParametersList() noexcept {
    return {MeshDiscretization::Parallel,
            MeshDiscretization::MeshFileName,
            MeshDiscretization::MeshReadMode,
            MeshDiscretization::NumberOfUniformRefinements,
            MeshDiscretization::Materials,
            MeshDiscretization::Boundaries,
            MeshDiscretization::Points,
            MeshDiscretization::PointsSets,
            MeshDiscretization::GeneralVerbosityLevel};
  }  // end of getParametersList

  MeshDiscretization::MeshDiscretization(mgis::Context& ctx,
                                         const Parameters& params)
      : pimpl(std::make_unique<Implementation>(ctx, params)) {
  }  // end of MeshDiscretization

  MeshDiscretization::MeshDiscretization(std::shared_ptr<Mesh<true>> m)
      : pimpl(std::make_unique<Implementation>(m)) {
  }  // end of MeshDiscretization

  MeshDiscretization::MeshDiscretization(std::shared_ptr<Mesh<false>> m)
      : pimpl(std::make_unique<Implementation>(m)) {
  }  // end of MeshDiscretization

  MeshDiscretization::MeshDiscretization(MeshDiscretization&&) noexcept =
      default;

  MeshDiscretization::MeshDiscretization(const MeshDiscretization&) noexcept =
      default;

  std::shared_ptr<Mesh<true>>
  MeshDiscretization::getMutableParallelMeshPointer() const noexcept {
    return this->pimpl->getMutableMeshPointer<true>();
  }  // end of getMutableSequentialMeshPointer

  std::shared_ptr<Mesh<false>>
  MeshDiscretization::getMutableSequentialMeshPointer() const noexcept {
    return this->pimpl->getMutableMeshPointer<false>();
  }  // end of getMutableSequentialMeshPointer

  std::shared_ptr<const Mesh<true>> MeshDiscretization::getParallelMeshPointer()
      const noexcept {
    return this->pimpl->getMeshPointer<true>();
  }  // end of getParallelMeshPointer

  std::shared_ptr<const Mesh<false>>
  MeshDiscretization::getSequentialMeshPointer() const noexcept {
    return this->pimpl->getMeshPointer<false>();
  }  // end of getSequentialMeshPointer

  //   std::shared_ptr<SubMesh<true>> MeshDiscretization::getParallelSubMesh(
  //       Context& ctx, const Parameter& p) const noexcept {
  //     return this->pimpl->getSubMesh<true>(ctx, p);
  //   }  // end of getParallelSubMesh
  //
  //   std::shared_ptr<SubMesh<false>> MeshDiscretization::getSequentialSubMesh(
  //       Context& ctx, const Parameter& p) const noexcept {
  //     return this->pimpl->getSubMesh<false>(ctx, p);
  //   }  // end of getSequentialSubMesh

  bool MeshDiscretization::describesAParallelComputation() const noexcept {
    return this->pimpl->describesAParallelComputation();
  }  // end of describesAParallelComputation

  bool MeshDiscretization::setMaterialsNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    return this->pimpl->setMaterialsNames(ctx, ids);
  }  // end of setMaterialsNames

  bool MeshDiscretization::setBoundariesNames(
      Context& ctx, const std::map<size_type, std::string>& ids) noexcept {
    return this->pimpl->setBoundariesNames(ctx, ids);
  }  // end of setBoundariesNames

  std::optional<std::string> MeshDiscretization::getMaterialName(
      Context& ctx, const size_type id) const noexcept {
    return this->pimpl->getMaterialName(ctx, id);
  }  // end of getMaterialName

  std::optional<std::string> MeshDiscretization::getBoundaryName(
      Context& ctx, const size_type id) const noexcept {
    return this->pimpl->getBoundaryName(ctx, id);
  }  // end of getBoundaryName

  std::optional<std::vector<size_type>>
  MeshDiscretization::getMaterialsIdentifiers(
      Context& ctx, const Parameter& p) const noexcept {
    return this->pimpl->getMaterialsIdentifiers(ctx, p);
  }  // end of getMaterialsIdentifiers

  std::optional<std::vector<size_type>>
  MeshDiscretization::getBoundariesIdentifiers(
      Context& ctx, const Parameter& p) const noexcept {
    return this->pimpl->getBoundariesIdentifiers(ctx, p);
  }  // end of getBoundariesIdentifiers

  std::optional<size_type> MeshDiscretization::getMaterialIdentifier(
      Context& ctx, const Parameter& p) const noexcept {
    return this->pimpl->getMaterialIdentifier(ctx, p);
  }  // end of getMaterialIdentifier

  std::optional<size_type> MeshDiscretization::getBoundaryIdentifier(
      Context& ctx, const Parameter& p) const noexcept {
    return this->pimpl->getBoundaryIdentifier(ctx, p);
  }  // end of getBoundaryIdentifier

  std::map<size_type, std::string> MeshDiscretization::getMaterialsNames()
      const noexcept {
    return this->pimpl->getMaterialsNames();
  }  // end of getMaterialsNames

  std::map<size_type, std::string> MeshDiscretization::getBoundariesNames()
      const noexcept {
    return this->pimpl->getBoundariesNames();
  }  // end of getBoundariesNames

#ifdef MGIS_HAVE_TFEL

  OptionalReference<const std::map<std::string, Point<2>, std::less<>>>
  MeshDiscretization::getPoints2D(Context& ctx) const noexcept {
    return this->pimpl->getPoints2D(ctx);
  }  // end of getPoints2D

  OptionalReference<const std::map<std::string, Point<3>, std::less<>>>
  MeshDiscretization::getPoints3D(Context& ctx) const noexcept {
    return this->pimpl->getPoints3D(ctx);
  }  // end of getPoints3D

  OptionalReference<
      const std::map<std::string, std::vector<Point<2>>, std::less<>>>
  MeshDiscretization::getPointsSets2D(Context& ctx) const noexcept {
    return this->pimpl->getPointsSets2D(ctx);
  }  // end of getPointsSets2D

  OptionalReference<
      const std::map<std::string, std::vector<Point<3>>, std::less<>>>
  MeshDiscretization::getPointsSets3D(Context& ctx) const noexcept {
    return this->pimpl->getPointsSets3D(ctx);
  }  // end of getPointsSets3D

  bool MeshDiscretization::addPoint(Context& ctx,
                                    std::string_view n,
                                    const Point<2>& pt) noexcept {
    return this->pimpl->addPoint(ctx, n, pt);
  }  // end of addPoint

  bool MeshDiscretization::addPoint(Context& ctx,
                                    std::string_view n,
                                    const Point<3>& pt) noexcept {
    return this->pimpl->addPoint(ctx, n, pt);
  }  // end of addPoint

  std::optional<Point<2>> MeshDiscretization::getPoint2D(
      Context& ctx, std::string_view n) const noexcept {
    return this->pimpl->getPoint2D(ctx, n);
  }  // end of getPoint2D

  std::optional<Point<3>> MeshDiscretization::getPoint3D(
      Context& ctx, std::string_view n) const noexcept {
    return this->pimpl->getPoint3D(ctx, n);
  }  // end of getPoint3D

  bool MeshDiscretization::addPointsSet(
      Context& ctx,
      std::string_view n,
      const std::vector<Point<2>>& pts) noexcept {
    return this->pimpl->addPointsSet(ctx, n, pts);
  }  // end of addPointsSet

  bool MeshDiscretization::addPointsSet(
      Context& ctx,
      std::string_view n,
      const std::vector<Point<3>>& pts) noexcept {
    return this->pimpl->addPointsSet(ctx, n, pts);
  }  // end of addPointsSet

  OptionalReference<const std::vector<Point<2>>>
  MeshDiscretization::getPointsSet2D(Context& ctx,
                                     std::string_view n) const noexcept {
    return this->pimpl->getPointsSet2D(ctx, n);
  }  // end of getPointsSet2D

  OptionalReference<const std::vector<Point<3>>>
  MeshDiscretization::getPointsSet3D(Context& ctx,
                                     std::string_view n) const noexcept {
    return this->pimpl->getPointsSet3D(ctx, n);
  }  // end of getPointsSet3D

#endif /* MGIS_HAVE_TFEL */

  MeshDiscretization::~MeshDiscretization() = default;

  const mfem::Array<size_type>& getMaterialsAttributes(
      const MeshDiscretization& m) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& mesh = m.getMesh<true>();
      return mesh.attributes;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    const auto& mesh = m.getMesh<false>();
    return mesh.attributes;
  }  // end of getMaterialsAttributes

  const mfem::Array<size_type>& getBoundariesAttributes(
      const MeshDiscretization& m) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      const auto& mesh = m.getMesh<true>();
      return mesh.bdr_attributes;
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    const auto& mesh = m.getMesh<false>();
    return mesh.bdr_attributes;
  }  // end of getBoundariesAttributes

  size_type getSpaceDimension(const MeshDiscretization& m) noexcept {
    if (m.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return m.getMesh<true>().SpaceDimension();
#else  /* MFEM_USE_MPI */
      reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return m.getMesh<false>().SpaceDimension();
  }  // end of getSpaceDimension

  template <>
  bool getInformation<MeshDiscretization>(
      Context&, std::ostream& os, const MeshDiscretization& m) noexcept {
    const auto& mnames = m.getMaterialsNames();
    os << "# Mesh\n\n"
       << "- space dimension: " << getSpaceDimension(m);
    if (!mnames.empty()) {
      os << "\n\n## Materials\n";
      for (const auto& [id, n] : mnames) {
        os << "\n- '" << n << "' associated with identifier (" << id << ")";
      }
    }
    const auto& bnames = m.getBoundariesNames();
    if (!bnames.empty()) {
      os << "\n\n## Boundaries\n";
      for (const auto& [id, n] : bnames) {
        os << "\n- '" << n << "' associated with identifier (" << id << ")";
      }
    }
    return true;
  }  // end of info

  bool operator==(const MeshDiscretization& lhs,
                  const MeshDiscretization& rhs) noexcept {
    const auto parallel = lhs.describesAParallelComputation();
    if (parallel != rhs.describesAParallelComputation()) {
      return false;
    }
    if (parallel) {
#ifdef MFEM_USE_MPI
      return (&(lhs.getMesh<true>())) == (&(rhs.getMesh<true>()));
#else
      reportUnsupportedParallelComputations();
#endif
    }
    return (&(lhs.getMesh<false>())) == (&(rhs.getMesh<false>()));
  }  // end of operator==

  bool operator!=(const MeshDiscretization& lhs,
                  const MeshDiscretization& rhs) noexcept {
    return !(lhs == rhs);
  }  // end of operator !=

  std::vector<size_type> getMaterialsIdentifiers(attributes::Throwing,
                                                 const MeshDiscretization& m,
                                                 const Parameter& p) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return m.getMaterialsIdentifiers(ctx, p) | or_raise;
  }

  std::vector<size_type> getBoundariesIdentifiers(attributes::Throwing,
                                                  const MeshDiscretization& m,
                                                  const Parameter& p) {
    auto ctx = Context{};
    auto or_raise = ctx.getThrowingFailureHandler();
    return m.getBoundariesIdentifiers(ctx, p) | or_raise;
  }  // end of getBoundariesIdentifiers

#ifdef MFEM_USE_MPI

  MPI_Comm getMPICommunicator(const MeshDiscretization& m) noexcept {
    const auto parallel = m.describesAParallelComputation();
    if (parallel) {
      return m.getMesh<true>().GetComm();
    }
    return MPI_COMM_WORLD;
  }  // end of getMPICommunicator

  bool isMainProcess(const MeshDiscretization& m) noexcept {
    int rank = 0;
    MPI_Comm_rank(getMPICommunicator(m), &rank);
    return rank == 0;
  }  // isMainProcess

#endif MFEM_USE_MPI

#ifdef MGIS_HAVE_TFEL

  template <>
  MFEM_MGIS_EXPORT std::optional<Point<2>> makePoint<2>(
      Context& ctx, const MeshDiscretization& m, const Parameter& p) noexcept {
    const auto opts = m.getPoints<2>(ctx);
    if (isInvalid(opts)) {
      return {};
    }
    return makePoint<2>(ctx, *opts, p);
  }  // end of makePoint<2>

  template <>
  MFEM_MGIS_EXPORT std::optional<Point<3>> makePoint<3>(
      Context& ctx, const MeshDiscretization& m, const Parameter& p) noexcept {
    const auto opts = m.getPoints<3>(ctx);
    if (isInvalid(opts)) {
      return {};
    }
    return makePoint<3>(ctx, *opts, p);
  }  // end of makePoint<3>

  template <>
  MFEM_MGIS_EXPORT std::optional<std::vector<Point<2>>> makePointsSet<2>(
      Context& ctx, const MeshDiscretization& m, const Parameter& p) noexcept {
    const auto opointsSets = m.getPointsSets<2>(ctx);
    if (isInvalid(opointsSets)) {
      return {};
    }
    const auto opts = m.getPoints<2>(ctx);
    if (isInvalid(opts)) {
      return {};
    }
    return makePointsSet<2>(ctx, *opointsSets, *opts, p);
  }  // end of makePointsSet<2>

  template <>
  MFEM_MGIS_EXPORT std::optional<std::vector<Point<3>>> makePointsSet<3>(
      Context& ctx, const MeshDiscretization& m, const Parameter& p) noexcept {
    const auto opointsSets = m.getPointsSets<3>(ctx);
    if (isInvalid(opointsSets)) {
      return {};
    }
    const auto opts = m.getPoints<3>(ctx);
    if (isInvalid(opts)) {
      return {};
    }
    return makePointsSet<3>(ctx, *opointsSets, *opts, p);
  }  // end of makePointsSet<3>

  template <>
  MFEM_MGIS_EXPORT std::optional<std::vector<Point<2>>> makePointsOnCurve<2>(
      Context& ctx, const MeshDiscretization& m, const Parameters& p) noexcept {
    const auto opts = m.getPoints<2>(ctx);
    if (isInvalid(opts)) {
      return {};
    }
    return makePointsOnCurve<2>(ctx, *opts, p);
  }  // end of makePointsOnCurve<2>

  template <>
  MFEM_MGIS_EXPORT std::optional<std::vector<Point<3>>> makePointsOnCurve<3>(
      Context& ctx, const MeshDiscretization& m, const Parameters& p) noexcept {
    const auto opts = m.getPoints<3>(ctx);
    if (isInvalid(opts)) {
      return {};
    }
    return makePointsOnCurve<3>(ctx, *opts, p);
  }  // end of makePointsOnCurve<3>

#endif /* MGIS_HAVE_TFEL */

}  // end of namespace mfem_mgis
