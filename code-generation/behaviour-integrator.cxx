#include <map>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <functional>
#include <string_view>
#include <filesystem>
#include <ginac/ginac.h>

using size_type = size_t;

auto icste = GiNaC::symbol("icste");
// sqrt(GiNaC::numeric(2)) / 2;

GiNaC::ex simplify(const GiNaC::ex e) { return expand(e); }

template <typename Exception = std::runtime_error>
[[noreturn]] void raise();

template <typename Exception = std::runtime_error, typename... Args>
[[noreturn]] void raise(Args &&...);

template <typename Exception>
void raise() {
  Exception e;
  throw(std::move(e));
}  // end of raise

template <typename Exception, typename... Args>
void raise(Args &&...a) {
  Exception e(std::forward<Args...>(a...));
  throw(std::move(e));
}  // end of raise

std::string makeUpperCase(const std::string &n) {
  std::string s(n);
  std::string::const_iterator p = n.begin();
  std::string::iterator p2 = s.begin();
  for (; p != n.end(); ++p, ++p2) {
    *p2 = static_cast<char>(toupper(*p));
  }
  return s;
}  // end of makeUpperCase

struct BehaviourIntegratorDescription {
  using BMatrixGenerator = std::pair<bool, GiNaC::matrix> (*)(
      std::ostream &, const std::string &, const bool);
  std::string name;
  std::string hypothesis;
  //! \brief name of the generator
  std::string generator_name;
  std::string unknown_name;
  BMatrixGenerator generator;
  /*!
   * \brief boolean stating if the values of the unknowns are passed to the
   * behaviour as an external state variable
   */
  bool requires_unknown_value_as_external_state_variable = false;
  bool isotropic = true;
};  // end of struct BehaviourIntegratorDescription

static bool isAxisymmetricalHypothesis(const std::string &h) {
  if ((h == "Axisymmetrical") ||
      (h == "AxisymmetricalGeneralisedPlaneStrain") ||
      (h == "AxisymmetricalGeneralisedPlaneStress")) {
    return true;
  }
  return false;
}  // end of isAxisymmetricalHypothesis

//! \return true if the integrator is based on a standard finite strain base
static bool isFiniteStrainIntegrator(const BehaviourIntegratorDescription &d) {
  return d.generator_name == "StandardFiniteStrainMechanics";
}  // end of isFiniteStrainIntegrator

/*!
 * \return the modelling hypothesis as in the text of the documentation
 * \param[in] h: modelling hypothesis
 */
static std::string getHypothesisDescription(const std::string &h) {
  if (h == "PlaneStrain") {
    return "plane strain";
  } else if (h == "PlaneStress") {
    return "plane stress";
  } else if (h == "Tridimensional") {
    return "tridimensional";
  }
  raise("getHypothesisDescription: unsupported hypothesis '" + h + "'");
}  // end of getHypothesisDescription

/*!
 * \return the modelling hypothesis as an enumeration value of MGIS
 * \param[in] h: modelling hypothesis
 */
static std::string getHypothesisEnumeration(const std::string &h) {
  if (h == "PlaneStrain") {
    return "PLANESTRAIN";
  } else if (h == "PlaneStress") {
    return "PLANESTRESS";
  } else if (h == "Tridimensional") {
    return "TRIDIMENSIONAL";
  }
  raise("getHypothesisEnumeration: unsupported hypothesis '" + h + "'");
}  // end of getHypothesisEnumeration

/*!
 * \brief write the file description
 * \param[out] os: output stream
 * \param[in] cn: class name
 * \param[in] a: action performed by the file, `declares` or `implements`
 */
static void writeFileDescription(std::ostream &os,
                                 const std::string &cn,
                                 const std::string &a) {
  const auto l = " * `" + cn + "` class";
  os << "/*!\n"
     << " * \\file\n"
     << " * \\brief  This file " << a << " the\n";
  if (l.size() <= 80) {
    os << l << '\n';
  } else {
    os << " * `" << cn << "`\n"
       << " * class\n";
  }
  os << " */\n\n";
}  // end of writeFileDescription

static void writeHeaderGuardAtBeginnnigOfFile(std::ostream &os,
                                              const std::string &cn,
                                              const std::string &ftype) {
  os << "#ifndef LIB_MFEM_MGIS_" << makeUpperCase(cn) << "_" << ftype << '\n'
     << "#define LIB_MFEM_MGIS_" << makeUpperCase(cn) << "_" << ftype << "\n\n";
}  // end of writeHeaderGuardAtBeginnnigOfFile

static void writeHeaderGuardAtEndOfFile(std::ostream &os,
                                        const std::string &cn,
                                        const std::string &ftype) {
  os << "#endif /* LIB_MFEM_MGIS_" << makeUpperCase(cn) << "_" << ftype
     << "*/\n";
}  // end of writeHeaderGuardAtEndOfFile

static std::vector<GiNaC::symbol> makeShapeFunctionDerivatives(
    std::ostream &os, const std::string &idx, const size_type d, const bool b) {
  std::vector<GiNaC::symbol> dNs;
  for (size_type i = 0; i != d; ++i) {
    const auto dN = "dN" + idx + "_" + std::to_string(i);
    if (b) {
      os << "const auto " << dN << " = dN(n" << idx << ", " << i << ");\n";
    }
    dNs.emplace_back(dN);
  }
  return dNs;
}  // end of makeShapeFunctionDerivative

static GiNaC::matrix makeVectorOfSymbols(const std::string &name,
                                         const size_type s) {
  GiNaC::matrix v(s, 1);
  for (size_type i = 0; i != s; ++i) {
    v(i, 0) = GiNaC::symbol(name + '[' + std::to_string(i) + "]");
  }
  return v;
}

static GiNaC::matrix makeMatrixOfSymbols(const std::string &name,
                                         const size_type n,
                                         const size_type m,
                                         const size_type o = 0) {
  GiNaC::matrix v(n, m);
  for (size_type i = 0; i != n; ++i) {
    for (size_type j = 0; j != m; ++j) {
      const auto offset = std::to_string(i * m + j + o);
      v(i, j) = GiNaC::symbol(name + '[' + offset + ']');
    }
  }
  return v;
}

std::string getMatrixComponentId(const std::string &n,
                                 const size_type i,
                                 const size_type j) {
  return n + "_" + std::to_string(i) + "_" + std::to_string(j);
}

std::string generateUnknownOffset(const std::string &nid, const size_type idx) {
  if (idx == 0) {
    return nid;
  } else if (idx == 1) {
    return nid + " + nnodes";
  }
  return nid + " + " + std::to_string(idx) + " * nnodes";
}

void generateUnknownOffsets(std::ostream &os,
                            const std::string &nid,
                            const size_type d) {
  if (d == 1) {
    return;
  }
  for (size_type i = 0; i != d; ++i) {
    os << "const auto " << nid << "_" << i << " = "
       << generateUnknownOffset(nid, i) << ";\n";
  }
}

/*!
 * \brief return the offset associated with an unknown component
 * \param[in] nid: node id
 * \param[in] c: component
 * \param[in] s: size of the unknowns
 * \param[in] b: if true, assume that the offsets of the unknowns have already
 * been computed.
 * \return the offset
 */
std::string getUnknownOffset(const std::string &nid,
                             const size_type c,
                             const size_type s,
                             const bool b = true) {
  if (!b) {
    return generateUnknownOffset(nid, c);
  } else {
    if (s == 1) {
      return nid;
    }
    return nid + "_" + std::to_string(c);
  }
}

GiNaC::matrix makeVectorOfUnknowns(std::ostream &os,
                                   const std::string &name,
                                   const std::string &nid,
                                   const size_type s,
                                   const bool b) {
  GiNaC::matrix v(s, 1);
  for (size_type i = 0; i != s; ++i) {
    const auto u = name + "_" + std::to_string(i);
    os << "const auto " << u << " = " << name << "["
       << getUnknownOffset(nid, i, s, b) << "];\n";
    v(i, 0) = GiNaC::symbol(u);
  }
  return v;
}

void generateUpdateGradient(std::ostream &os,
                            const BehaviourIntegratorDescription &d) {
  os << "inline void\n"
     << d.name << "::updateGradients(std::span<real> &g,\n"
     << "const mfem::Vector &u,\n"
     << "const mfem::DenseMatrix &dN,\n"
     << "const size_type ni) noexcept {\n";
  const auto [b, Bi] = d.generator(os, "i", true);
  if (Bi.cols() != 1) {
    os << "const auto nnodes = dN.NumRows();\n";
  }
  const auto u = makeVectorOfUnknowns(os, "u", "ni", Bi.cols(), false);
  const auto g = Bi.mul(u);
  for (size_type i = 0; i != Bi.rows(); ++i) {
    os << "g[" << i << "] += " << simplify(g(i, 0)) << ";\n";
  }
  os << "} // end of updateGradients\n\n";
}

void generateUpdateInnerForces(std::ostream &os,
                               const BehaviourIntegratorDescription &d) {
  os << "inline void\n"
     << d.name << "::updateInnerForces(mfem::Vector &Fe,\n"
     << "const std::span<const real> &s,\n"
     << "const mfem::DenseMatrix &dN,\n"
     << "const real w,\n"
     << "const size_type ni) const noexcept {\n";
  const auto [bi, Bi] = d.generator(os, "i", true);
  const auto S = makeVectorOfSymbols("s", Bi.rows());
  const auto Fe = transpose(Bi).mul(S);
  if (Bi.cols() != 1) {
    os << "const auto nnodes = dN.NumRows();\n";
  }
  generateUnknownOffsets(os, "ni", Bi.cols());
  for (size_type i = 0; i != Bi.cols(); ++i) {
    const auto ni = getUnknownOffset("ni", i, Bi.cols());
    os << "Fe[" << ni << "] += "
       << "w * (" << simplify(Fe(i, 0)) << ");\n";
  }
  os << "} // end of updateInnerForces\n\n";
}  // end of generateUpdateInnerForces

void generateUpdateStiffnessMatrix(std::ostream &os,
                                   const BehaviourIntegratorDescription &d) {
  if (d.requires_unknown_value_as_external_state_variable) {
    os << "inline void\n"
       << d.name << "::updateStiffnessMatrix(mfem::DenseMatrix &Ke,\n"
       << "const std::span<const real> &Kip,\n"
       << "const mfem::Vector &N,\n"
       << "const mfem::DenseMatrix &dN,\n"
       << "const real w,\n"
       << "const size_type ni) const noexcept {\n";
  } else {
    os << "inline void\n"
       << d.name << "::updateStiffnessMatrix(mfem::DenseMatrix &Ke,\n"
       << "const std::span<const real> &Kip,\n"
       << "const mfem::DenseMatrix &dN,\n"
       << "const real w,\n"
       << "const size_type ni) const noexcept {\n";
  }
  os << "const auto nnodes = dN.NumRows();\n";
  const auto [bi, Bi] = d.generator(os, "i", true);
  const auto K = makeMatrixOfSymbols("Kip", Bi.rows(), Bi.rows());
  generateUnknownOffsets(os, "ni", Bi.cols());
  os << "for (size_type nj = 0; nj != nnodes; ++nj) {\n";
  const auto [bj, Bj] = d.generator(os, "j", true);
  const auto Ke = transpose(Bi).mul(K).mul(Bj);
  generateUnknownOffsets(os, "nj", Bj.cols());
  if (d.requires_unknown_value_as_external_state_variable) {
    if (Bi.cols() != 1) {
      raise("generateUpdateStiffnessMatrix: invalid call");
    }
    const auto ds_ddu =
        makeMatrixOfSymbols("Kip", Bi.rows(), 1, Bi.rows() * Bi.rows());
    const auto dFe_ddu = transpose(Bi).mul(ds_ddu);
    for (size_type i = 0; i != Bi.cols(); ++i) {
      const auto ni = getUnknownOffset("ni", i, Bi.cols());
      for (size_type j = 0; j != Bi.cols(); ++j) {
        const auto nj = getUnknownOffset("nj", j, Bj.cols());
        const auto Ke1v = simplify(Ke(i, j));
        const auto Ke2v = simplify(dFe_ddu(i, 0));
        os << "Ke(" << ni << ", " << nj << ") += "
           << "w * (" << simplify(Ke1v) << " + (" << Ke2v << ") * "
           << "N[" << nj << "]);\n";
      }
    }
  } else {
    for (size_type i = 0; i != Bi.cols(); ++i) {
      const auto ni = getUnknownOffset("ni", i, Bi.cols());
      for (size_type j = 0; j != Bi.cols(); ++j) {
        const auto nj = getUnknownOffset("nj", j, Bj.cols());
        os << "Ke(" << ni << ", " << nj << ") += "
           << "w * (" << simplify(Ke(i, j)) << ");\n";
      }
    }
  }
  os << "} // end of for (size_type nj = 0; nj != nnodes; ++nj)\n"
     << "} // end of updateStiffnessMatrix\n\n";
}  // end of generateUpdateStiffnessMatrix

std::pair<bool, GiNaC::matrix> makePlaneStrainSmallStrainMechanicsBMatrix(
    std::ostream &os, const std::string &nid, const bool b) {
  auto B = GiNaC::matrix(4, 2);
  const auto Bn = "B" + nid;
  const auto dN = [&nid](const size_type i) {
    return "dN(n" + nid + ", " + std::to_string(i) + ")";
  };
  const auto generate = [&os, &B, &Bn, &b](const size_type i, const size_type j,
                                           const std::string &v) {
    B(i, j) = GiNaC::symbol(getMatrixComponentId(Bn, i, j));
    if (b) {
      os << "const auto " << B(i, j) << " = " << v << ";\n";
    }
  };
  generate(0, 0, dN(0));
  generate(1, 1, dN(1));
  generate(3, 0, dN(1) + " * icste");
  generate(3, 1, dN(0) + " * icste");
  return std::make_pair(false, B);
}

std::pair<bool, GiNaC::matrix> makePlaneStressSmallStrainMechanicsBMatrix(
    std::ostream &os, const std::string &nid, const bool b) {
  return makePlaneStrainSmallStrainMechanicsBMatrix(os, nid, b);
}

std::pair<bool, GiNaC::matrix> makeTridimensionalSmallStrainMechanicsBMatrix(
    std::ostream &os, const std::string &nid, const bool b) {
  auto B = GiNaC::matrix(6, 3);
  const auto Bn = "B" + nid;
  const auto dN = [&nid](const size_type i) {
    return "dN(n" + nid + ", " + std::to_string(i) + ")";
  };
  const auto generate = [&os, &B, &Bn, &b](const size_type i, const size_type j,
                                           const std::string &v) {
    B(i, j) = GiNaC::symbol(getMatrixComponentId(Bn, i, j));
    if (b) {
      os << "const auto " << B(i, j) << " = " << v << ";\n";
    }
  };
  generate(0, 0, dN(0));
  generate(1, 1, dN(1));
  generate(2, 2, dN(2));
  generate(3, 0, dN(1) + "* icste");
  generate(3, 1, dN(0) + "* icste");
  generate(4, 0, dN(2) + "* icste");
  generate(4, 2, dN(0) + "* icste");
  generate(5, 1, dN(2) + "* icste");
  generate(5, 2, dN(1) + "* icste");
  return std::make_pair(false, B);
}

std::pair<bool, GiNaC::matrix> makePlaneStrainFiniteStrainMechanicsBMatrix(
    std::ostream &os, const std::string &idx, const bool b) {
  const auto d = size_type{2};
  const auto dN = makeShapeFunctionDerivatives(os, idx, d, b);
  auto B = GiNaC::matrix(5, d);
  B(0, 0) = dN[0];
  B(1, 1) = dN[1];
  B(3, 0) = dN[1];  // dux_dy
  B(4, 1) = dN[0];  // duy_dx
  return std::make_pair(false, B);
}

std::pair<bool, GiNaC::matrix> makePlaneStressFiniteStrainMechanicsBMatrix(
    std::ostream &os, const std::string &idx, const bool b) {
  return makePlaneStrainFiniteStrainMechanicsBMatrix(os, idx, b);
}

std::pair<bool, GiNaC::matrix> makeTridimensionalFiniteStrainMechanicsBMatrix(
    std::ostream &os, const std::string &idx, const bool b) {
  const auto d = size_type{3};
  const auto dN = makeShapeFunctionDerivatives(os, idx, d, b);
  auto B = GiNaC::matrix(9, d);
  B(0, 0) = dN[0];
  B(1, 1) = dN[1];
  B(2, 2) = dN[2];
  B(3, 0) = dN[1];  // dux_dy
  B(4, 1) = dN[0];  // duy_dx
  B(5, 0) = dN[2];  // dux_dz
  B(6, 2) = dN[0];  // duz_dx
  B(7, 1) = dN[2];  // duy_dz
  B(8, 2) = dN[1];  // duz_dy
  return std::make_pair(false, B);
}

template <unsigned short dime>
std::pair<bool, GiNaC::matrix> makeHeatTransferBMatrix(std::ostream &os,
                                                       const std::string &idx,
                                                       const bool b) {
  const auto dN = makeShapeFunctionDerivatives(os, idx, dime, b);
  auto B = GiNaC::matrix(dime, 1);
  for (size_type i = 0; i != dime; ++i) {
    B(i, 0) = dN[i];
  }
  return std::make_pair(false, B);
}

void generateHeaderFile(std::ostream &os,
                        const BehaviourIntegratorDescription &d) {
  const auto fs = isFiniteStrainIntegrator(d);
  const auto fs_base =
      d.hypothesis + "StandardFiniteStrainMechanicsBehaviourIntegratorBase";
  const auto [b, B] = d.generator(os, "", false);
  writeFileDescription(os, d.name, "declares");
  writeHeaderGuardAtBeginnnigOfFile(os, d.name, "HXX");
  os << "#include <array>\n"
     << "#include <mfem/linalg/densemat.hpp>\n"
     << "#include \"MFEMMGIS/Config.hxx\"\n"
     << "#include \"MFEMMGIS/BehaviourIntegratorTraits.hxx\"\n"
     << "#include \"MFEMMGIS/StandardBehaviourIntegratorCRTPBase.hxx\"\n";
  if (fs) {
    os << "#include \"MFEMMGIS/" << fs_base << ".hxx\"\n";
  }
  os << '\n'
     << "namespace mfem_mgis {\n"
     << '\n'
     << "// forward declaration\n"
     << "struct FiniteElementDiscretization;\n"
     << '\n'
     << "// forward declaration\n"
     << "struct " << d.name << ";\n"
     << '\n'
     << "/*!\n"
     << " * \\brief specialisation of the `BehaviourIntegratorTraits` class "
     << "for the\n"
     << " * `" << d.name << "`\n"
     << " * behaviour integrator\n"
     << " */\n"
     << "template<>\n"
     << "struct BehaviourIntegratorTraits<" << d.name << ">{\n"
     << "//! \\brief number of components of the unknowns\n"
     << "static constexpr size_type unknownsSize = " << B.cols() << ";\n"
     << "//! \\brief if the computation of the gradients requires the shape "
     << "functions\n"
     << "static constexpr bool gradientsComputationRequiresShapeFunctions"
     << " = false;\n"
     << "/*!\n"
     << " * \\brief if the computation of the gradients requires the "
     << "derivatives of\n"
     << " * the shape functions\n"
     << " */\n"
     << "static constexpr bool "
     << "gradientsComputationRequiresShapeFunctionsDerivatives = true;\n"
     << "//! \\brief if the external state variables are updated from the "
     << "unknowns\n";
  if (d.requires_unknown_value_as_external_state_variable) {
    os << "static constexpr bool "
          "updateExternalStateVariablesFromUnknownsValues = true;\n";
  } else {
    os << "static constexpr bool "
          "updateExternalStateVariablesFromUnknownsValues = false;\n";
  }
  os << "}; // end of struct BehaviourIntegratorTraits<>\n"
     << '\n'
     << "/*!\n"
     << " * \\brief behaviour integrator for "
     << (d.isotropic ? "isotropic" : "orthotropic");
  if (d.generator_name == "StandardSmallStrainMechanics") {
    os << " small strain mechanical\n"
       << " * behaviours";
  } else if (d.generator_name == "StandardFiniteStrainMechanics") {
    os << " finite strain mechanical\n"
       << " * behaviours";
  } else if (d.generator_name == "StationaryNonLinearHeatTransfer") {
    os << " stationary non linear heat\n"
       << " * transfer behaviours";
  } else {
    raise("generateHeaderFile: unsupported generator '" + d.generator_name +
          "'");
  }
  os << " under the " << getHypothesisDescription(d.hypothesis)
     << " hypothesis\n"
     << " */\n"
     << "struct MFEM_MGIS_EXPORT " << d.name;
  // the isotropic finite strain integrators are the base classes of the
  // Faltus2026 regularized integrators
  if (!(fs && d.isotropic)) {
    os << " final";
  }
  os << "\n: StandardBehaviourIntegratorCRTPBase<" << d.name << ">";
  if (fs) {
    os << ", " << fs_base;
  }
  os << " {\n";
  if (d.generator_name == "StandardSmallStrainMechanics") {
    os << "/*!\n"
       << " * \\brief a constant value used for the computation of\n"
       << " * symmetric tensors\n"
       << " */\n"
       << "static constexpr const auto icste = real{0.70710678118654752440};\n";
  }
  if (d.isotropic) {
    os << "//! \\brief a dummy structure\n"
       << "struct RotationMatrix {};\n";
  } else {
    os << "//! \\brief a simple alias\n"
       << "using RotationMatrix = std::array<real, 9u>;\n";
  }
  os << "/*!\n"
     << " * \\brief constructor\n"
     << " * \\param[in] fed: finite element discretization.\n"
     << " * \\param[in] m: material attribute.\n"
     << " * \\param[in] b_ptr: behaviour\n"
     << " */\n"
     << d.name << "(const FiniteElementDiscretization &fed,\n"
     << "const size_type m,\n"
     << "std::unique_ptr<const Behaviour> b_ptr);\n";
  if (d.requires_unknown_value_as_external_state_variable) {
    os << "/*!\n"
       << " * \\brief method called at the beginning of each resolution\n"
       << " *\n"
       << " * In addition to the base class setup, store a pointer to the "
       << "values of\n"
       << " * the `" << d.unknown_name
       << "` external state variable, which must be defined,\n"
       << " * mutable and not uniform.\n"
       << " *\n"
       << " * \\param[in, out] ctx: execution context\n"
       << " * \\param[in] t: time at the beginning of the time step\n"
       << " * \\param[in] dt: time increment\n"
       << " * \\return true on success\n"
       << " */\n"
       << "[[nodiscard]] bool setup(Context &ctx, const real t, "
       << "const real dt) noexcept override;\n"
       << "/*!\n"
       << " * \\brief update the external state variable\n"
       << " * corresponding to the unknowns.\n"
       << " *\n"
       << " * \\param[in] u: unknowns\n"
       << " * \\param[in] N: values of the shape functions\n"
       << " * \\param[in] o: offset associated with the\n"
       << " * current integration point\n"
       << " */\n"
       << "void updateExternalStateVariablesFromUnknownsValues(\n"
       << "const mfem::Vector &u, const mfem::Vector &N, const size_type o);\n";
  }
  os << "/*!\n"
     << " * \\brief get the rotation matrix at the given integration point\n"
     << " * \\return the rotation matrix associated with the given "
     << "integration\n"
     << " * point\n"
     << " * \\param[in] i: offset of the integration point\n"
     << " */\n"
     << "inline RotationMatrix getRotationMatrix(const size_type i) const;\n"
     << '\n'
     << "/*!\n"
     << " * \\brief rotate the gradients in the material frame\n"
     << " * \\param[in, out] g: gradients\n"
     << " * \\param[in] r: rotation matrix\n";
  if (d.isotropic) {
    os << " * \\note does nothing\n";
  }
  os << " */\n"
     << "inline void rotateGradients(std::span<real> g, "
     << "const RotationMatrix &r);\n"
     << '\n'
     << "/*!\n"
     << " * \\brief rotate the thermodynamic forces in the global frame\n";
  if (d.isotropic) {
    os << " * \\param[in] s: thermodynamic forces\n";
  } else {
    os << " * \\param[in] s: thermodynamic forces in the material frame\n";
  }
  os << " * \\param[in] r: rotation matrix\n"
     << " * \\return the thermodynamic forces in the global frame\n";
  if (d.isotropic) {
    os << " * \\note returns s unchanged\n"
       << " */\n"
       << "inline std::span<const real>\n";
  } else {
    os << " */\n"
       << "inline std::array<real, " << B.rows() << ">\n";
  }
  os << "rotateThermodynamicForces(std::span<const real> s, "
     << "const RotationMatrix &r);\n"
     << '\n'
     << "/*!\n"
     << " * \\brief rotate the tangent operator blocks in the global frame\n"
     << " * \\param[in, out] Kip: tangent operator blocks\n"
     << " * \\param[in] r: rotation matrix\n";
  if (d.isotropic) {
    os << " * \\note does nothing\n";
  }
  os << " */\n"
     << "inline void rotateTangentOperatorBlocks(std::span<real> Kip,\n"
     << "const RotationMatrix &r);\n"
     << '\n'
     << "const mfem::IntegrationRule &getIntegrationRule(\n"
     << "    const mfem::FiniteElement &e,\n"
     << "    const mfem::ElementTransformation &tr) const override;\n"
     << '\n'
     << "real getIntegrationPointWeight(mfem::ElementTransformation &tr,\n"
     << "                               const mfem::IntegrationPoint &ip)\n"
     << "                              const noexcept override;\n"
     << '\n'
     << "bool integrate(const mfem::FiniteElement &e,\n"
     << "               mfem::ElementTransformation &tr,\n"
     << "               const mfem::Vector &u,\n"
     << "               const IntegrationType it) override;\n"
     << '\n'
     << "void updateResidual(mfem::Vector &Fe,\n"
     << "                    const mfem::FiniteElement &e,\n"
     << "                    mfem::ElementTransformation &tr,\n"
     << "                    const mfem::Vector &u) override;\n"
     << '\n'
     << "void updateJacobian(mfem::DenseMatrix &Ke,\n"
     << "                    const mfem::FiniteElement &e,\n"
     << "                    mfem::ElementTransformation &tr,\n"
     << "                    const mfem::Vector &u) override;\n"
     << '\n'
     << "void computeInnerForces(mfem::Vector &Fe,\n"
     << "                        const mfem::FiniteElement &e,\n"
     << "                        mfem::ElementTransformation &tr) override;\n"
     << '\n'
     << "//! \\brief destructor\n"
     << "~" << d.name << "() override;\n"
     << '\n'
     << "protected:\n"
     << "//! \\brief allow the CRTP base class to access the protected "
     << "members\n"
     << "friend struct StandardBehaviourIntegratorCRTPBase<" << d.name << ">;\n"
     << "/*!\n"
     << " * \\brief select the integration rule for the given element and "
     << "element\n"
     << " * transformation\n"
     << " * \\param[in] e: element\n"
     << " * \\param[in] t: element transformation\n"
     << " * \\return the integration rule\n"
     << " */\n"
     << "static const mfem::IntegrationRule &selectIntegrationRule(\n"
     << "    const mfem::FiniteElement &e,\n"
     << "    const mfem::ElementTransformation &t);\n"
     << "/*!\n"
     << " * \\brief build the quadrature space for the given material\n"
     << " * \\param[in] fed: finite element discretization.\n"
     << " * \\param[in] m: material attribute.\n"
     << " * \\return the partial quadrature space\n"
     << " */\n"
     << "static std::shared_ptr<const PartialQuadratureSpace> "
     << "buildQuadratureSpace(const FiniteElementDiscretization &fed,\n"
     << "                     const size_type m);\n";
  if (!fs) {
    const auto heat = d.generator_name == "StationaryNonLinearHeatTransfer";
    const auto g = std::string(heat ? "temperature gradient" : "strain");
    const auto s = std::string(heat ? "heat flux" : "stress");
    os << "/*!\n"
       << " * \\brief update the " << g << " with the contribution of the\n"
       << " * given node\n"
       << " * \\param[in, out] g: " << g << "\n"
       << " * \\param[in] u: nodal "
       << (heat ? "temperatures" : "displacements") << "\n"
       << " * \\param[in] dN: derivatives of the shape functions\n"
       << " * \\param[in] ni: node index\n"
       << " */\n"
       << "void updateGradients(std::span<real> &g,\n"
       << "                     const mfem::Vector &u,\n"
       << "                     const mfem::DenseMatrix &dN,\n"
       << "                     const size_type ni) noexcept;\n"
       << "/*!\n"
       << " * \\brief update the inner forces of the given node with\n"
       << " * the contribution of the " << s << " of an integration point.\n"
       << " *\n"
       << " * \\param[in, out] Fe: inner forces\n"
       << " * \\param[in] s: " << s << "\n"
       << " * \\param[in] dN: derivatives of the shape functions\n"
       << " * \\param[in] w: weight of the integration point\n"
       << " * \\param[in] ni: node index\n"
       << " */\n"
       << "void updateInnerForces(mfem::Vector &Fe,\n"
       << "                       const std::span<const real> &s,\n"
       << "                       const mfem::DenseMatrix &dN,\n"
       << "                       const real w,\n"
       << "                       const size_type ni) const noexcept;\n"
       << "/*!\n"
       << " * \\brief update the stiffness matrix of the given node\n"
       << " * with the contribution of the consistent tangent operator of an\n"
       << " * integration point.\n"
       << " *\n"
       << " * \\param[in, out] Ke: stiffness matrix\n"
       << " * \\param[in] Kip: consistent tangent operator of the integration "
       << "point\n";
    if (d.requires_unknown_value_as_external_state_variable) {
      os << " * \\param[in] N: shape function values\n";
    }
    os << " * \\param[in] dN: derivatives of the shape functions\n"
       << " * \\param[in] w: weight of the integration point\n"
       << " * \\param[in] ni: node index\n"
       << " */\n"
       << "void updateStiffnessMatrix(mfem::DenseMatrix &Ke,\n"
       << "                           const std::span<const real> &Kip,\n";
    if (d.requires_unknown_value_as_external_state_variable) {
      os << "                           const mfem::Vector &N,\n";
    }
    os << "                           const mfem::DenseMatrix &dN,\n"
       << "                           const real w,\n"
       << "                           const size_type ni) const noexcept;\n";
  }
  os << '\n';
  if (d.requires_unknown_value_as_external_state_variable) {
    if (B.cols() == 1) {
      os << "/*!\n"
         << " * \\brief pointer to the external state variable\n"
         << " * associated with the unknown\n"
         << " */\n"
         << "real *uesv = nullptr;\n";
    } else {
      raise("generateHeaderFile: unsupported case");
    }
  }
  os << "};  // end of struct " << d.name << '\n'
     << '\n'
     << "}  // end of namespace mfem_mgis\n"
     << '\n';
  writeHeaderGuardAtEndOfFile(os, d.name, "HXX");
}  // end of generateHeaderFile

void generateSourceFile(std::ostream &os,
                        const BehaviourIntegratorDescription &d) {
  const auto [b, B] = d.generator(os, "", false);
  writeFileDescription(os, d.name, "implements");
  os << "#include <algorithm>\n"
     << "#include \"MGIS/Behaviour/Behaviour.hxx\"\n"
     << "#include \"MFEMMGIS/" << d.name << ".hxx\"\n\n"
     << "namespace mfem_mgis {\n\n";
  os << "const mfem::IntegrationRule &\n"
     << "" << d.name << "::selectIntegrationRule(\n"
     << "const mfem::FiniteElement &e,"
     << "const mfem::ElementTransformation &t) {\n"
     << "  const auto order = 2 * t.OrderGrad(&e);\n";
  if (isAxisymmetricalHypothesis(d.hypothesis)) {
    os << "  return mfem::IntRules.Get(e.GetGeomType(), order + 1);\n";
  } else {
    os << "  return mfem::IntRules.Get(e.GetGeomType(), order);\n";
  }
  os << "}\n"
     << '\n'
     << "std::shared_ptr<const PartialQuadratureSpace>\n"
     << "" << d.name << "::buildQuadratureSpace(\n"
     << "    const FiniteElementDiscretization &fed, const size_type m) {\n"
     << "  auto selector = [](const mfem::FiniteElement &e,\n"
     << "                     const mfem::ElementTransformation &tr)\n"
     << "      -> const mfem::IntegrationRule & {\n"
     << "    return selectIntegrationRule(e, tr);\n"
     << "  };  // end of selector\n"
     << "  return std::make_shared<PartialQuadratureSpace>(fed, m, "
        "selector);\n"
     << "}  // end of buildQuadratureSpace\n"
     << '\n'
     << d.name << "::" << d.name << "(\n"
     << "        const FiniteElementDiscretization &fed,\n"
     << "        const size_type m,\n"
     << "        std::unique_ptr<const Behaviour> b_ptr)\n"
     << "    : StandardBehaviourIntegratorCRTPBase<" << d.name << ">(\n"
     << "          buildQuadratureSpace(fed, m), std::move(b_ptr)) {\n";
  if (isFiniteStrainIntegrator(d)) {
    os << "this->checkIfAFiniteStrainBehaviourIsDeclared(throwing);\n";
  }
  if (d.isotropic) {
    os << "this->checkBehaviourSymmetry(throwing, Behaviour::ISOTROPIC);\n";
  } else {
    os << "this->checkBehaviourSymmetry(throwing, Behaviour::ORTHOTROPIC);\n";
  }
  os << "this->checkHypothesis(throwing, Hypothesis::"
     << getHypothesisEnumeration(d.hypothesis) << ");\n";
  if (d.generator_name == "StandardSmallStrainMechanics") {
    // sizes of the strain, of the stress and of the tangent operator block
    const auto n = B.rows();
    os << "this->checkBehaviourVariablesSizes(throwing, {" << n << "}, {" << n
       << "}, {{" << n << ", " << n << "}});\n";
  }
  os << "}  // end of " << d.name << '\n'
     << '\n'
     << "real " << d.name << "::getIntegrationPointWeight"
     << "(mfem::ElementTransformation &tr,\n"
     << " const mfem::IntegrationPoint &ip) const noexcept{\n";
  if (isAxisymmetricalHypothesis(d.hypothesis)) {
    os << "constexpr const real two_pi = 2 * 3.14159265358979323846;\n"
       << "return two_pi * ip.x * ip.weight * tr.Weight();\n";
  } else {
    os << "return ip.weight * tr.Weight();\n";
  }
  os << "}\n"
     << "const mfem::IntegrationRule &\n"
     << d.name << "::getIntegrationRule(\n"
     << "const mfem::FiniteElement &e, "
     << "const mfem::ElementTransformation &t)const  {\n"
     << "return " << d.name << "::selectIntegrationRule(e, t);\n"
     << "}\n";
  if (d.requires_unknown_value_as_external_state_variable) {
    const auto &n = d.unknown_name;
    os << "\n"
       << "bool " << d.name << "::setup(Context &ctx, const real t, "
       << "const real dt) noexcept {\n"
       << "using MutableFieldHolder =\n"
       << "mgis::behaviour::MaterialStateManager::MutableFieldHolder;\n"
       << "if (!BehaviourIntegratorBase::setup(ctx, t, dt)) {\n"
       << "return false;\n"
       << "}\n";
    if (B.cols() == 1) {
      os << "const auto &pev = this->s1.external_state_variables.find("
         << "\"" << n << "\");\n"
         << "if (pev == this->s1.external_state_variables.end()) {\n"
         << "return ctx.registerErrorMessage(\n"
         << "\"external state variable '" << n << "' is not defined\");\n"
         << "}\n"
         << "if (std::holds_alternative<std::monostate>(pev->second)) {\n"
         << "return ctx.registerErrorMessage(\n"
         << "\"external state variable '" << n << "' is not initialized\");\n"
         << "} else if (std::holds_alternative<std::span<const real>>"
         << "(pev->second)) {\n"
         << "return ctx.registerErrorMessage(\n"
         << "\"external state variable '" << n
         << "' shall not be immutable\");\n"
         << "} else {\n"
         << "if (!std::holds_alternative<MutableFieldHolder>(pev->second)) {\n"
         << "return ctx.registerErrorMessage(\n"
         << "\"unsupported field holder for external state variable \"\n"
         << "\"'" << n << "'\");\n"
         << "}\n"
         << "auto &value = std::get<MutableFieldHolder>(pev->second).value;\n"
         << "if (std::holds_alternative<std::span<real>>(value)) {\n"
         << "this->uesv = std::get<std::span<real>>(value).data();\n"
         << "} else if (std::holds_alternative<std::vector<real>>(value)) {\n"
         << "this->uesv = std::get<std::vector<real>>(value).data();\n"
         << "} else {\n"
         << "return ctx.registerErrorMessage(\n"
         << "\"external state variable '" << n << "' shall not be uniform\");\n"
         << "}\n"
         << "}\n"
         << "return true;\n"
         << "}  // end of setup\n";
    } else {
      raise("generateSourceFile: unsupported case");
    }
    os << "\n"
       << "void " << d.name
       << "::updateExternalStateVariablesFromUnknownsValues(\n"
       << "const mfem::Vector &u, const mfem::Vector &N, "
       << "const size_type o){\n";
    if (B.cols() == 1) {
      os << "auto &ev = this->uesv[o];\n"
         << "ev = 0;\n"
         << "for(size_type i=0; i != u.Size(); ++i){\n"
         << "  ev += u[i] * N[i];"
         << "}\n";
    } else {
      raise(
          "unsupport case: multicomponents unknowns as external state "
          "variables is not supported");
    }
    os << "}\n\n";
  }
  if (d.isotropic) {
    os << d.name << "::RotationMatrix\n"
       << d.name << "::getRotationMatrix(const size_type) const{\n"
       << "return RotationMatrix{};\n"
       << "} // end of getRotationMatrix\n"
       << '\n';
  } else {
    os << d.name << "::RotationMatrix\n"
       << d.name << "::getRotationMatrix(const size_type i) const{\n"
       << "return this->get_rotation_fct_ptr(this->r2D, this->r3D, i);\n"
       << "} // end of getRotationMatrix\n"
       << '\n';
  }
  if (d.isotropic) {
    os << "void " << d.name
       << "::rotateGradients(std::span<real>, const RotationMatrix &){\n"
       << "} // end of rotateGradients\n";
  } else {
    os << "void " << d.name
       << "::rotateGradients(std::span<real> g, const RotationMatrix &r){\n"
       << "this->b.rotate_gradients_ptr(g.data(), g.data(), r.data());\n"
       << "} // end of rotateGradients\n";
  }
  os << '\n';
  if (d.isotropic) {
    os << "std::span<const real>\n"
       << d.name << "::rotateThermodynamicForces(std::span<const real> s, "
       << "const RotationMatrix &){\n"
       << "return s;\n"
       << "}\n\n";
  } else {
    os << "std::array<real, " << B.rows() << ">\n"
       << d.name << "::rotateThermodynamicForces(std::span<const real> s, "
       << "const RotationMatrix &r){\n"
       << "std::array<real, " << B.rows() << "> rs;\n"
       << "std::copy(s.begin(), s.end(), rs.begin());\n"
       << "this->b.rotate_thermodynamic_forces_ptr(rs.data(), rs.data(), "
          "r.data());\n"
       << "return rs;\n"
       << "}\n\n";
  }
  if (d.isotropic) {
    os << "void " << d.name
       << "::rotateTangentOperatorBlocks(std::span<real>,\n"
       << "const RotationMatrix &){\n"
       << "}\n";
  } else {
    os << "void " << d.name
       << "::rotateTangentOperatorBlocks(std::span<real> Kip,\n"
       << "const RotationMatrix &r){\n"
       << "this->b.rotate_tangent_operator_blocks_ptr(Kip.data(), Kip.data(), "
          "r.data());\n"
       << "}\n";
  }
  os << '\n';
  if (!isFiniteStrainIntegrator(d)) {
    generateUpdateGradient(os, d);
    generateUpdateInnerForces(os, d);
    generateUpdateStiffnessMatrix(os, d);
  }
  os << "bool " << d.name << "::integrate(const mfem::FiniteElement &e,\n"
     << "                                 mfem::ElementTransformation &tr,\n"
     << "                                 const mfem::Vector &u,\n"
     << "                                 const IntegrationType it) {\n"
     << "  return this->implementIntegrate(e, tr, u, it);\n"
     << "}  // end of integrate\n"
     << '\n'
     << "void " << d.name << "::updateResidual(mfem::Vector &Fe,\n"
     << "                         const mfem::FiniteElement &e,\n"
     << "                         mfem::ElementTransformation &tr,\n"
     << "                         const mfem::Vector &u) {\n"
     << "  this->implementUpdateResidual(Fe, e, tr, u);\n"
     << "}  // end of updateResidual\n"
     << '\n'
     << "  void " << d.name << "::updateJacobian(mfem::DenseMatrix &Ke,\n"
     << "                             const mfem::FiniteElement &e,\n"
     << "                             mfem::ElementTransformation &tr,\n"
     << "                             const mfem::Vector &) {\n"
     << "  this->implementUpdateJacobian(Ke, e, tr);\n"
     << "}  // end of updateJacobian\n"
     << '\n'
     << "void " << d.name << "::computeInnerForces(mfem::Vector &Fe,\n"
     << "                         const mfem::FiniteElement &e,\n"
     << "                         mfem::ElementTransformation &tr) {\n"
     << "  this->implementComputeInnerForces(Fe, e, tr);\n"
     << "}  // end of computeInnerForces\n"
     << '\n'
     << d.name << "::~" << d.name << "() = default;\n"
     << '\n'
     << "}  // end of namespace mfem_mgis\n";
}  // end of generateSourceFile

template <typename T>
T convert(const std::string_view &);

template <>
bool convert(const std::string_view &s) {
  if (s == "false") {
    return false;
  }
  if (s != "true") {
    raise("convert: can't convert '" + std::string(s) + "' to boolean value");
  }
  return true;
}  // end of convert

template <>
std::string convert(const std::string_view &s) {
  return std::string{s};
}  // end of convert

/*!
 * \brief structure in charge of parsing the command line arguments
 */
struct CommandLineArgumentParser {
  //! \brief a simple alias
  using CallBack = std::function<void(std::string_view)>;
  /*!
   * \brief parse the command line arguments
   * \param[in] argc: number of command line arguments
   * \param[in] argv: arguments
   */
  void parse(const int argc, const char *const *const argv);
  /*!
   * \brief add support for a new command line argument without option.
   * \param[in] s: short command line argument name
   * \param[in] a: long command line argument name
   * \param[in] c: action associated with the command line argument
   * \return this object
   */
  CommandLineArgumentParser &registerCommandLineArgument(
      const char s, const char *const a, const std::function<void()> &c);
  /*!
   * \brief add support for a new command line argument without option.
   * \param[in] a: long command line argument name
   * \param[in] c: action associated with the command line argument
   * \return this object
   */
  CommandLineArgumentParser &registerCommandLineArgument(
      const char *const a, const std::function<void()> &c);
  /*!
   * \brief add support for a new command line argument.
   * \param[in] s: short command line argument name
   * \param[in] a: long command line argument name
   * \param[in] c: action associated with the command line argument
   * \return this object
   */
  CommandLineArgumentParser &registerCommandLineArgument(const char s,
                                                         const char *const a,
                                                         const CallBack &c);
  /*!
   * \brief add support for a new command line argument.
   * \param[in] a: long command line argument name
   * \param[in] c: action associated with the command line argument
   * \return this object
   */
  CommandLineArgumentParser &registerCommandLineArgument(const char *const a,
                                                         const CallBack &c);
  /*!
   * \brief add support for a new command line argument which aims at setting a
   * value.
   * \param[out] v: value set by the treatment of the command line argument
   * \param[in] a: long command line argument name
   * \return this object
   */
  template <typename T>
  CommandLineArgumentParser &registerCommandLineArgument(T &v,
                                                         const char *const a);
  /*!
   * \brief add support for a new command line argument which aims at setting a
   * value.
   * \param[out] v: value set by the treatment of the command line argument
   * \param[in] s: short command line argument name
   * \param[in] a: long command line argument name
   * \return this object
   */
  template <typename T>
  CommandLineArgumentParser &registerCommandLineArgument(T &v,
                                                         const char s,
                                                         const char *const a);

 private:
  //! \brief registered call backs associated with short command line arguments
  std::map<char, CallBack> short_command_line_arguments_callbacks;
  //! \brief registered call backs associated with long command line arguments
  std::map<std::string, CallBack, std::less<>>
      long_command_line_arguments_callbacks;
};

void CommandLineArgumentParser::parse(const int argc,
                                      const char *const *const argv) {
  auto is_long_command_line_argument = [](const std::string_view &v) {
    if (v.size() < 3) {
      return false;
    }
    return (v[0] == '-') && (v[1] == '-');
  };
  auto pa = argv + 1;
  while (pa != argv + argc) {
    auto a = std::string_view{*pa};
    if (is_long_command_line_argument(a)) {
      a.remove_prefix(2);
      const auto pos = a.find_first_of('=');
      const auto n = a.substr(0, pos);
      const auto pc = this->long_command_line_arguments_callbacks.find(n);
      if (pc == this->long_command_line_arguments_callbacks.end()) {
        raise("invalid command line argument '" + std::string(n) + "'");
      }
      if (pos == std::string_view::npos) {
        pc->second("");
      } else {
        a.remove_prefix(pos + 1);
        if (a.empty()) {
          raise("no option given to the '" + std::string(n) +
                "' command line argument");
        }
        pc->second(a);
      }
    } else if (a[0] == '-') {
      if (a.size() != 2) {
        raise("invalid command line argument '" + std::string(a) + "'");
      }
      const auto n = a[1];
      const auto pc = this->short_command_line_arguments_callbacks.find(n);
      if (pc == this->short_command_line_arguments_callbacks.end()) {
        raise("invalid command line argument '" + std::string(1, n) + "'");
      }
      const auto pna = pa + 1;  // next argument
      if ((pna != argv + argc) && (*pna[0] != '-')) {
        ++pa;
        a = *pa;
        pc->second(a);
      } else {
        pc->second("");
      }
    } else {
      raise("invalid command line argument '" + std::string(a) + "'");
    }
    ++pa;
  }
}

CommandLineArgumentParser &
CommandLineArgumentParser::registerCommandLineArgument(
    const char *const a, const std::function<void()> &c) {
  auto action = [n = std::string(a), c](std::string_view o) {
    if (!o.empty()) {
      raise("command line argument '" + n + "' does not accept any option");
    }
    c();
  };
  this->long_command_line_arguments_callbacks.insert({a, action});
  return *this;
}  // end of registerCommandLineArgument

CommandLineArgumentParser &
CommandLineArgumentParser::registerCommandLineArgument(
    const char s, const char *const a, const std::function<void()> &c) {
  auto action = [n = std::string(a), s, c](std::string_view o) {
    if (!o.empty()) {
      raise("command line argument '" + n + "' (" + std::string(1, s) +
            ") does not accept any option");
    }
    c();
  };
  this->short_command_line_arguments_callbacks.insert({s, action});
  this->long_command_line_arguments_callbacks.insert({a, action});
  return *this;
}  // end of registerCommandLineArgument

CommandLineArgumentParser &
CommandLineArgumentParser::registerCommandLineArgument(const char *const a,
                                                       const CallBack &c) {
  auto action = [n = std::string(a), c](std::string_view o) {
    if (o.empty()) {
      raise("no option given to the '" + n + "' command line argument");
    }
    c(o);
  };
  this->long_command_line_arguments_callbacks.insert({a, action});
  return *this;
}  // end of registerCommandLineArgument

CommandLineArgumentParser &
CommandLineArgumentParser::registerCommandLineArgument(const char s,
                                                       const char *const a,
                                                       const CallBack &c) {
  auto action = [n = std::string(a), c](std::string_view o) {
    if (o.empty()) {
      raise("no option given to the '" + n + "' command line argument");
    }
    c(o);
  };
  this->short_command_line_arguments_callbacks.insert({s, action});
  this->long_command_line_arguments_callbacks.insert({a, action});
  return *this;
}  // end of registerCommandLineArgument

template <typename T>
CommandLineArgumentParser &
CommandLineArgumentParser::registerCommandLineArgument(T &v,
                                                       const char s,
                                                       const char *const a) {
  auto action = [&v, this](std::string_view o) { v = convert<T>(o); };
  this->registerCommandLineArgument(s, a, action);
  return *this;
}  // end of registerCommandLineArgument

template <typename T>
CommandLineArgumentParser &
CommandLineArgumentParser::registerCommandLineArgument(T &v,
                                                       const char *const a) {
  auto action = [&v, this](std::string_view s) { v = convert<T>(s); };
  this->registerCommandLineArgument(a, action);
  return *this;
}  // end of registerCommandLineArgument

int main(const int argc, const char *const *const argv) {
  using Generators =
      std::map<std::string, BehaviourIntegratorDescription::BMatrixGenerator,
               std::less<>>;
  std::map<std::string, Generators, std::less<>> generators;
  auto insert = [&generators](
                    const char *const h, const char *const n,
                    const BehaviourIntegratorDescription::BMatrixGenerator &g) {
    generators[h].insert({n, g});
  };
  insert("Tridimensional", "StandardFiniteStrainMechanics",
         &makeTridimensionalFiniteStrainMechanicsBMatrix);
  insert("Tridimensional", "StandardSmallStrainMechanics",
         &makeTridimensionalSmallStrainMechanicsBMatrix);
  insert("PlaneStrain", "StandardFiniteStrainMechanics",
         &makePlaneStrainFiniteStrainMechanicsBMatrix);
  insert("PlaneStrain", "StandardSmallStrainMechanics",
         &makePlaneStrainSmallStrainMechanicsBMatrix);
  insert("PlaneStress", "StandardFiniteStrainMechanics",
         &makePlaneStressFiniteStrainMechanicsBMatrix);
  insert("PlaneStress", "StandardSmallStrainMechanics",
         &makePlaneStressSmallStrainMechanicsBMatrix);
  //
  insert("PlaneStrain", "StationaryNonLinearHeatTransfer",
         &makeHeatTransferBMatrix<2u>);
  insert("PlaneStress", "StationaryNonLinearHeatTransfer",
         &makeHeatTransferBMatrix<2u>);
  insert("Tridimensional", "StationaryNonLinearHeatTransfer",
         &makeHeatTransferBMatrix<3u>);
  //
  auto d = BehaviourIntegratorDescription{};
  //
  auto generator = std::string{};
  auto generate_header_file = false;
  auto generate_source_file = false;
  CommandLineArgumentParser args_parser;
  args_parser
      .registerCommandLineArgument(
          's', "symmetry",
          [&d](std::string_view o) {
            if (o == "orthotropic") {
              d.isotropic = false;
            } else if (o == "isotropic") {
              d.isotropic = true;
            } else {
              raise(
                  "invalid value for the --symmetry (-s) command line "
                  "argument. "
                  "Expected 'isotropic' or 'orthotropic'.");
            }
          })
      .registerCommandLineArgument("header-file",
                                   [&]() { generate_header_file = true; })
      .registerCommandLineArgument("source-file",
                                   [&]() { generate_source_file = true; })
      .registerCommandLineArgument("generator",
                                   [&generator](std::string_view g) {
                                     if (!generator.empty()) {
                                       raise("generator already set");
                                     }
                                     generator = g;
                                   })
      .registerCommandLineArgument("hypothesis",
                                   [&d](std::string_view h) {
                                     if (!d.hypothesis.empty()) {
                                       raise("hypothesis already set");
                                     }
                                     d.hypothesis = h;
                                   })
      .registerCommandLineArgument("unknown-name",
                                   [&d](std::string_view n) {
                                     if (!d.unknown_name.empty()) {
                                       raise("unknown name already set");
                                     }
                                     d.unknown_name = n;
                                   })
      .registerCommandLineArgument(
          "requires-unknown-values-as-external-state-variable", [&d]() {
            d.requires_unknown_value_as_external_state_variable = true;
          });
  args_parser.parse(argc, argv);
  //
  if (generator.empty()) {
    raise("no generator set");
  }
  if (d.hypothesis.empty()) {
    raise("no hypothesis set");
  }
  const auto ph = generators.find(d.hypothesis);
  if (ph == generators.end()) {
    raise("unsupported hypothesis '" + d.hypothesis + "'");
  }
  const auto pg = ph->second.find(generator);
  if (pg == ph->second.end()) {
    raise("undefined generator '" + generator + "' for hypothesis '" +
          d.hypothesis + "'");
  }
  d.generator = pg->second;
  d.generator_name = generator;
  // source code generation
  if ((!generate_header_file) && (!generate_source_file)) {
    raise("no file to be generated");
  }
  d.name =
      d.isotropic
          ? "Isotropic" + d.hypothesis + generator + "BehaviourIntegrator"
          : "Orthotropic" + d.hypothesis + generator + "BehaviourIntegrator";
  if (d.requires_unknown_value_as_external_state_variable) {
    if (d.unknown_name.empty()) {
      raise("unknown name must be specified\n");
    }
  }
  if (generate_header_file) {
    std::ofstream out(d.name + ".hxx");
    generateHeaderFile(out, d);
  }
  if (generate_source_file) {
    std::ofstream out(d.name + ".cxx");
    generateSourceFile(out, d);
  }
  return EXIT_SUCCESS;
}
