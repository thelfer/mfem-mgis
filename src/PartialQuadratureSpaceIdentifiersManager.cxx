/*!
 * \file   PartialQuadratureSpaceIdentifiersManager.cxx
 * \brief  This file implements the `PartialQuadratureSpaceIdentifiersManager`
 * class
 * \author Thomas Helfer
 * \date   29/03/2026
 */

#include <cassert>
#include "MFEMMGIS/PartialQuadratureSpaceIdentifiersManager.hxx"

namespace mfem_mgis {

  [[nodiscard]] static size_type getIdentifier(
      std::vector<std::vector<std::shared_ptr<const PartialQuadratureSpace>>>&
          identifiers,
      const std::shared_ptr<const PartialQuadratureSpace>& s) noexcept {
    // check if the quadrature space has already been declared
    for (size_type id = 0; const auto& spaces : identifiers) {
      for (const auto& s2 : spaces) {
        if (s == s2) {
          return id;
        }
      }
    }
    // add the quadrature space
    for (size_type id = 0; auto& spaces : identifiers) {
      assert(!spaces.empty());  // by design, spaces can't be empty
      if (areEquivalent(*s, *(spaces.at(0)))) {
        spaces.push_back(s);
        return id;
      }
    }
    // the quadrature space is not equivalent to any of the previously registred
    // space
    const auto id = identifiers.size();
    identifiers.push_back(
        std::vector<std::shared_ptr<const PartialQuadratureSpace>>(1, s));
    return id;
  }  // end of getIdentifier

  PartialQuadratureSpaceIdentifiersManager::
      PartialQuadratureSpaceIdentifiersManager(
          const MeshDiscretization& m) noexcept
      : mesh(m) {}  // end of PartialQuadratureSpaceIdentifiersManager

  std::optional<size_type>
  PartialQuadratureSpaceIdentifiersManager::getIdentifier(
      Context& ctx,
      const std::shared_ptr<const PartialQuadratureSpace>& s) const noexcept {
    using PerMaterialIdentifiersList =
        std::vector<std::vector<std::shared_ptr<const PartialQuadratureSpace>>>;
    if (s.get() == nullptr) {
      return ctx.registerErrorMessage("invalid partial quadrature space");
    }
    if (s->getMeshDiscretization() != this->mesh) {
      return ctx.registerErrorMessage(
          "partial quadrature space is not defined on the proper mesh");
    }
    const auto mid = s->getId();
    auto p = this->ids.find(mid);
    if (p == this->ids.end()) {
      p = this->ids.insert({mid, PerMaterialIdentifiersList{}}).first;
    }
    return ::mfem_mgis::getIdentifier(p->second, s);
  }  // end of getIdendifier

  bool PartialQuadratureSpaceIdentifiersManager::areEquivalent(
      const std::shared_ptr<const PartialQuadratureSpace>& qspace1,
      const std::shared_ptr<const PartialQuadratureSpace>& qspace2)
      const noexcept {
    auto ctx = Context{};
    const auto oid1 = this->getIdentifier(ctx, qspace1);
    const auto oid2 = this->getIdentifier(ctx, qspace2);
    if (!areValid(oid1, oid2)) {
      return false;
    }
    return *oid1 == *oid2;
  }  // end of areEquivalent

  PartialQuadratureSpaceIdentifiersManager::
      ~PartialQuadratureSpaceIdentifiersManager() noexcept = default;

}  // end of namespace mfem_mgis
