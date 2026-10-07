/*!
 * \file   src/PartialQuadratureFunctionsSet.cxx
 * \brief  This file implements the `PartialQuadratureFunctionsSet` class
 * \author Thomas Helfer
 * \date   02/06/2025
 */

#include <algorithm>
#include "MFEMMGIS/PartialQuadratureFunctionsSet.hxx"

namespace mfem_mgis {

  static std::vector<std::shared_ptr<PartialQuadratureFunction>>
  buildPartialQuadratureFunctionsSet(
      attributes::Throwing,
      const std::vector<std::shared_ptr<const PartialQuadratureSpace>>& qspaces,
      const mfem_mgis::size_type n) {
    auto functions = std::vector<std::shared_ptr<PartialQuadratureFunction>>{};
    functions.reserve(qspaces.size());
    for (const auto& qspace : qspaces) {
      functions.push_back(
          std::make_shared<PartialQuadratureFunction>(qspace, n));
    }
    return functions;
  }  // end of functions

  PartialQuadratureFunctionsSet::PartialQuadratureFunctionsSet(
      const std::vector<std::shared_ptr<const PartialQuadratureSpace>>& qspaces,
      const mfem_mgis::size_type n)
      : PartialQuadratureFunctionsSet(
            buildPartialQuadratureFunctionsSet(throwing, qspaces, n)) {}

  PartialQuadratureFunctionsSet::PartialQuadratureFunctionsSet(
      const std::vector<std::shared_ptr<PartialQuadratureFunction>>& functions)
      : std::vector<std::shared_ptr<PartialQuadratureFunction>>(functions) {
    auto locations = std::vector<LocationIdentifier>{};
    locations.reserve(this->size());
    for (const auto& fptr : *this) {
      if (fptr.get() == nullptr) {
        raise("invalid function");
      }
      const auto& qspace = fptr->getPartialQuadratureSpace();
      const auto l = qspace.getLocation();
      if (std::find(locations.begin(), locations.end(), l) != locations.end()) {
        raise("multiple functions defined on " + getLocationDescription(l));
      }
      locations.push_back(l);
    }
  }  // end of PartialQuadratureFunctionsSet

  std::vector<std::shared_ptr<const PartialQuadratureFunction>>
  PartialQuadratureFunctionsSet::getFunctions() const {
    return std::vector<std::shared_ptr<const PartialQuadratureFunction>>(
        this->begin(), this->end());
  }  // end of getFunctions

  const std::vector<std::shared_ptr<PartialQuadratureFunction>>&
  PartialQuadratureFunctionsSet::getFunctions() {
    return *this;
  }  // end of getFunctions

  std::vector<LocationIdentifier> PartialQuadratureFunctionsSet::getLocations()
      const noexcept {
    auto locations = std::vector<LocationIdentifier>{};
    locations.reserve(this->size());
    for (const auto& fptr : *this) {
      locations.push_back(fptr->getPartialQuadratureSpace().getLocation());
    }
    return locations;
  }  // end of getMaterialIdentifiers

  std::shared_ptr<PartialQuadratureFunction> PartialQuadratureFunctionsSet::get(
      Context& ctx, const LocationIdentifier l) noexcept {
    for (const auto& fptr : *this) {
      const auto fl = fptr->getPartialQuadratureSpace().getLocation();
      if (fl == l) {
        return fptr;
      }
    }
    return ctx.registerErrorMessage("no function associated with " +
                                    getLocationDescription(l) + " found");
  }  // end of get

  std::shared_ptr<const PartialQuadratureFunction>
  PartialQuadratureFunctionsSet::get(
      Context& ctx, const LocationIdentifier l) const noexcept {
    for (const auto& fptr : *this) {
      const auto fl = fptr->getPartialQuadratureSpace().getLocation();
      if (l == fl) {
        return fptr;
      }
    }
    return ctx.registerErrorMessage("no function associated with " +
                                    getLocationDescription(l) + " found");
  }  // end of get

  bool PartialQuadratureFunctionsSet::update(Context& ctx, UpdateFunction& f) {
    for (const auto& fptr : *this) {
      if (!f(ctx, *fptr)) {
        return false;
      }
    }
    return true;
  }

  void PartialQuadratureFunctionsSet::update(UpdateFunction2& f) {
    for (const auto& fptr : *this) {
      f(*fptr);
    }
  }  // end of update

}  // end of namespace mfem_mgis
