/*!
 * \file   src/PartialQuadratureFunctionsSet.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/06/2025
 */

#include <algorithm>
#include "MFEMMGIS/PartialQuadratureFunctionsSet.hxx"

namespace mfem_mgis {

  static std::vector<std::shared_ptr<PartialQuadratureFunction>>
  buildPartialQuadratureFunctionsSet(
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
            buildPartialQuadratureFunctionsSet(qspaces, n)) {}

  PartialQuadratureFunctionsSet::PartialQuadratureFunctionsSet(
      const std::vector<std::shared_ptr<PartialQuadratureFunction>>& functions)
      : std::vector<std::shared_ptr<PartialQuadratureFunction>>(functions) {
    auto mids = std::vector<mfem_mgis::size_type>{};
    mids.reserve(this->size());
    for (const auto& fptr : *this) {
      if (fptr.get() == nullptr) {
        raise("invalid function");
      }
      const auto& qspace = fptr->getPartialQuadratureSpace();
      if (std::find(mids.begin(), mids.end(), qspace.getId()) != mids.end()) {
        raise("multiple functions defined on material '" +
              std::to_string(qspace.getId()) + "'");
      }
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

  std::vector<size_type> PartialQuadratureFunctionsSet::getMaterialIdentifiers()
      const {
    auto mids = std::vector<mfem_mgis::size_type>{};
    mids.reserve(this->size());
    for (const auto& fptr : *this) {
      mids.push_back(fptr->getPartialQuadratureSpace().getId());
    }
    return mids;
  }  // end of getMaterialIdentifiers

  std::shared_ptr<PartialQuadratureFunction> PartialQuadratureFunctionsSet::get(
      Context& ctx, const mfem_mgis::size_type m) {
    for (const auto& fptr : *this) {
      const auto mid = fptr->getPartialQuadratureSpace().getId();
      if (m == mid) {
        return fptr;
      }
    }
    return ctx.registerErrorMessage("no function associated with material '" +
                                    std::to_string(m) + "' found");
  }  // end of get

  std::shared_ptr<const PartialQuadratureFunction>
  PartialQuadratureFunctionsSet::get(Context& ctx,
                                     const mfem_mgis::size_type m) const {
    for (const auto& fptr : *this) {
      const auto mid = fptr->getPartialQuadratureSpace().getId();
      if (m == mid) {
        return fptr;
      }
    }
    return ctx.registerErrorMessage("no function associated with material '" +
                                    std::to_string(m) + "' found");
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