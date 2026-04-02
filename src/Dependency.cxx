/*!
 * \file   Dependency.cxx
 * \brief  This file implements the DependencyBase
 * \date   01/04/2026
 */

#include "MFEMMGIS/Provider.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Dependency.hxx"

namespace mfem_mgis {

  DependencyBase::DependencyBase(std::string_view n,
                                 const DependencyStatus s) noexcept
      : status(s), name(n) {}  // end of DependencyBase

  DependencyBase::DependencyBase(const DependencyBase &) noexcept = default;

  DependencyBase::DependencyBase(DependencyBase &&) noexcept = default;

  DependencyBase::~DependencyBase() noexcept = default;

  std::string DependencyBase::getSpecificationsAsString() const noexcept {
    const auto nc = this->getNumberOfComponents();
    auto s = std::string{"("};
    if (nc.has_value()) {
      s += std::to_string(*nc);
    } else {
      s += "unspecified number of components";
    }
    s += ')';
    return s;
  }  // end of getSpecificationsAsString

  bool DependencyBase::hasConcreteSpecifications() const noexcept {
    return this->number_of_components.has_value();
  }  // end of hasConcreteSpecifications

  bool DependencyBase::matchesSpecifications(
      Context &ctx, const size_type nc) const noexcept {
    if (this->number_of_components.has_value()) {
      if (*(this->number_of_components) != nc) {
        return ctx.registerErrorMessage(
            "inconsistent number of components for " + this->getDescription());
      }
    }
    return true;
  }  // end of matchesSpecifications

  bool DependencyBase::checkSpecifications(Context &ctx,
                                           const size_type nc) const noexcept {
    if (!this->hasConcreteSpecifications()) {
      return ctx.registerErrorMessage(
          "no concrete specifications for dependency '" + this->name + "'");
    }
    return this->matchesSpecifications(ctx, nc);
  }  // end of checkSpecifications

  const std::string &DependencyBase::getName() const noexcept {
    return this->name;
  }  // end of getName

  std::optional<size_type> DependencyBase::getNumberOfComponents()
      const noexcept {
    return this->number_of_components;
  }  // end of getNumRows

  bool DependencyBase::setProvider(Context &ctx, const Provider &p) noexcept {
    if ((this->hasProvider()) && (this->provider != &p)) {
      return ctx.registerErrorMessage(
          "'" + p.getIdentifier() +
          "' can't be declared as the provider of the " +
          this->getDescription() + " as another provider named '" +
          this->provider->getIdentifier() + "' has already been declared");
    }
    this->provider = &p;
    return true;
  }  // end of setProvider

  bool DependencyBase::hasProvider() const noexcept {
    return this->provider != nullptr;
  }  // end of hasProvider

  OptionalReference<const Provider> DependencyBase::getProvider(
      Context &ctx) const noexcept {
    if (!this->hasProvider()) {
      return ctx.registerErrorMessage("no provider defined for " +
                                      this->getDescription());
    }
    return {this->provider};
  }  // end of getProvider

  bool DependencyBase::isRequired() const noexcept {
    return this->status == DependencyStatus::REQUIRED;
  }  // end of isRequired

  bool DependencyBase::isOptional() const noexcept {
    return this->status == DependencyStatus::OPTIONAL;
  }  // end of isOptional

  bool DependencyBase::checkAndUpdate(Context &ctx,
                                       const DependencyBase &d) noexcept {
    if (d.number_of_components.has_value()) {
      if (this->number_of_components.has_value()) {
        if (*(this->number_of_components) != *(d.number_of_components)) {
          return ctx.registerErrorMessage("inconsistent number of components");
        }
      } else {
        this->number_of_components = *(d.number_of_components);
      }
    }
    return true;
  }  // end of update

  bool DependencyBase::checkAndUpdate(Context &ctx,
                                      const size_type nc) noexcept {
    if (this->number_of_components.has_value()) {
      if (*(this->number_of_components) != nc) {
        return ctx.registerErrorMessage("inconsistent number of rows");
      }
    }
    this->number_of_components = nc;
    return true;
  }  // end of checkAndUpdate_

  bool DependencyBase::checkAndUpdate(Context &ctx,
                                      std::string_view n,
                                      const size_type nc) noexcept {
    if (this->name != n) {
      return ctx.registerErrorMessage("inconsistent name");
    }
    return this->checkAndUpdate(ctx, nc);
  }  // end of checkAndUpdate

  bool QPDependency::reportProviderRequiresQuadratureIdToBeDefined(
      Context &ctx, const QPDependency &d, const std::string &n) noexcept {
    return ctx.registerErrorMessage("the provider '" + n +
                                    "' can only provide the dependency "
                                    "on integration points '" +
                                    d.getName() + "' on material '" +
                                    std::to_string(d.getMaterialIdentifier()) +
                                    "' if the quadradure is defined.");
  }  // end of reportProviderRequiresQuadratureIdToBeDefined

  QPDependency::QPDependency(const size_type m,
                             std::string_view n,
                             const DependencyStatus s) noexcept
      : DependencyBase(n, s), material_identifier(m) {}  // end of QPDependency

  QPDependency::QPDependency(std::shared_ptr<const PartialQuadratureSpace> s,
                             std::string_view n,
                             const DependencyStatus ds) noexcept
      : DependencyBase(n, ds),
        material_identifier(s->getId()),
        qspace(s) {}  // end of QPDependency

  QPDependency::QPDependency(const QPDependency &) noexcept = default;

  QPDependency::QPDependency(QPDependency &&) noexcept = default;

  size_type QPDependency::getMaterialIdentifier() const noexcept {
    return this->material_identifier;
  }  // end of getMaterialIdentifier

  OptionalReference<const PartialQuadratureSpace>
  QPDependency::getPartialQuadratureSpace() const noexcept {
    return {this->qspace.get()};
  }  // end of getQuadratureId

  std::shared_ptr<const PartialQuadratureSpace>
  QPDependency::getPartialQuadratureSpacePointer() const noexcept {
    return this->qspace;
  }  // end of getQuadratureId

  bool QPDependency::setPartialQuadratureSpace(
      Context &ctx, std::shared_ptr<const PartialQuadratureSpace> s) noexcept {
    if (s.get() == nullptr) {
      return ctx.registerErrorMessage("invalid partial quadrature space");
    }
    if (s->getId() != this->material_identifier) {
      return ctx.registerErrorMessage("inconsistent material identifier");
    }
    if (isValid(this->qspace)) {
      return ctx.registerErrorMessage("partial quadrature space already defined");
    }
    this->qspace = s;
    return true;
  }

  bool QPDependency::setDuplicate(Context &ctx) noexcept {
    if (this->is_duplicate) {
      return ctx.registerErrorMessage(this->getDescription() +
                                      " can't be marked duplicate twice");
    }
    this->is_duplicate = true;
    return true;
  }  // end of setDuplicate

  bool QPDependency::isDuplicate() const noexcept {
    return this->is_duplicate;
  }  // end of isDuplicate

  std::string QPDependency::getDescription() const noexcept {
    const auto d = [this] {
      if (this->isRequired()) {
        return "required dependency '" + this->getName() + "' on material '" +
               std::to_string(this->material_identifier) + "'";
      }
      return "optional dependency '" + this->getName() + "' on material '" +
             std::to_string(this->material_identifier) + "'";
    }();
    if (isValid(this->qspace)) {
      return d + " for specified quadrature space";
    }
    return d + " on unspecified quadrature";
  }  // end of getDescription

  bool QPDependency::checkAndUpdate(Context &ctx, const size_type nc) noexcept {
    return DependencyBase::checkAndUpdate(ctx, nc);
  }  // end of checkAndUpdate

  bool QPDependency::checkAndUpdate(Context &ctx,
                                    const QPDependency &d) noexcept {
    return DependencyBase::checkAndUpdate(ctx, d);
  }  // end of checkAndUpdate

  QPDependency::~QPDependency() noexcept = default;

}  // end of namespace mfem_mgis
