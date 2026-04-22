/*!
 * \file   DependenciesManager.cxx
 * \brief  This file implements the `DependenciesManager` class
 * \author Thomas Helfer
 * \date   02/04/2026
 */

#include <algorithm>
#include "MFEMMGIS/Provider.hxx"
#include "MFEMMGIS/Dependency.hxx"
#include "MFEMMGIS/DependenciesManager.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PartialQuadratureSpaceIdentifiersManager.hxx"

namespace mfem_mgis {

  std::string DependenciesManager::getLocationDescription(
      const QPDependency &d, const TimeStepStage s) noexcept {
    const auto tss = [&s]() -> std::string {
      if (s == ets) {
        return "end";
      }
      return "beginning";
    };
    const auto qspace = d.getPartialQuadratureSpace();
    auto r = std::string{};
    if (isValid(qspace)) {
      r += "material '" + qspace->getMaterialName() + "' ";
      r += "for given quadrature ";
    } else {
      const auto mid = d.getMaterialIdentifier();
      r += "material '" + std::to_string(mid) + "' ";
      r += "on unspecified quadrature ";
    }
    return r + "at the " + tss() + " of the time step";
  }  // end of getLocationDescription

  std::vector<QPDependency> &DependenciesManager::getLocalQPDependenciesManager(
      const size_type m, const TimeStepStage s) noexcept {
    if (s == bts) {
      return this->registeredQPDependencies[0][m];
    }
    return this->registeredQPDependencies[1][m];
  }  // end of getLocalQPDependenciesManager

  DependenciesManager::DependenciesManager(
      const PartialQuadratureSpaceIdentifiersManager &m) noexcept
      : qids(m) {}  // end of DependenciesManager

  //   //   bool DependenciesManager::declareDependency(
  //   //       Context &ctx,
  //   //       const MeshSet &m,
  //   //       const QuadratureSpace &qspace,
  //   //       const TimeStepStage s,
  //   //       const QPEvaluatorDescription &d) noexcept {
  //   //     return this->declareDependency(ctx, s, QPDependency{m, qspace,
  //   d});
  //   //   }  // end of declareDependency

  bool DependenciesManager::declareDependency(Context &ctx,
                                              const TimeStepStage s,
                                              const QPDependency &d) noexcept {
    if (d.isDuplicate()) {
      return ctx.registerErrorMessage(
          d.getDescription() +
          " can't be registered as it is flagged as being a duplicate");
    }
    auto nd = d;
    auto &ldm =
        this->getLocalQPDependenciesManager(nd.getMaterialIdentifier(), s);
    // we first make a check regardless of the quadrature id
    auto p = std::find_if(ldm.begin(), ldm.end(), [&nd](const auto &dep) {
      return nd.getName() == dep.getName();
    });
    if (p != ldm.end()) {
      auto report_inconsistency_check = [&ctx, &s, &nd]() -> bool {
        return ctx.registerErrorMessage(
            "inconsistent dependency declaration '" + nd.getName() + "' on " +
            DependenciesManager::getLocationDescription(nd, s));
      };
      // consistency checks
      if (!p->checkAndUpdate(ctx, nd)) {
        return report_inconsistency_check();
      }
      if (!nd.checkAndUpdate(ctx, *p)) {
        return report_inconsistency_check();
      }
    }
    const auto pd = [this, &nd, &ldm] {
      const auto qspace = nd.getPartialQuadratureSpacePointer();
      if (isValid(qspace)) {
        return std::find_if(ldm.begin(), ldm.end(),
                            [this, &nd, qspace](const auto &dep) {  //
                              const auto q = dep.getPartialQuadratureSpacePointer();
                              if (isInvalid(q)) {
                                return false;
                              }
                              return (nd.getName() == dep.getName()) &&
                                     (this->qids.areEquivalent(qspace, q));
                            });
      } else {
        return std::find_if(ldm.begin(), ldm.end(), [&nd](const auto &dep) {  //
          const auto q = dep.getPartialQuadratureSpacePointer();
          if (isValid(q)) {
            return false;
          }
          return nd.getName() == dep.getName();
        });
      }
    }();
    if (pd == ldm.end()) {
      ctx.debug("declaring dependency at integration points '" + nd.getName() +
                "' on " + DependenciesManager::getLocationDescription(nd, s));
      ldm.push_back(nd);
    }
    return true;
  }  // end of declareDependency

  bool DependenciesManager::setProvider(
      Context &ctx,
      const Provider &pr,
      const QPDependency &d,
      std::shared_ptr<const PartialQuadratureSpace> qspace,
      const size_type nc,
      const TimeStepStage ts) noexcept {
    if (isInvalid(qspace)) {
      return ctx.registerErrorMessage("invalid partial quadrature space");
    }
    if (!d.matchesSpecifications(ctx, nc)) {
      return false;
    }
    const auto &m = d.getMaterialIdentifier();
    const auto dqspace = d.getPartialQuadratureSpacePointer();
    if (isValid(dqspace)) {
      if (!this->qids.areEquivalent(dqspace, qspace)) {
        return ctx.registerErrorMessage(
            "inconsistent partial quadrature spaces");
      }
    }
    auto &dependencies_container =
        [&ts, this]() -> std::map<size_type, std::vector<QPDependency>> & {
      if (ts == bts) {
        return this->registeredQPDependencies[0];
      }
      return this->registeredQPDependencies[1];
    }();
    if (!dependencies_container.contains(m)) {
      return ctx.registerErrorMessage(
          "dependency is not handled by this manager");
    }
    auto &ldm = dependencies_container.at(m);
    auto pd = ldm.end();
    for (auto p = ldm.begin(); p != ldm.end(); ++p) {
      if (p->isDuplicate()) {
        // duplicated dependencies arise when dependencies with initially
        // no specific partial quadrature space are resolved
        continue;
      }
      auto &ld = *p;
      if (d.getName() == ld.getName()) {
        if (!ld.checkAndUpdate(ctx, nc)) {
          return false;
        }
        const auto oldqspace = ld.getPartialQuadratureSpacePointer();
        if ((isValid(oldqspace)) && (isValid(dqspace))) {
          if (this->qids.areEquivalent(oldqspace, qspace)) {
            if (pd != ldm.end()) {
              // using message error to signal an internal error that is
              // currently not being testable. This is a paranoïc check
              mfem_mgis::abort("internal error: " + d.getDescription() +
                               " is registered twice");
            }
            if (!ld.setProvider(ctx, pr)) {
              return false;
            }
            pd = p;
          }
        } else if ((isInvalid(oldqspace)) && (isInvalid(dqspace))) {
          if (pd != ldm.end()) {
            // Using message error to signal an internal error that is currently
            // not being testable. This is a paranoïc check
            mfem_mgis::abort("internal error: " + d.getDescription() +
                             " is registered twice");
          }
          if (!ld.setProvider(ctx, pr)) {
            return false;
          }
          pd = p;
        }
      }
    }
    if (pd == ldm.end()) {
      return ctx.registerErrorMessage(
          "dependency is not handled by this manager");
    }
    if (isInvalid(dqspace)) {
      if (!pd->setPartialQuadratureSpace(ctx, qspace)) {
        return false;
      }
      // here we have potentially created a duplicate with another dependency
      // who have specified its partial quadrature space.
      // This duplication may create a problem, so we will mark our dependency
      // as duplicate if this case happens
      for (auto p = ldm.begin(); p != ldm.end(); ++p) {
        if ((p == pd) || (p->getName() != pd->getName()) ||
            (p->getMaterialIdentifier() != pd->getMaterialIdentifier())) {
          continue;
        }
        if (isInvalid(p->getPartialQuadratureSpace())) {
          // if p has no specified quadrature id, this means that p and pd
          // points to the same dependency this shall not happen, because the
          // declareDepency method should have avoided this case
          mfem_mgis::abort("internal error: " + pd->getDescription() +
                           " declared twice");
        }
        // we don't check if *p is a duplicate.
        // if it is, the following code will be executed by the original anyway
        if (this->qids.areEquivalent(p->getPartialQuadratureSpacePointer(),
                                     qspace)) {
          if (!pd->setDuplicate(ctx)) {
            return false;
          }
          break;
        }
      }
    }
    return true;
  }  // end of setProvider

  DependenciesManager::DependenciesAnalysisOutput
  DependenciesManager::analyseDependencies(
      const AnalyseDependenciesFilter f) const noexcept {
    auto r = DependenciesAnalysisOutput{};
    auto exe = [&f]<typename DependencyType>(
                   const std::map<size_type, std::vector<DependencyType>> &in) {
      auto out = std::vector<DependencyType>{};
      for (const auto &[m, deps] : in) {
        static_cast<void>(m);
        for (const auto &d : deps) {
          if (d.hasProvider()) {
            continue;
          }
          if (d.isRequired()) {
            if ((f == AnalyseDependenciesFilter::ONLY_REQUIRED) ||
                (f == AnalyseDependenciesFilter::ALL)) {
              out.push_back(d);
            }
          } else {
            if ((f == AnalyseDependenciesFilter::ONLY_OPTIONAL) ||
                (f == AnalyseDependenciesFilter::ALL)) {
              out.push_back(d);
            }
          }
        }
      }
      return out;
    };
    r.missingQPDependencies_bts = exe(this->registeredQPDependencies[0]);
    r.missingQPDependencies_ets = exe(this->registeredQPDependencies[1]);
    return r;
  }  // end of analyseDependencies

  bool DependenciesManager::resolveDependencies(
      Context &ctx, QPEvaluatorsFactory &f) const noexcept {
    auto r = DependenciesAnalysisOutput{};
    auto exe = [&ctx, &f](
                   const std::map<size_type, std::vector<QPDependency>> &in,
                   const TimeStepStage ts) -> bool {
      for (const auto &[m, deps] : in) {
        static_cast<void>(m);
        // At this stage, the list of dependencies on input can contain
        // duplicates. This is due to the that some dependencies did initially
        // specify their quadrature id. if we keep duplicates
        // `resolveDependency` will be called twice, we will results in
        // multiples registrations of the ip evaluator generator, leading to an
        // error
        for (const auto &d : deps) {
          if (!d.hasProvider()) {
            if (d.isRequired()) {
              return ctx.registerErrorMessage(d.getDescription() +
                                              " has no provider");
            }
            continue;
          }
          if (d.isDuplicate()) {
            // skipping duplicated dependency
            // the resolveDepency method will be called on the "original"
            // one
            continue;
          }
          const auto p = d.getProvider(ctx);
          if (isInvalid(p)) {
            return false;
          }
          if (!(*p)->resolveDependency(ctx, f, d, ts)) {
            return false;
          }
        }
      }
      return true;
    };
    if (!exe(this->registeredQPDependencies[0], bts)) {
      return false;
    }
    if (!exe(this->registeredQPDependencies[1], ets)) {
      return false;
    }
    return true;
  }  // end of resolveDependencies

  std::pair<size_type, std::string> getDescription(
      const DependenciesManager::DependenciesAnalysisOutput &a) noexcept {
    auto deps = std::string{};
    auto ndeps = size_type{};
    auto get_ip_dependencies = [&deps, &ndeps](const auto &mds,
                                               const TimeStepStage ts) {
      for (const auto &d : mds) {
        deps += "\n- '" + d.getName() + "' " + d.getSpecificationsAsString() +
                " on " + DependenciesManager::getLocationDescription(d, ts);
        ++ndeps;
      }
    };
    get_ip_dependencies(a.missingQPDependencies_bts, bts);
    get_ip_dependencies(a.missingQPDependencies_ets, ets);
    return {ndeps, deps};
  }  // end of getDescription

}  // end of namespace mfem_mgis
