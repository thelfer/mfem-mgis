/*!
 * \file   src/CouplingSchemeBase.cxx
 * \brief  This file implements the `CouplingSchemeBase` class
 * \date   05/12/2022
 */

#include "MFEMMGIS/AbstractModel.hxx"
#include "MFEMMGIS/CouplingSchemeBase.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  CouplingSchemeBase::getParametersDescription() noexcept {
    //    auto d = getCouplingItemParametersDescription();
    //     d.insert(
    //         {"PrintResourcesUsage",
    //          "boolean stating if a coarse grain profiling of resources usage
    //          shall " "be displayed on the standard output at each iteration
    //          (false by " "default)"});
    //     d.insert({"CouplingItems", "list of coupling items"});
    return {};
  }  // end of getParametersDescription

  CouplingSchemeBase::ContextState CouplingSchemeBase::update(
      Context &ctx, AbstractCouplingItem &m) noexcept {
    ContextState s;
    s.verbosity_level = ctx.getVerbosityLevel();
    s.log_stream = ctx.getLogStreamPointer();
    auto log = m.getLogStreamPointer();
    if (log.get() != nullptr) {
      ctx.setLogStream(log);
    }
    const auto l = m.getVerbosityLevel();
    // only update the verbosity if the coupling item is more verbose
    if (l >= ctx.getVerbosityLevel()) {
      ctx.setVerbosityLevel(l);
    }
    return s;
  }  // end of update

  void CouplingSchemeBase::restore(Context &ctx,
                                   const ContextState &s) noexcept {
    ctx.setVerbosityLevel(s.verbosity_level);
    ctx.setLogStream(s.log_stream);
  }  // end of restore

  CouplingSchemeBase::CouplingSchemeBase() = default;

  std::vector<std::string> CouplingSchemeBase::getLocations() const noexcept {
    return {};
  }  // end of getMeshSetsNames

  VerbosityLevel CouplingSchemeBase::getVerbosityLevel() const noexcept {
    if (this->verbosity_level.has_value()) {
      return *(this->verbosity_level);
    }
    return ::mfem_mgis::getDefaultVerbosityLevel();
  }  // end of getVerbosityLevel

  void CouplingSchemeBase::setVerbosityLevel(const VerbosityLevel l) noexcept {
    this->verbosity_level = l;
  }  // end of setVerbosityLevel

  void CouplingSchemeBase::setLogStream(
      std::shared_ptr<std::ostream> l) noexcept {
    this->log_stream = l;
  }  // end of setLogStream

  std::shared_ptr<std::ostream>
  CouplingSchemeBase::getLogStreamPointer() noexcept {
    return this->log_stream;
  }  // end of getLogStreamPointer

  std::vector<const Provider *> CouplingSchemeBase::getProviders() noexcept {
    std::vector<const Provider *> providers;
    for (const auto &i : this->items) {
      if (auto p = std::dynamic_pointer_cast<Provider>(i); isValid(p)) {
        providers.push_back(p.get());
      }
      if (auto p = std::dynamic_pointer_cast<AbstractCouplingScheme>(i);
          isValid(p)) {
        const auto nproviders = p->getProviders();
        providers.insert(providers.end(), nproviders.begin(), nproviders.end());
      }
    }
    return providers;
  }  // end of getProviders

  //   bool CouplingSchemeBase::add(Context &ctx,
  //                                const Parameters &parameters) noexcept {
  //     if (parameters.size() != 1) {
  //       return ctx.registerErrorMessage(
  //           "expecting a dictionary with one element");
  //     }
  //     const auto &ci = *(parameters.begin());
  //     const auto &t = ci.first;
  //     // description of the coupling item
  //     const auto &ocid = get<Parameters>(ctx, ci.second);
  //     if (ocid->size() != 1) {
  //       return ctx.registerErrorMessage(
  //           "expecting a dictionary with one element to define a '" + t +
  //           "'");
  //     }
  //     const auto &d = *(ocid->begin());
  //     const auto &n = d.first;
  //     if (!is<Parameters>(d.second)) {
  //       return ctx.registerErrorMessage("expecting a dictionary to define the
  //       '" +
  //                                       t + "' name '" + n + "'");
  //     }
  //     const auto olparams = get<Parameters>(ctx, d.second);
  //     if (isInvalid(olparams)) {
  //       return {};
  //     }
  //     return this->addCouplingItem(ctx, t, n, *olparams);
  //   }  // end of add

  //   bool CouplingSchemeBase::addCouplingItem(
  //       Context &ctx,
  //       std::string_view t,
  //       std::string_view n,
  //       const Parameters &parameters) noexcept {
  //     if ((t == "Model") || (t == "model")) {
  //       return this->addModel(ctx, n, parameters);
  //     }
  //     return ctx.registerErrorMessage("unsupported type '" + std::string{t} +
  //                                     "'");
  //   }  // end of addCouplingItem

  bool CouplingSchemeBase::addCouplingItem(
      Context &ctx, std::shared_ptr<AbstractCouplingItem> i) noexcept {
    if (i.get() == nullptr) {
      return ctx.registerErrorMessage("invalid coupling item");
    }
    this->items.push_back(i);
    return true;
  }

  //   bool CouplingSchemeBase::addModel(Context &ctx,
  //                                     std::string_view n,
  //                                     const Parameters &parameters) noexcept
  //                                     {
  //     const auto &f = ModelFactory::get();
  //     return this->addModel(
  //         ctx, f.create(ctx, std::string{n}, this->physicalSystem_,
  //         parameters));
  //   }  // end of addModel

  bool CouplingSchemeBase::addModel(Context &ctx,
                                    std::shared_ptr<AbstractModel> m) noexcept {
    if (m.get() == nullptr) {
      return ctx.registerErrorMessage("invalid model");
    }
    return this->addCouplingItem(ctx, m);
  }  // end of addModel

  //   bool CouplingSchemeBase::declareDependencies(
  //       Context &ctx, DependenciesManager &dm) const noexcept {
  //     for (const auto &i : this->items) {
  //       auto cs = update(ctx, *i);
  //       const auto s = i->declareDependencies(ctx, dm);
  //       restore(ctx, cs);
  //       if (!s) {
  //         ctx.debug("* declareDependencies failed for '" +
  //                   getShortDescription(*i) + "'");
  //         return false;
  //       }
  //     }
  //     return true;
  //   }  // end of declareDependencies

  bool CouplingSchemeBase::performInitializationTaksAtTheBeginningOfTheTimeStep(
      Context &ctx, const TimeStep &ts) noexcept {
    for (const auto &i : this->items) {
      auto cs = update(ctx, *i);
      ctx.log(verboseLevel2,
              "* calling performInitializationTaksAtTheBeginningOfTheTimeStep "
              "for '" +
                  getShortDescription(*i) + "'");
      const auto s =
          i->performInitializationTaksAtTheBeginningOfTheTimeStep(ctx, ts);
      restore(ctx, cs);
      if (!s) {
        ctx.debug(
            "* performInitializationTaksAtTheBeginningOfTheTimeStep failed for "
            "'" +
            getShortDescription(*i) + "'");
        return false;
      }
    }
    return true;
  }  // end of performInitializationTaksAtTheBeginningOfTheTimeStep

  bool CouplingSchemeBase::executeInitialPostProcessingTasks(
      Context &ctx, const real t) noexcept {
    for (const auto &i : this->items) {
      auto cs = update(ctx, *i);
      ctx.log(verboseLevel2,
              "* calling executeInitialPostProcessingTasks for '" +
                  getShortDescription(*i) + "'");
      auto r = i->executeInitialPostProcessingTasks(ctx, t);
      restore(ctx, cs);
      if (!r) {
        ctx.debug("* executeInitialPostProcessingTasks failed for '" +
                  getShortDescription(*i) + "'");
        return false;
      }
    }
    return true;
  }  // end of executeInitialPostProcessingTasks

  std::optional<real> CouplingSchemeBase::getNextTimeIncrement(
      Context &ctx, const real t, const real te) const noexcept {
    auto dt = te - t;
    for (const auto &i : this->items) {
      auto cs = update(ctx, *i);
      const auto odt = i->getNextTimeIncrement(ctx, t, te);
      restore(ctx, cs);
      if (isInvalid(odt)) {
        ctx.debug("* computing the next time increment failed for '" +
                  getShortDescription(*i) + "'");
        return {};
      }
      dt = std::min(dt, *odt);
    }
    return dt;
  }  // end of getNextTimeIncrement

  bool CouplingSchemeBase::executePostProcessingTasks(Context &ctx,
                                                      const TimeStep & ts,
                                                      const bool b) noexcept {
    for (const auto &i : this->items) {
      auto cs = update(ctx, *i);
      ctx.log(verboseLevel2, "* calling executePostProcessingTasks for '" +
                                 getShortDescription(*i) + "'");
      auto r = i->executePostProcessingTasks(ctx, ts, b);
      restore(ctx, cs);
      if (!r) {
        ctx.debug("* executePostProcessingTasks failed for '" +
                  getShortDescription(*i) + "'");
        return false;
      }
    }
    return true;
  }  // end of executePostProcessingTasks

  bool CouplingSchemeBase::update(Context &ctx) noexcept {
    for (const auto &i : this->items) {
      auto cs = update(ctx, *i);
      const auto r = i->update(ctx);
      restore(ctx, cs);
      if (!r) {
        return false;
      }
    }
    return true;
  }  // end of update

  bool CouplingSchemeBase::revert(Context &ctx) noexcept {
    for (const auto &i : this->items) {
      auto cs = update(ctx, *i);
      const auto r = i->revert(ctx);
      restore(ctx, cs);
      if (!r) {
        return false;
      }
    }
    return true;
  }  // end of revert

  std::string CouplingSchemeBase::getCouplingItemsDescription() const noexcept {
    if (this->items.empty()) {
      return {};
    }
    auto d =
        std::string{"This coupling scheme contains the following items:\n"};
    for (const auto &i : this->items) {
      Context ctx;
      const auto od = i->describe(ctx, false, {});
      if (isInvalid(od)) {
        d += "\n- invalid description";
      } else {
        d += "\n- " + *od;
      }
    }
    return d;
  }  // end of getCouplingItemsDescription_

  CouplingSchemeBase::~CouplingSchemeBase() noexcept = default;

}  // namespace mfem_mgis
