/*!
 * \file   manta/physical_system/physical_system.cpp
 * \brief  This file implements the `PhysicalSystem` class
 * \date   05/12/2022
 */

#include <sstream>
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractModel.hxx"
#include "MFEMMGIS/AbstractCouplingScheme.hxx"
#include "MFEMMGIS/LoopCouplingScheme.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"

namespace mfem_mgis {

  PhysicalSystem::PhysicalSystem(const MeshDiscretization &m) noexcept
      : mesh(m) {}  // end of PhysicalSystem

  std::optional<std::string> PhysicalSystem::describe(
      Context &ctx, const bool b, const Parameters &parameters) const noexcept {
    const auto options = std::map<std::string, std::string>{
        {"HelpOptions", "return the options of the describe method"},
        {"Mesh", "description of the mesh"},
        {"CouplingScheme", "description of the coupling scheme"},
        {"Loadings", "list of the loadings"},
        {"LoadingEvolutions", "list of the loading evolutions"},
        {"PostProcessings", "list of post-processings"},
        {"Header", "boolean stating if headers shall be displayed"},
        {"HeaderLevel", "level of headers (integer)"}};
    if (!checkParameters(ctx, parameters, options)) {
      return {};
    }
    const auto oho = get_if<bool>(ctx, parameters, "helpOptions", false);
    if (isInvalid(oho)) {
      return {};
    }
    if (*oho) {
      auto d = std::string{};
      for (auto first = true; const auto &[k, o] : options) {
        if (k == "HelpOptions") {
          continue;
        }
        if (!first) {
          d += '\n';
        }
        d += "- " + k + ": " + o;
        first = false;
      }
      return d;
    }
    //     auto [omesh, meshParameters] =
    //         [&ctx, b, &parameters]() -> std::pair<std::optional<bool>,
    //         Parameters> {
    //       if (contains(parameters, "mesh")) {
    //         if (is<Parameters>(parameters, "mesh")) {
    //           return {true, get<Parameters>(parameters, "mesh")};
    //         }
    //       }
    //       return {get_if<bool>(ctx, parameters, "mesh", b), Parameters{}};
    //     }();
    const auto ocsd = get_if<bool>(ctx, parameters, "CouplingScheme", b);
    //     const auto oloadings = get_if<bool>(ctx, parameters, "loadings", b);
    //     const auto oloadingEvolutions =
    //         get_if<bool>(ctx, parameters, "loadingsEvolutions", b);
    //     const auto opostProcessings =
    //         get_if<bool>(ctx, parameters, "postProcessings", b);
    const auto oheader = get_if<bool>(ctx, parameters, "Header", b);
    const auto oheaderLevel =
        get_if<size_type>(ctx, parameters, "HeaderLevel", 1);
    if (areInvalid(/* omesh, */ ocsd, /* oloadings, oloadingEvolutions,
                                         opostProcessings, */
                   oheader, oheaderLevel)) {
      return {};
    }
    if (*oheaderLevel < 1) {
      return ctx.registerErrorMessage("invalid header level ('" +
                                      std::to_string(*oheaderLevel) + "')");
    }
    const auto nsections =
        /* *omesh + */ *ocsd
        /* + *oloadings + *oloadingEvolutions + *opostProcessings */;
    auto d = std::string{};
    auto add = [&d](std::string_view o) {
      if (o.empty()) {
        return;
      }
      if (!d.empty()) {
        d += "\n\n";
      }
      d += o;
    };
    auto header = [nsections, &add, b = *oheader,
                   lvl = *oheaderLevel](const std::string &h) {
      if (b) {
        if (nsections != 0) {
          add(std::string(lvl, '#') + " " + h);
        }
      }
    };
    //     // mesh
    //     if (*omesh) {
    //       meshParameters["headerLevel"] = *oheaderLevel + 1;
    //       const auto omd = this->getMesh().describe(ctx, true,
    //       meshParameters); if (isInvalid(omd)) {
    //         return {};
    //       }
    //       header("Description of the mesh");
    //       if (omd->empty()) {
    //         add("No mesh description available.");
    //       } else {
    //         add(*omd);
    //       }
    //     }
    // coupling scheme
    if (*ocsd) {
      header("Description of the coupling scheme");
      if (this->coupling_scheme.get() != nullptr) {
        auto csd = "The physical system is based on the '" +
                   this->coupling_scheme->getName() + "' coupling scheme.";
        const auto ocsd2 = this->coupling_scheme->describe(
            ctx, true, {{"ShortDescription", false}});
        if (isInvalid(ocsd2)) {
          return {};
        }
        add(csd);
        add(*ocsd2);
      } else {
        add("No coupling scheme defined.");
      }
    }
    //     // loadings
    //     if (*oloadings) {
    //       header("Description of the loadings");
    //       if (this->loadings_.empty()) {
    //         add("No loadings defined.");
    //       } else {
    //         add("The following loadings are defined:\n");
    //         for (const auto &l : this->loadings_) {
    //           d += "\n- ";
    //           d += l->getShortDescription();
    //         }
    //       }
    //     }
    //     // loading evolutions
    //     if (*oloadingEvolutions) {
    //       header("Description of the loading evolutions");
    //       const auto levs =
    //           this->loadingEvolutions_->getRegisteredLoadingEvolutions();
    //       if (levs.empty()) {
    //         add("No loading evolution defined.");
    //       } else {
    //         add("The following loading evolutions are defined:\n");
    //         for (const auto &ln : levs) {
    //           const auto ol = this->loadingEvolutions_->get(ctx, ln);
    //           if (isInvalid(ol)) {
    //             return {};
    //           }
    //           d += "\n- ";
    //           d += ln;
    //           d += " of type '" + (*ol)->getName() + "'";
    //         }
    //       }
    //     }
    //     // post-processings
    //     if (*opostProcessings) {
    //       header("Description of the post-processings");
    //       if (this->post_processings.empty()) {
    //         add("No post-processing defined.");
    //       } else {
    //         add("The following post-processings are defined:\n");
    //         for (const auto &p : this->post_processings) {
    //           d += "\n- ";
    //           d += p->getName();
    //         }
    //       }
    //     }
    return d;
  }  // end of describe

  MeshDiscretization PhysicalSystem::getMeshDiscretization() const noexcept {
    return this->mesh;
  }  // end of getMeshDiscretization

  bool PhysicalSystem::isCouplingSchemeDefined() const noexcept {
    return this->coupling_scheme.get() != nullptr;
  }  // end of isCouplingSchemeDefined

  //   bool PhysicalSystem::setCouplingScheme(Context &ctx,
  //                                          std::string_view n,
  //                                          const Parameters &p) noexcept {
  //     if (this->isCouplingSchemeDefined()) {
  //       return ctx.registerErrorMessage("coupling scheme already defined");
  //     }
  //     return this->setCouplingScheme(ctx, buildCouplingScheme(ctx, *this, n,
  //     p));
  //   }  // end of setCouplingScheme

  bool PhysicalSystem::setCouplingScheme(
      Context &ctx, std::shared_ptr<AbstractCouplingScheme> c) noexcept {
    if (this->isCouplingSchemeDefined()) {
      return ctx.registerErrorMessage("coupling scheme already defined");
    }
    if (c.get() == nullptr) {
      return ctx.registerErrorMessage("invalid coupling scheme specified");
    }
    if (c->getMeshDiscretization() != this->mesh) {
      return ctx.registerErrorMessage("inconsistent meshes");
    }
    this->coupling_scheme = std::move(c);
    return true;
  }  // end of setCouplingScheme

  //   bool PhysicalSystem::setModel(Context &ctx,
  //                                 std::string_view n,
  //                                 const Parameters &p) noexcept {
  //     using namespace ::manta::internal;
  //     if (this->isCouplingSchemeDefined()) {
  //       return ctx.registerErrorMessage("coupling scheme already defined");
  //     }
  //     return this->setModel(ctx, buildModel(ctx, *this, n, p));
  //   }  // end of setModel

  bool PhysicalSystem::setModel(Context &ctx,
                                std::shared_ptr<AbstractModel> m) noexcept {
    if (this->isCouplingSchemeDefined()) {
      return ctx.registerErrorMessage("coupling scheme already defined");
    }
    if (m.get() == nullptr) {
      return ctx.registerErrorMessage("invalid model specified");
    }
    if (m->getMeshDiscretization() != this->mesh) {
      return ctx.registerErrorMessage("inconsistent meshes");
    }
    auto c = make_shared<LoopCouplingScheme>(ctx, this->mesh);
    if (isInvalid(c)) {
      return false;
    }
    if (!c->addModel(ctx, m)) {
      return false;
    }
    return this->setCouplingScheme(ctx, c);
  }  // end of setModel

  bool PhysicalSystem::updateLoadingsAtTheBeginningOfTheTimeStep(
      Context &, const TimeStep &) noexcept {
    //     for (auto &l : this->loadings_) {
    //       if (!l->updateLoading(ctx)) {
    //         return false;
    //       }
    //     }
    return true;
  }  // end of updateLoadingsAtTheBeginningOfTheTimeStep

  bool PhysicalSystem::performInitializationTaksAtTheBeginningOfTheTimeStep(
      Context &ctx, const TimeStep &ts) noexcept {
    if (isInvalid(this->coupling_scheme)) {
      return ctx.registerErrorMessage("no coupling scheme defined");
    }
    ctx.debug(
        "PhysicalSystem::performInitializationTaksAtTheBeginningOfTheTimeStep");
    return this->coupling_scheme
        ->performInitializationTaksAtTheBeginningOfTheTimeStep(ctx, ts);
  }  // end of performInitializationTaksAtTheBeginningOfTheTimeStep

  std::optional<real> PhysicalSystem::getNextTimeIncrement(
      Context &ctx, const real t, const real te) const noexcept {
    if (isInvalid(this->coupling_scheme)) {
      return ctx.registerErrorMessage("no coupling scheme defined");
    }
    return this->coupling_scheme->getNextTimeIncrement(ctx, t, te);
  }  // end of getNextTimeIncrement

  std::pair<ExitStatus, std::optional<ComputeNextStateOutput>>
  PhysicalSystem::computeNextState(Context &ctx, const TimeStep &ts) noexcept {
    if (isInvalid(this->coupling_scheme)) {
      std::ignore = ctx.registerErrorMessage("no coupling scheme defined");
      return {ExitStatus::unrecoverableError, {}};
    }
    return this->coupling_scheme->computeNextState(ctx, ts);
  }  // end of computeNextState

  bool PhysicalSystem::executeInitialPostProcessingTasks(
      Context &ctx, const real t) noexcept {
    if (isInvalid(this->coupling_scheme)) {
      return ctx.registerErrorMessage("no coupling scheme defined");
    }
#pragma message("HERE")
    //     for (const auto &p : this->post_processings) {
    //       ctx.log(
    //           verboseLevel2,
    //           "* calling executeInitialPostProcessingTasks on post-processing
    //           '" +
    //               p->getName() + "'");
    //       auto r = p->executeInitialPostProcessingTasks(ctx);
    //       if (!r) {
    //         ctx.debug(
    //             "* executeInitialPostProcessingTasks failed for
    //             post-processing '" + p->getName() + "'");
    //         return false;
    //       }
    //     }
    return this->coupling_scheme->executeInitialPostProcessingTasks(ctx, t);
  }  // end of executeInitialPostProcessingTasks

  bool PhysicalSystem::executePostProcessingTasks(Context &ctx,
                                                  const TimeStep &ts,
                                                  const bool b) {
    if (isInvalid(this->coupling_scheme)) {
      return ctx.registerErrorMessage("no coupling scheme defined");
    }
#pragma message("HERE")
    //     for (const auto &p : this->post_processings) {
    //       ctx.log(verboseLevel2,
    //               "* calling executePostProcessingTasks on post-processing '"
    //               +
    //                   p->getName() + "'");
    //       auto r = p->executePostProcessingTasks(ctx, b);
    //       if (!r) {
    //         ctx.debug("* executePostProcessingTasks failed for
    //         post-processing '" +
    //                   p->getName() + "'");
    //         return false;
    //       }
    //     }
    return this->coupling_scheme->executePostProcessingTasks(ctx, ts, b);
  }  // end of executePostProcessingTasks

  bool PhysicalSystem::update(Context &ctx) noexcept {
    if (isInvalid(this->coupling_scheme)) {
      return ctx.registerErrorMessage("no coupling scheme defined");
    }
    return this->coupling_scheme->update(ctx);
  }  // end of update

  bool PhysicalSystem::revert(Context &ctx) noexcept {
    if (isInvalid(this->coupling_scheme)) {
      return ctx.registerErrorMessage("no coupling scheme defined");
    }
    return this->coupling_scheme->revert(ctx);
  }  // end of revert

  bool PhysicalSystem::addPostProcessing(
      Context &ctx, std::string_view n, const Parameters &parameters) noexcept {
#pragma message("HERE")
    //     const auto &f = PostProcessingFactory::get();
    //     return this->addPostProcessing(
    //         ctx, f.create(ctx, std::string{n}, *this, parameters));
    return false;
  }  // end of addPostProcessing

  bool PhysicalSystem::addPostProcessing(
      Context &ctx, std::shared_ptr<AbstractPostProcessing> p) noexcept {
#pragma message("HERE")
    //     if (p.get() == nullptr) {
    //       return ctx.registerErrorMessage("invalid post-processing");
    //     }
    //     for (const auto &lp : this->post_processings) {
    //       if (lp.get() == p.get()) {
    //         // post-processing already registered
    //         return true;
    //       }
    //     }
    //     this->post_processings.push_back(p);
    return true;
  }  // end of addPostProcessing

  PhysicalSystem::~PhysicalSystem() noexcept = default;

}  // namespace mfem_mgis
