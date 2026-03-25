/*!
 * \file   src/AbstractCouplingItem.cxx
 * \brief  This file declares the `AbstractCouplingItem` class
 * \date   06/12/2022
 */

#include <fstream>
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/AbstractCouplingItem.hxx"

namespace mfem_mgis {

  const std::string AbstractCouplingItem::verbosityLevelParameter =
      "VerbosityLevel";
  const std::string AbstractCouplingItem::logFileParameter = "LogFile";

  AbstractCouplingItem::~AbstractCouplingItem() = default;

  std::string getShortDescription(const AbstractCouplingItem &i) noexcept {
    auto id = i.getName();

    const auto locations = i.getLocations();
    if (!locations.empty()) {
      id += " on [";
      auto pe = locations.end();
      for (auto p = locations.begin(); p != pe;) {
        id += *p;
        if (++p != pe) {
          id += ", ";
        }
      }
      id += "]";
    }
    return id;
  }  // end of getShortDescription

  std::map<std::string, std::string>
  getCouplingItemParametersDescription() noexcept {
    return {{AbstractCouplingItem::verbosityLevelParameter, "verbosity level"},
            {AbstractCouplingItem::logFileParameter, "log file"}};
  }  // end of getCouplingItemParametersDescription

  bool handleCouplingItemParameters(Context &ctx,
                                    AbstractCouplingItem &i,
                                    const Parameters &params) {
    if (contains(params, AbstractCouplingItem::verbosityLevelParameter)) {
      const auto ol = get<std::string>(
          ctx, params, AbstractCouplingItem::verbosityLevelParameter);
      if (isInvalid(ol)) {
        return false;
      }
      const auto olvl = convertToVerbosityLevel(ctx, *ol);
      if (isInvalid(olvl)) {
        return false;
      }
      i.setVerbosityLevel(*olvl);
    }
    if (contains(params, AbstractCouplingItem::logFileParameter)) {
      const auto of =
          get<std::string>(ctx, params, AbstractCouplingItem::logFileParameter);
      if (isInvalid(of)) {
        return false;
      }
      auto os = std::make_shared<std::ofstream>(*of);
      if (!(*os)) {
        return ctx.registerErrorMessage("unable to open file '" + *of + "'");
      }
      i.setLogStream(os);
    }
    return true;
  }  // end of handleCouplingItemParameters

}  // end of namespace mfem_mgis
