/*!
 * \file   src/AbstractBoundaryCondition.cxx
 * \brief
 * \author Thomas Helfer
 * \date   27/09/2024
 */

#include <tuple>
#include "MFEMMGIS/AbstractBoundaryCondition.hxx"

namespace mfem_mgis {

  bool AbstractBoundaryCondition::setup(Context& ctx,
                                        const real t,
                                        const real dt) noexcept {
    try {
      this->setup(t, dt);
    } catch (...) {
      std::ignore = registerExceptionInErrorBacktrace(ctx);
      return false;
    }
    return true;
  }  // end of setup

  void AbstractBoundaryCondition::setup(const real, const real) {
  }  // end of setup

  AbstractBoundaryCondition::~AbstractBoundaryCondition() = default;

}  // end of namespace mfem_mgis
