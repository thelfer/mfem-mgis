/*!
 * \file   MFEMMGIS/DefaultTimeIncrementComputer.cxx
 * \brief  This file implements the `DefaultTimeIncrementComputer` class
 * \date   08/12/2023
 */

#include "MFEMMGIS/DefaultTimeIncrementComputer.hxx"

namespace mfem_mgis {

  DefaultTimeIncrementComputer::DefaultTimeIncrementComputer() noexcept =
      default;

  std::optional<real> DefaultTimeIncrementComputer::getNextTimeIncrement(
      Context &, const real t, const real se) const noexcept {
    return se - t;
  }  // end of getNextTimeStep

  DefaultTimeIncrementComputer::~DefaultTimeIncrementComputer() noexcept =
      default;

}  // end of namespace mfem_mgis
