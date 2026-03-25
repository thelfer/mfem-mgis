/*!
 * \file   MFEMMGIS/TimeIncrementComputerBase.cxx
 * \brief  This class implements the `TimeIncrementComputerBase` class
 * \date   04/12/2023
 */

#include "MFEMMGIS/TimeIncrementComputerBase.hxx"

namespace mfem_mgis {

  TimeIncrementComputerBase::TimeIncrementComputerBase() noexcept = default;

  bool TimeIncrementComputerBase::initialize(Context &) noexcept {
    return true;
  }  // end of initialize

  bool TimeIncrementComputerBase::prepareNextTimeStep(Context &) noexcept {
    return true;
  }  // end of prepareNextTimeStep

  TimeIncrementComputerBase::~TimeIncrementComputerBase() noexcept = default;

}  // end of namespace mfem_mgis
