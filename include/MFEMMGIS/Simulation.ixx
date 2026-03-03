/*!
 * \file   MFEMMGIS/Simulation.ixx
 * \brief
 * \author Thomas Helfer
 * \date   03/03/2026
 */

#ifndef LIB_MFEM_MGIS_SIMULATION_IXX
#define LIB_MFEM_MGIS_SIMULATION_IXX

namespace mfem_mgis {

  inline ExitStatus::ExitStatus() noexcept
      : status(Status{.value = Status::success}) {}  // end of ExitStatus

  inline ExitStatus::ExitStatus(const Status s) noexcept
      : status(s) {}  // end of ExitStatus

  inline ExitStatus::ExitStatus(const bool s) noexcept
      : status(s ? Status{.value = Status::success}
                 : Status{.value = Status::unrecoverableError}) {
  }  // end of ExitStatus

  inline void ExitStatus::setStatus(const bool s) noexcept {
    this->operator=(ExitStatus(s));
  }

  inline bool ExitStatus::shallContinue() const noexcept {
    return (this->status == ExitStatus::success) ||
           (this->status == ExitStatus::unreliableResults);
  }  // end of shallContinue

  inline bool ExitStatus::shallStop() const noexcept {
    return (this->status == ExitStatus::recoverableError) ||
           (this->status == ExitStatus::unrecoverableError);
  }  // end of shallStop

  inline bool ExitStatus::update(const ExitStatus &o) noexcept {
    this->status = this->status > o.status ? this->status : o.status;
    return this->shallContinue();
  }  // end of update

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_SIMULATION_IXX */
