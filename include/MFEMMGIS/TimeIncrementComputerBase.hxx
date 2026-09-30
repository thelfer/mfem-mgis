/*!
 * \file   MFEMMGIS/TimeIncrementComputerBase.hxx
 * \brief  This file declares the `TimeIncrementComputerBase` class
 * \date   04/12/2023
 */

#ifndef LIB_MFEM_MGIS_TIMEINCREMENTCOMPUTERBASE_HXX
#define LIB_MFEM_MGIS_TIMEINCREMENTCOMPUTERBASE_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractTimeIncrementComputer.hxx"

namespace mfem_mgis {

  //! \brief a common class for most time step computers
  struct MFEM_MGIS_EXPORT TimeIncrementComputerBase
      : AbstractTimeIncrementComputer {
    //! \brief constructor
    TimeIncrementComputerBase() noexcept;
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \return true
     */
    [[nodiscard]] bool initialize(Context& ctx) noexcept override;
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \return true
     */
    [[nodiscard]] bool prepareNextTimeStep(Context& ctx) noexcept override;
    //! \brief destructor
    ~TimeIncrementComputerBase() noexcept override;
  };

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_TIMEINCREMENTCOMPUTERBASE_HXX */
