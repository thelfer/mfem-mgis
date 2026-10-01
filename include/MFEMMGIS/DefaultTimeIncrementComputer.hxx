/*!
 * \file   MFEMMGIS/DefaultTimeIncrementComputer.hxx
 * \brief  This file declares the `DefaultTimeIncrementComputer` class
 * \date   08/12/2023
 */

#ifndef LIB_MFEM_MGIS_DEFAULTTIMEINCREMENTCOMPUTER_HXX
#define LIB_MFEM_MGIS_DEFAULTTIMEINCREMENTCOMPUTER_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/TimeIncrementComputerBase.hxx"

namespace mfem_mgis {

  //! \brief the default time step computer
  struct DefaultTimeIncrementComputer : TimeIncrementComputerBase {
    //! \brief constructor
    DefaultTimeIncrementComputer() noexcept;
    /*!
     * \return the time remaining until the end of the temporal sequence
     * \param[in, out] ctx: execution context
     * \param[in] t: current time in the temporal sequence
     * \param[in] te: end of the temporal sequence
     */
    std::optional<real> getNextTimeIncrement(
        Context& ctx, const real t, const real te) const noexcept override;
    //! \brief destructor
    ~DefaultTimeIncrementComputer() noexcept override;
  };  // end of DefaultTimeIncrementComputer

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_DEFAULTTIMEINCREMENTCOMPUTER_HXX */
