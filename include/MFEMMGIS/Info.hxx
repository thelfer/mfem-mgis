/*!
 * \file   MFEMMGIS/Info.hxx
 * \brief  A header declaring some useful functions to retrieve information
 * about various objects in MFEM/MGIS
 * \author Thomas Helfer
 * \date   24/02/2026
 */

#ifndef LIB_MFEM_MGIS_INFO_HXX
#define LIB_MFEM_MGIS_INFO_HXX

#include <iosfwd>
#include <type_traits>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief print information in the given log stream
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] os: output stream
   * \param[in] t: object for which information is requested
   * \return true on success
   *
   * \note by default, getInformation does nothing. This function is meant to be
   * specialized.
   */
  template <typename T>
  [[nodiscard]] bool getInformation(Context& ctx,
                                    std::ostream& os,
                                    const T& t) noexcept;

  /*!
   * \brief print information in the given stream
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] os: output stream
   * \param[in] t: object for which information is requested
   * \return true on success
   */
  template <typename T>
  bool info(Context& ctx, std::ostream& os, const T& t) noexcept;
  /*!
   * \brief print information in the log stream of the execution context
   *
   * \param[in, out] ctx: execution context
   * \param[in] t: object for which information is requested
   * \return true on success
   */
  template <typename T>
  bool info(Context& ctx, const T& t) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Info.ixx"

#endif /* LIB_MFEM_MGIS_INFO_HXX */
