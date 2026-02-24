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
   * \param[out] os: output stream
   * \param[in, out] t: object for which information are requested
   *
   * \note by default, getInformation does nothing. This function is meant to be
   * overloaded.
   */
  template <typename T>
  void getInformation(std::ostream&, const T&) noexcept;

  /*!
   * \brief print information in the default log stream
   *
   * \param[in, out] t: object for which information are requested
   */
  template <typename T>
  void info(const T&) noexcept;
  /*!
   * \brief print information in the given stream
   *
   * \param[out] os: output stream
   * \param[in, out] t: object for which information are requested
   */
  template <typename T>
  void info(std::ostream&, const T&) noexcept;
  /*!
   * \brief print information in the log stream associated with the context
   *
   * \param[in, out] ctx: execution context
   * \param[in, out] t: object for which information are requested
   */
  template <typename T>
  void info(Context&, const T&) noexcept;

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Info.ixx"

#endif /* LIB_MFEM_MGIS_INFO_HXX */
