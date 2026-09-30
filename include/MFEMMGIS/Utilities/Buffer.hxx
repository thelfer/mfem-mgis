/*!
 * \file   MFEMMGIS/Buffer.hxx
 * \brief
 * \author Thomas Helfer
 * \date   30/04/2025
 */

#ifndef LIB_MFEMMGIS_BUFFER_HXX
#define LIB_MFEMMGIS_BUFFER_HXX

#include <span>
#include <array>
#include <vector>
#include <type_traits>
#include <MFEMMGIS/Config.hxx>

namespace mfem_mgis {

  /*!
   * \brief buffer of real values.
   * A `std::vector` if the extent is dynamic, a `std::array` otherwise.
   */
  template <size_type Extent = dynamic_extent>
  using Buffer =
      std::conditional_t<Extent == dynamic_extent,
                         std::vector<real>,
                         std::array<real, static_cast<std::size_t>(Extent)>>;

  /*!
   * \brief create a span on a buffer
   * \param[in] b: buffer
   * \return a span on the given buffer
   */
  template <std::size_t Extent>
  auto makeSpan(const std::array<real, Extent>& b) noexcept {
    if constexpr (Extent == std::dynamic_extent) {
      return std::span<const real>(b);
    } else {
      return std::span<const real, Extent>(b);
    }
  }  // end of makeSpan

  /*!
   * \brief create a span on a buffer
   * \param[in] b: buffer
   * \return a span on the given buffer
   */
  inline auto makeSpan(const std::vector<real>& b) noexcept {
    return std::span<const real>(b);
  }  // end of makeSpan

}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_BUFFER_HXX */
