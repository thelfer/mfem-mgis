/*!
 * \file   MFEMMGIS/Info.ixx
 * \brief
 * \author Thomas Helfer
 * \date   24/02/2026
 */

#ifndef LIB_MFEM_MGIS_INFO_IXX
#define LIB_MFEM_MGIS_INFO_IXX

namespace mfem_mgis {

  template <typename T>
  bool getInformation(Context&, std::ostream&, const T&) noexcept {
    return true;
  }  // end of getInformation

  template <typename T>
  bool info(Context& ctx, std::ostream& os, const T& t) noexcept {
    return ::mfem_mgis::getInformation(ctx, os, t);
  }  // end of info

  template <typename T>
  bool info(Context& ctx, const T& t) noexcept {
    return ::mfem_mgis::getInformation(ctx, ctx.log(), t);
  }  // end of info

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_INFO_IXX */
