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
  void getInformation(std::ostream&, const T&) noexcept {
  }  // end of getInformation

  template <typename T>
  void info(const T& t) noexcept {
    ::mfem_mgis::getInformation(getDefaultLogStream(), t);
  }  // end of info

  template <typename T>
  void info(std::ostream& os, const T& t) noexcept {
    ::mfem_mgis::getInformation(os, t);
  }  // end of info

  template <typename T>
  void info(Context& ctx, const T& t) noexcept {
    ::mfem_mgis::getInformation(ctx.log(), t);
  }  // end of info

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_INFO_IXX */
