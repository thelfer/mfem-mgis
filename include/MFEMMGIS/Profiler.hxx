/*!
 * \file   MFEMMGIS/Profiler.hxx
 * \brief  This file declares the profiling utilities and output manager
 */

#ifndef LIB_MFEMMGIS_PROFILER_HXX
#define LIB_MFEMMGIS_PROFILER_HXX

#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include <mfem.hpp>
#include "MGIS/Context.hxx"

namespace mfem_mgis {
  namespace Profiler {
    namespace Utils {
      bool is_master();
      double reduce_max(double a_duration);

      double sum(double in);
      int sum(int in);
      int64_t sum(int64_t in);

      template <typename Arg>
      void Message(Arg a_msg) {
        if (is_master()) {
          mfem::out << a_msg << std::endl;
        }
      }

      template <typename Arg, typename... Args>
      void Message(Arg a_msg, Args... a_msgs) {
        if (is_master()) {
          mfem::out << a_msg << " ";
          Message(a_msgs...);
        }
      }
    }  // namespace Utils

    namespace OutputManager {
      std::string build_name();
      void printTimeTable(const mgis::Context& ctx);
      void writeFile(const mgis::Context& ctx);
      void writeFile(const mgis::Context& ctx, std::string a_name);
    }  // namespace OutputManager

  }  // namespace Profiler
}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PROFILER_HXX */