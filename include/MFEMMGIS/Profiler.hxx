/*!
 * \file   MFEMMGIS/Profiler.hxx
 * \brief  Profiling utilities and output manager.
 */

#ifndef LIB_MFEMMGIS_PROFILER_HXX
#define LIB_MFEMMGIS_PROFILER_HXX

#include <iostream>
#include <string>
#include <mfem.hpp>

namespace mgis {
  struct Context;
}

namespace mfem_mgis {
  namespace Profiler {
    namespace Utils {

      /*!
       * \brief Checks if the current MPI rank is the master.
       * \return True if master or running sequentially.
       */
      bool is_master();

      /*!
       * \brief Computes the maximum value across all MPI ranks.
       */
      double reduce_max(double a_duration);

      /*!
       * \brief Computes the sum across all MPI ranks.
       */
      double sum(double in);
      int sum(int in);
      int64_t sum(int64_t in);

      /*!
       * \brief Prints a single message strictly on the master rank.
       */
      template <typename Arg>
      void Message(Arg a_msg) {
        if (is_master()) {
          mfem::out << a_msg << std::endl;
        }
      }

      /*!
       * \brief Prints multiple arguments strictly on the master rank.
       */
      template <typename Arg, typename... Args>
      void Message(Arg a_msg, Args... a_msgs) {
        if (is_master()) {
          mfem::out << a_msg << " ";
          Message(a_msgs...);
        }
      }
    }  // namespace Utils

    namespace OutputManager {

      /*!
       * \brief Builds the default profiling output file name.
       */
      std::string build_name();

      /*!
       * \brief Prints the profiling table to the console.
       * \param[in] ctx The MGIS context.
       */
      void printTimeTable(const mgis::Context& ctx);

      /*!
       * \brief Writes profiling data to an auto-generated file.
       * \param[in] ctx The MGIS context.
       */
      void writeFile(const mgis::Context& ctx);

      /*!
       * \brief Writes profiling data to a specified file.
       * \param[in] ctx The MGIS context.
       * \param[in] a_name Output file name.
       */
      void writeFile(const mgis::Context& ctx, std::string a_name);

    }  // namespace OutputManager
  }  // namespace Profiler
}  // namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PROFILER_HXX */