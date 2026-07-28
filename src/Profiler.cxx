/*!
 * \file   src/Profiler.cxx
 * \brief  This file implements the profiling output utilities
 */

#include "MGIS/Context.hxx"
#include "MGIS/Profiling.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <mfem.hpp>

#ifdef MFEM_USE_MPI
#include "mpi.h"
#endif

namespace mfem_mgis {
  namespace Profiler {
    namespace Utils {
      constexpr int master = 0;
      constexpr int cWidth = 20;

      int get_mpi_size() {
#ifdef MFEM_USE_MPI
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return size;
#else
        return 1;
#endif
      }

      double sum(double in) {
        double res = in;
#ifdef MFEM_USE_MPI
        if (get_mpi_size() > 1) {
            MPI_Reduce(&in, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
#endif
        return res;
      }

      int sum(int in) {
        int res = in;
#ifdef MFEM_USE_MPI
        if (get_mpi_size() > 1) {
            MPI_Reduce(&in, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
#endif
        return res;
      }

      int64_t sum(int64_t in) {
        int64_t res = in;
#ifdef MFEM_USE_MPI
        if (get_mpi_size() > 1) {
            MPI_Reduce(&in, &res, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        }
#endif
        return res;
      }

      bool is_master() {
#ifdef MFEM_USE_MPI
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return (rank == master);
#else
        return true;
#endif
      }

      double reduce_max(double a_duration) {
#ifdef MFEM_USE_MPI
        double global = a_duration;
        if (get_mpi_size() > 1) {
            MPI_Reduce(&a_duration, &global, 1, MPI_DOUBLE, MPI_MAX, master, MPI_COMM_WORLD);
        }
        return global;
#else
        return a_duration;
#endif
      }

    }  // namespace Utils

    namespace OutputManager {

      std::string build_name() {
        std::string base_name = "mfem-mgis";
        int mpiSize = Utils::get_mpi_size();
        
        if (mpiSize > 1) {
            return base_name + "." + std::to_string(mpiSize) + ".perf";
        } else {
            int nthreads = 0;
#if defined(_OPENMP)
#pragma omp parallel
            { nthreads = omp_get_num_threads(); }
#endif
            return base_name + "." + std::to_string(nthreads) + ".perf";
        }
      }

      static void printReplicate(int begin, int end, std::string motif) {
        for (int i = begin; i < end; i++) mfem::out << motif;
      }

      static int get_max_length(const mgis::ProfilingData& node, int level) {
        int length = level * 3 + node.name.size();
        int max_len = length;
        for (const auto& child : node.children) {
          max_len = std::max(max_len, get_max_length(*child, level + 1));
        }
        return max_len;
      }

      static void printBanner(int shift, int mpi_size) {
        if (!Utils::is_master()) return;
        
        mfem::out << std::endl;
        mfem::out << "Glossary: " << std::endl;
        mfem::out << "NS: Newton Solver Class" << std::endl;
        mfem::out << "NLEPIB: Non Linear Evolution Problem Implementation Base Class" << std::endl;
        mfem::out << "FED: Finite Element Discretization Class" << std::endl;
        mfem::out << std::endl;

        std::vector<std::string> headers;
        if (mpi_size > 1) {
            headers = {"number Of Calls", "min(s)", "mean(s)", "max(s)", "part(%)", "imb(%)"};
        } else {
            headers = {"number Of Calls", "time(s)", "part(%)"};
        }

        int nCols = headers.size();
        std::string start_name = " |-- start timetable ";
        mfem::out << start_name;
        int end = shift + nCols * (Utils::cWidth + 1) + 1;
        printReplicate(start_name.size(), end, "-");
        mfem::out << "|\n";
        
        std::string name = " |    name";
        mfem::out << name;
        printReplicate(name.size(), shift + 1, " ");
        for (int i = 0; i < nCols; i++) {
          mfem::out << "|";
          int size = headers[i].size();
          printReplicate(0, Utils::cWidth - size - 1, " ");
          mfem::out << headers[i] << " ";
        }
        mfem::out << "|\n |";
        printReplicate(2, end, "-");
        mfem::out << "|\n";
      }

      static void printNode(const mgis::ProfilingData& node, int level, double root_time, int shift, int mpi_size) {
        int nCols = (mpi_size > 1) ? 6 : 3;
        std::vector<std::string> cValue(nCols);
        
        if (Utils::is_master()) {
          mfem::out << " | ";
          int currentShift = 3;
          for (int i = 0; i < level - 1; i++) {
            mfem::out << "   ";
            currentShift += 3;
          }
          if (level > 0) {
            mfem::out << "|--> ";
            currentShift += 5;
          } else {
            mfem::out << "> ";
            currentShift += 2;
          }
          mfem::out << node.name;
          currentShift += node.name.size();
          printReplicate(currentShift, shift, " ");
          cValue[0] = std::to_string(node.calls);
        }

        double local_time = (level == 0) ? root_time : node.time_in_seconds;

        if (mpi_size > 1) {
#ifdef MFEM_USE_MPI
          std::vector<double> list;
          if (Utils::is_master()) list.resize(mpi_size);
          
          MPI_Gather(&local_time, 1, MPI_DOUBLE, list.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

          if (Utils::is_master()) {
            const auto [min, max] = std::minmax_element(list.begin(), list.end());
            auto global_max = *max;
            auto global_min = *min;
            auto sum = std::accumulate(list.begin(), list.end(), 0.0);
            auto global_mean = sum / double(mpi_size);
            auto part_time = (global_max / root_time) * 100.0;

            auto fmt = [](double val) {
              std::ostringstream out;
              out << std::fixed << std::setprecision(6) << val;
              return out.str();
            };

            cValue[1] = fmt(global_min);
            cValue[2] = fmt(global_mean);
            cValue[3] = fmt(global_max);
            cValue[4] = fmt(part_time) + "%";
            cValue[5] = (global_mean > 0) ? fmt(((global_max / global_mean) - 1.0) * 100.0) + "%" : "0.000000%";
          }
#endif
        } else {
          if (Utils::is_master()) {
            auto fmt = [](double val) {
              std::ostringstream out;
              out << std::fixed << std::setprecision(6) << val;
              return out.str();
            };
            cValue[1] = fmt(local_time);
            cValue[2] = fmt((local_time / root_time) * 100.0) + "%";
          }
        }

        if (Utils::is_master()) {
          for (int i = 0; i < nCols; i++) {
            mfem::out << "|";
            int _size = cValue[i].size();
            printReplicate(0, Utils::cWidth - _size - 1, " ");
            mfem::out << cValue[i] << " ";
          }
          mfem::out << "|\n";
        }

        for (const auto& child : node.children) {
          printNode(*child, level + 1, root_time, shift, mpi_size);
        }
      }

      void printTimeTable(const mgis::Context& ctx) {
        if (!ctx.isProfilingEnabled()) return;

        int mpi_size = Utils::get_mpi_size();
        const auto& root = ctx.getProfilingResultTree();
        double root_time = Utils::reduce_max(root.time_in_seconds);
        
        if (root_time <= 0.0) {
          for (const auto& child : root.children) {
            root_time += Utils::reduce_max(child->time_in_seconds);
          }
        }

        if (root_time <= 0.0) root_time = 1e-9; 

        int shift = get_max_length(root, 0) + 6;

        printBanner(shift, mpi_size);
        printNode(root, 0, root_time, shift, mpi_size);
        
        if (Utils::is_master()) {
          int nCols = (mpi_size > 1) ? 6 : 3;
          int end = shift + nCols * (Utils::cWidth + 1) + 1;
          std::string end_name = " |-- end timetable ";
          mfem::out << end_name;
          printReplicate(end_name.size(), end, "-");
          mfem::out << "|\n";
        }
      }

      static void writeNode(const mgis::ProfilingData& node, int level, double root_time, std::ofstream& file) {
        std::string space;
        std::string motif = "   ";

        for (int i = 0; i < level; i++) space += motif;

        const auto max_time = (level == 0) ? root_time : Utils::reduce_max(node.time_in_seconds);

        if (Utils::is_master()) {
          file << space << node.name << " " << node.calls
               << " " << max_time << " " << (max_time / root_time) * 100
               << std::endl;
        }

        for (const auto& child : node.children) {
          writeNode(*child, level + 1, root_time, file);
        }
      }

      void writeFile(const mgis::Context& ctx) {
        writeFile(ctx, build_name());
      }

      void writeFile(const mgis::Context& ctx, std::string a_name) {
        if (!ctx.isProfilingEnabled()) return;

        std::ofstream myFile(a_name, std::ofstream::out);
        const auto& root = ctx.getProfilingResultTree();
        auto rootTime = Utils::reduce_max(root.time_in_seconds);
        
        if (rootTime <= 0.0) {
          for (const auto& child : root.children) {
            rootTime += Utils::reduce_max(child->time_in_seconds);
          }
        }

        if (rootTime <= 0.0) rootTime = 1e-9;

        writeNode(root, 0, rootTime, myFile);
      }

    }  // namespace OutputManager
  }  // namespace Profiler
}  // namespace mfem_mgis