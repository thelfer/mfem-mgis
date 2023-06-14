#pragma once

#include <iostream>
#include <chrono>
#include <vector>
#include <omp.h>
#include <fstream>
#include <algorithm>
#include <cassert>

#include <mfem.hpp>
// variables
#ifdef MFEM_USE_MPI
#include <mpi.h>
const size_t cWidth = 20;
const size_t nColumns = 6;
const std::string cName[nColumns] = {
    "number Of Calls", "min(s)",  "mean(s)",
    "max(s)",          "part(%)", "imb(%)"};  // [1-Imax/Imean]%
#else
const size_t cWidth = 20;
const size_t nColumns = 3;
const std::string cName[nColumns] = {"number Of Calls", "max(s),", "part(%)"};
#endif

enum enumTimeSection { CURRENT, ROOT };

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
    };  // namespace Utils

    namespace timers {
      class ProfilerTimeSection {
        using duration = std::chrono::duration<double>;

       public:
        ProfilerTimeSection();
        ProfilerTimeSection(std::string name, ProfilerTimeSection* mother);
        ProfilerTimeSection* find(std::string name);

        // printer functions
        //
        void printReplicate(size_t begin, size_t end, std::string motif);
        void space();
        void column();
        void endline();
        void printBanner(size_t shift);
        void printEnding(size_t shift);
        duration* get_ptr_duration();
        void print(size_t shift, double runtime);

        // accessor
        //
        std::string getName();
        std::size_t get_iteration();
        std::size_t get_level();
        std::vector<ProfilerTimeSection*>& get_daughter();
        ProfilerTimeSection* get_mother();
        double get_duration();

       private:
        std::string m_name;
        std::size_t m_iteration;
        std::size_t m_level;
        std::vector<ProfilerTimeSection*> m_daughter;
        ProfilerTimeSection* m_mother;
        duration m_duration;
      };

      void init_timers();
      void print_and_write_timers();
      void print_timers();
      void write_timers();

      template <enumTimeSection T>
      ProfilerTimeSection*& get_timer() {
        static ProfilerTimeSection* __current;
        return __current;
      }

      template <typename Lambda>
      double chrono_section(Lambda&& lambda) {
        using steady_clock = std::chrono::steady_clock;
        using time_point = std::chrono::time_point<steady_clock>;
        time_point tic, toc;
        tic = steady_clock::now();
        lambda();
        toc = steady_clock::now();
        auto measure = toc - tic;
        return measure.count();
      }
    };  // namespace timers

    namespace timer {
      using duration = std::chrono::duration<double>;
      using steady_clock = std::chrono::steady_clock;
      using time_point = std::chrono::time_point<steady_clock>;
      class TimeSection {
       public:
        TimeSection(duration* acc);
        void start();
        void end();
        ~TimeSection();

       private:
        time_point m_start;
        time_point m_stop;
        duration* m_duration;
      };

      template <enumTimeSection T>
      TimeSection*& get_timer() {
        static TimeSection* __timer;
        return __timer;
      }

      template <enumTimeSection T>
      void start_global_timer() {
        assert(T == enumTimeSection::ROOT);
        auto& timer = get_timer<T>();
        timer = new TimeSection(
            Profiler::timers::get_timer<T>()->get_ptr_duration());
        timer->start();  // reset start
      }

      template <enumTimeSection T>
      void end_global_timer() {
        assert(T == enumTimeSection::ROOT);
        auto timer = get_timer<T>();
        timer->end();
      }
    }  // namespace timer

    namespace OutputManager {
      using Profiler::timers::ProfilerTimeSection;

      std::string build_name();
      void printTimeTable();
      void writeFile();
      void writeFile(std::string a_name);

      template <typename Func, typename... Args>
      void recursive_call(Func& func, ProfilerTimeSection* ptr, Args&... arg) {
        func(ptr, arg...);
        auto& daughters = ptr->get_daughter();
        for (auto& it : daughters) recursive_call(func, it, arg...);
      }

      template <typename Func, typename Sort, typename... Args>
      void recursive_sorted_call(Func& func,
                                 Sort mySort,
                                 ProfilerTimeSection* ptr,
                                 Args&... arg) {
        func(ptr, arg...);
        auto& daughters = ptr->get_daughter();
        std::sort(daughters.begin(), daughters.end(), mySort);
        for (auto& it : daughters)
          recursive_sorted_call(func, mySort, it, arg...);
      }

    };  // namespace OutputManager
  };    // namespace Profiler

};  // namespace mfem_mgis

#define CatchTimeSection(X)                                                    \
  mfem_mgis::Profiler::timers::ProfilerTimeSection*& current =                            \
      mfem_mgis::Profiler::timers::get_timer<CURRENT>();                                  \
  assert(current != nullptr && "do not use an undefined ProfilerTimeSection"); \
  current = current->find(X);                                                  \
  mfem_mgis::Profiler::timer::TimeSection non_generic_name(current->get_ptr_duration());

#define CatchNestedTimeSection(X)                                              \
  current = mfem_mgis::Profiler::timers::get_timer<CURRENT>();                            \
  assert(current != nullptr && "do not use an undefined ProfilerTimeSection"); \
  current = current->find(X);                                                  \
  mfem_mgis::Profiler::timer::TimeSection non_nested_generic_name(                        \
      current->get_ptr_duration());
