/*!
 * \file   src/Profiler.cxx
 * \brief
 * \author Thomas Helfer
 * \date   02/04/2021
 */

#include <ctime>
#include <vector>
#include <numeric>
#include <ostream>
#include <algorithm>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Profiler.hxx"
#ifdef MFEM_USE_MPI
#include "mpi.h"
#endif /* MFEM_USE_MPI */

namespace mfem_mgis {

  /*!
   * \brief add a new measure
   * start : start of the measure
   * end   : end of the measure
   */
  static inline uint64_t get_measure(const timespec& start,
                                     const timespec& end) {
    /* http://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime
     */
    timespec temp;
    if ((end.tv_nsec - start.tv_nsec) < 0) {
      temp.tv_sec = end.tv_sec - start.tv_sec - 1;
      temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    } else {
      temp.tv_sec = end.tv_sec - start.tv_sec;
      temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
    return 1000000000 * temp.tv_sec + temp.tv_nsec;
  }  // end of add_measure

  /*!
   * print a time to the specified stream
   */
  static void print_time(std::ostream& os, const uint64_t time) {
    constexpr uint64_t musec_d = 1000;
    constexpr uint64_t msec_d = 1000 * musec_d;
    constexpr uint64_t sec_d = msec_d * 1000;
    constexpr uint64_t min_d = sec_d * 60;
    constexpr uint64_t hour_d = min_d * 60;
    constexpr uint64_t days_d = hour_d * 24;
    uint64_t t = time;
    const uint64_t ndays = t / days_d;
    t -= ndays * days_d;
    const uint64_t nhours = t / hour_d;
    t -= nhours * hour_d;
    const uint64_t nmins = t / min_d;
    t -= nmins * min_d;
    const uint64_t nsecs = t / sec_d;
    t -= nsecs * sec_d;
    const uint64_t nmsecs = t / msec_d;
    t -= nmsecs * msec_d;
    const uint64_t nmusecs = t / musec_d;
    t -= nmusecs * musec_d;
    if (ndays > 0) {
      os << ndays << "days ";
    }
    if (nhours > 0) {
      os << nhours << "hours ";
    }
    if (nmins > 0) {
      os << nmins << "mins ";
    }
    if (nsecs > 0) {
      os << nsecs << "secs ";
    }
    if (nmsecs > 0) {
      os << nmsecs << "msecs ";
    }
    if (nmusecs > 0) {
      os << nmusecs << "musecs ";
    }
    os << t << "nsecs";
  }  // end pf print

  Profiler& Profiler::getProfiler() {
    static Profiler p;
    return p;
  }  // end of getProfiler

  Profiler::Timer::Timer(Profiler::TimeSection& s) : ts(s) {
    ::clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(this->start));
  }  // end of Timer

  Profiler::Timer::~Timer() {
    ::clock_gettime(CLOCK_THREAD_CPUTIME_ID, &(this->end));
    this->ts.close(get_measure(this->start, this->end));
  }  // end of Timer

  Profiler::TimeSection::TimeSection(Profiler::TimeSection* p)
      : parent(p), measure(0) {}  // end of TimeSection

  void Profiler::TimeSection::close(const uint64_t m) {
    auto& p = Profiler::getProfiler();
    p.current = this->parent;
    this->measure += m;
  }  // end of TimeSection::close

  Profiler::Timer Profiler::TimeSection::getTimer(std::string_view n) {
    auto& p = Profiler::getProfiler();
    auto ps = this->subsections.find(n);
    if (ps == this->subsections.end()) {
      const auto k = std::string(n);
      ps = this->subsections.insert({k, std::make_unique<TimeSection>(this)})
               .first;
    }
    p.current = ps->second.get();
    return Timer(*(ps->second.get()));
  }  // end of getTimer

  void Profiler::TimeSection::print(std::ostream& os,
                                    const std::string& n,
                                    const std::string& s) const {
#ifdef MFEM_USE_MPI
    int rank;
    int gsize;
    std::vector<std::uint64_t> measures;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      MPI_Comm_size(MPI_COMM_WORLD, &gsize);
      measures.resize(gsize);
    }
    MPI_Gather(&(this->measure), 1, MPI_UINT64_T, measures.data(), 1,
               MPI_UINT64_T, 0, MPI_COMM_WORLD);
    if (rank == 0) {
      const auto r = std::minmax_element(measures.begin(), measures.end());
      const auto min = *(r.first);
      const auto max = *(r.second);
      const auto mean_value =
          std::accumulate(measures.begin(), measures.end(), std::uint64_t{}) /
          gsize;
      os << s << "- " << n << ": ";
      print_time(os, min);
      os << " ";
      print_time(os, max);
      os << " ";
      print_time(os, mean_value);
      os << '\n';
    }
#else  /* MFEM_USE_MPI */
    os << s << "- " << n << ": ";
    print_time(os, this->measure);
    os << '\n';
#endif /* MFEM_USE_MPI */
    for (const auto& ts : this->subsections) {
      ts.second->print(os, ts.first, s + "  ");
    }
  }  // end of print

  Profiler::Profiler() : main(), current(&main) {}  // end of Profiler

  Profiler::Timer Profiler::getTimer(std::string_view n) {
    return this->current->getTimer(n);
  }  // end of getTimer

  void Profiler::print(std::ostream& os) const {
    for (const auto& ts : this->main.subsections) {
      ts.second->print(os, ts.first, "");
    }
  }  // end of print

  Profiler::Timer getTimer(std::string_view n) {
    return Profiler::getProfiler().getTimer(n);
  }  // end of getTimer

}  // end of namespace mfem_mgis