#include <MFEMMGIS/Profiler.hxx>
#ifdef MFEM_USE_MPI
#include "mpi.h"
#endif
#include <numeric>

namespace mfem_mgis {
  namespace Profiler {
    namespace Utils {
      constexpr int master = 0;

      double sum(double in) {
        double res = 0;
#ifdef MFEM_USE_MPI
        MPI_Reduce(&in, &res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        return res;
      }

      int sum(int in) {
        int res = 0;
#ifdef MFEM_USE_MPI
        MPI_Reduce(&in, &res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        return res;
      }

      int64_t sum(int64_t in) {
        int64_t res = 0;
#ifdef MFEM_USE_MPI
        MPI_Reduce(&in, &res, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
        return res;
      };

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
        double global = 0.0;
        MPI_Reduce(&a_duration, &global, 1, MPI_DOUBLE, MPI_MAX, master,
                   MPI_COMM_WORLD);  // master rank is 0
        return global;
#else
        return a_duration;
#endif
      }

    };  // namespace Utils

    namespace timers {
      using duration = std::chrono::duration<double>;

      // constructor
      ProfilerTimeSection::ProfilerTimeSection()
          : m_daughter()  // only used for root
      {
        m_name = "root";
        m_iteration = 1;
        m_level = 0;
        m_mother = nullptr;
      }

      ProfilerTimeSection::ProfilerTimeSection(std::string name,
                                               ProfilerTimeSection* mother)
          : m_daughter(), m_duration(0) {
        m_name = name;
        m_iteration = 1;
        m_level = mother->m_level + 1;
        m_mother = mother;
      }

      ProfilerTimeSection* ProfilerTimeSection::find(std::string name) {
        assert(this != nullptr);
        for (auto it = m_daughter.begin(); it < m_daughter.end(); it++) {
          if ((*it)->m_name == name) {
            (*it)->m_iteration++;
            return (*it);
          }
        }
        ProfilerTimeSection* myTmp = new ProfilerTimeSection(name, this);
        m_daughter.push_back(myTmp);
        return myTmp;
      }

      void ProfilerTimeSection::printReplicate(int begin,
                                               int end,
                                               std::string motif) {
        for (int i = begin; i < end; i++) mfem::out << motif;
      }

      void ProfilerTimeSection::space() { mfem::out << " "; }

      void ProfilerTimeSection::column() { mfem::out << "|"; }

      void ProfilerTimeSection::endline() { mfem::out << std::endl; }

      void ProfilerTimeSection::printBanner(int shift) {
        if (m_name == "root") {
#ifndef MFEM_USE_MPI
          Profiler::Utils::Message(
              " MPI feature is disable for timers, if you use MPI please add "
              "-DMFEM_USE_MPI ");
#else
          if (Profiler::Utils::is_master()) {
            Profiler::Utils::Message(" MPI feature activated, rank 0:");
#endif
          mfem::out << std::endl;
          mfem::out << "Glossary: " << std::endl;
          mfem::out << "NS: Newton Solver Class" << std::endl;
          mfem::out << "NLEPIB: Non Linear Evolution Problem Implementation "
                       "Base Class"
                    << std::endl;
          mfem::out << "FED: Finite Element Discretization Class" << std::endl;
          mfem::out << std::endl;

          std::string start_name = " |-- start timetable ";
          mfem::out << start_name;
          int end = shift + nColumns * (cWidth + 1) + 1;
          printReplicate(start_name.size(), end, "-");
          column();
          endline();
          std::string name = " |    name";
          mfem::out << name;
          printReplicate(name.size(), shift + 1, " ");
          for (int i = 0; i < nColumns; i++) {
            column();
            int size = cName[i].size();
            printReplicate(0, (int(cWidth) - size - 1), " ");
            mfem::out << cName[i];
            space();
          }
          column();
          endline();
          space();
          column();
          printReplicate(2, end, "-");
          column();
          endline();
#ifdef MFEM_USE_MPI
        }
#endif
      }
    }  // namespace timers

    void ProfilerTimeSection::printEnding(int shift) {
      if (m_name == "root") {
        if (Profiler::Utils::is_master()) {
          shift += nColumns * (cWidth + 1) + 1;  // +1 for "|";
          std::string end_name = " |-- end timetable ";
          mfem::out << end_name;
          printReplicate(end_name.size(), shift, "-");
          column();
          endline();
        }
      }
    }

    duration* ProfilerTimeSection::get_ptr_duration() { return &m_duration; }

    void ProfilerTimeSection::print(int shift, double total_time) {
      assert(total_time >= 0);
      std::string cValue[nColumns];
      if (Profiler::Utils::is_master()) {
        int realShift = shift;
        space();
        column();
        space();
        int currentShift = 3;
        for (int i = 0; i < int(m_level) - 1; i++) {
          int spaceSize = 3;
          for (int j = 0; j < spaceSize; j++) space();
          currentShift += spaceSize;
        }
        if (m_level > 0) {
          mfem::out << "|--";
          currentShift += 3;
        }
        mfem::out << "> " << m_name;
        currentShift += m_name.size() + 1;
        printReplicate(currentShift, realShift, " ");

        cValue[0] = std::to_string(m_iteration);
      }
#ifdef MFEM_USE_MPI
      double local = m_duration.count();
      int size = -1;
      MPI_Comm_size(MPI_COMM_WORLD, &size);

      assert(size > 0);
      std::vector<double> list;

      if (Profiler::Utils::is_master()) list.resize(size);

      MPI_Gather(&local, 1, MPI_DOUBLE, list.data(), 1, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);  // master rank is 0

      if (Profiler::Utils::is_master()) {
        const auto [min, max] = std::minmax_element(list.begin(), list.end());
        auto global_max = *max;
        auto global_min = *min;
        auto sum = std::accumulate(list.begin(), list.end(), double(0.));
        auto global_mean = sum / double(size);
        auto part_time = (global_max / total_time) * 100;

        assert(global_mean >= 0);
        assert(global_min >= 0);
        assert(global_max >= 0);

        assert(global_max >= global_mean);
        assert(global_mean >= global_min);

        cValue[1] = std::to_string(global_min);
        cValue[2] = std::to_string(global_mean);
        cValue[3] = std::to_string(global_max);
        cValue[4] = std::to_string(part_time) + "%";
        cValue[5] = std::to_string((global_max / global_mean) - 1) + "%";
      }
#else
          cValue[1] = std::to_string(m_duration.count());
          cValue[2] = std::to_string((m_duration.count() / total_time) * 100);
#endif
      if (Profiler::Utils::is_master()) {
        for (int i = 0; i < nColumns; i++) {
          column();
          int _size = cValue[i].size();
          printReplicate(0, (int(cWidth) - _size - 1), " ");
          mfem::out << cValue[i];
          space();
        }
        column();
        endline();
      }
    }

    std::string ProfilerTimeSection::getName() { return m_name; }

    double ProfilerTimeSection::get_duration() { return m_duration.count(); }

    int ProfilerTimeSection::get_iteration() { return m_iteration; }

    int ProfilerTimeSection::get_level() { return m_level; }

    std::vector<ProfilerTimeSection*>& ProfilerTimeSection::get_daughter() {
      return m_daughter;
    }

    ProfilerTimeSection* ProfilerTimeSection::get_mother() { return m_mother; }

    void init_timers() {
      Profiler::Utils::Message(" Init timers ");
      ProfilerTimeSection*& root_timer_ptr =
          Profiler::timers::get_timer<ROOT>();
      assert(root_timer_ptr == nullptr);
      root_timer_ptr = new ProfilerTimeSection();
      ProfilerTimeSection*& current = Profiler::timers::get_timer<CURRENT>();
      current = root_timer_ptr;
      assert(current != nullptr);
      Profiler::timer::start_global_timer<ROOT>();
    }

    void print_and_write_timers() {
      Profiler::timer::end_global_timer<ROOT>();
      Profiler::OutputManager::writeFile();
      Profiler::OutputManager::printTimeTable();
    }

    void print_timers() {
      Profiler::timer::end_global_timer<ROOT>();
      Profiler::OutputManager::printTimeTable();
    }

    void write_timers() {
      Profiler::timer::end_global_timer<ROOT>();
      Profiler::OutputManager::writeFile();
    }
  };  // namespace Profiler

  namespace timer {
    TimeSection::TimeSection(duration* acc) {
      m_duration = acc;
      start();
    }

    inline void TimeSection::start() { m_start = steady_clock::now(); }

    inline void TimeSection::end() {
      assert(m_duration != nullptr && "duration has to be initialised");
      m_stop = steady_clock::now();
      *m_duration += m_stop - m_start;
      assert(m_duration->count() >= 0);
    }

    TimeSection::~TimeSection() {
      end();
      auto& current_timer = Profiler::timers::get_timer<CURRENT>();
      current_timer = current_timer->get_mother();
    }
  };  // namespace timer

  namespace OutputManager {
    std::string build_name() {
      std::string base_name = "mfem-mgis";
#ifdef MFEM_USE_MPI
      int mpiSize;
      MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
      std::string file_name =
          base_name + "." + std::to_string(mpiSize) + ".perf";
#else
          int nthreads = 0;
#if defined(_OPENMP)
#pragma omp parallel
          { nthreads = omp_get_num_threads(); }
#endif
          std::string file_name =
              base_name + "." + std::to_string(nthreads) + ".perf";
#endif
      return file_name;
    }

    void printTimeTable() {
      ProfilerTimeSection* root_timer = Profiler::timers::get_timer<ROOT>();
      double runtime = root_timer->get_duration();
      runtime =
          Profiler::Utils::reduce_max(runtime);  // if MPI, else return runtime

      auto my_print = [](ProfilerTimeSection* a_ptr, int a_shift,
                         double a_runtime) {
        a_ptr->print(a_shift, a_runtime);
      };

//#define SortPrintTimeTable 0
#ifdef SortPrintTimeTable
      auto sort_comp = [](ProfilerTimeSection* a_ptr,
                          ProfilerTimeSection* b_ptr) {
        return a_ptr->get_duration() > b_ptr->get_duration();
      };
#endif

      auto max_length = [](ProfilerTimeSection* a_ptr, int& a_count,
                           int& a_nbElem) {
        int length = a_ptr->get_level() * 3 + a_ptr->getName().size();
        a_count = std::max(a_count, length);
        a_nbElem++;
      };
      int count(0), nbElem(0);

      recursive_call(max_length, root_timer, count, nbElem);
      count += 6;
      root_timer->printBanner(count);
#ifdef SortPrintTimeTable
      recursive_sorted_call(my_print, sort_comp, root_timer, count, runtime);
#else
          recursive_call(my_print, root_timer, count, runtime);
#endif
      root_timer->printEnding(count);
    }

    void writeFile() {
      std::string name = build_name();
      writeFile(name);
    }

    void writeFile(std::string a_name) {
      using namespace Profiler::Utils;
      // using Profiler::Utils::reduce_max;

      std::ofstream myFile(a_name, std::ofstream::out);
      ProfilerTimeSection* root_timer = Profiler::timers::get_timer<ROOT>();
      auto rootTime = root_timer->get_duration();
      rootTime = Profiler::Utils::reduce_max(rootTime);
      auto my_write = [rootTime](ProfilerTimeSection* a_ptr,
                                 std::ofstream& a_file) {
        std::string space;
        std::string motif = "   ";

        for (int i = 0; i < a_ptr->get_level(); i++) space += motif;

        const auto max_time = reduce_max(a_ptr->get_duration());

        if (is_master()) {
          a_file << space << a_ptr->getName() << " " << a_ptr->get_iteration()
                 << " " << max_time << " " << (max_time / rootTime) * 100
                 << std::endl;
        }
      };

      recursive_call(my_write, root_timer, myFile);
    }
  }  // namespace OutputManager
};   // namespace mfem_mgis
}
;  // namespace Profiler
