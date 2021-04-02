/*!
 * \file   Profiler.hxx
 * \brief
 * \author Thomas Helfer
 * \date   02/04/2021
 */

#ifndef LIB_MFEM_MGIS_PROFILER_HXX
#define LIB_MFEM_MGIS_PROFILER_HXX

#include <map>
#include <memory>
#include <atomic>
#include <string>
#include <iosfwd>
#include <cstdint>
#include <string_view>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  /*!
   * \brief structure used to profile the execution
   */
  struct MFEM_MGIS_EXPORT  Profiler {
    // foward declaration
    struct TimeSection;
    /*!
     * \brief class used to measure the time spent in a given portion of the
     * code.
     */
    struct Timer {
      //! \brief constructor
      Timer(TimeSection&);
      //! \brief move constructor
      Timer(Timer&&);
      //! \brief destructor
      ~Timer();

     private:
      //! \brief time section
      TimeSection& ts;
      //! \brief start
      timespec start;
      //! \brief end
      timespec end;
    }; // end of struct Timer
    /*!
     * \brief structure handling the measurement of a portion of the code
     */
    struct TimeSection {
      /*!
       * \brief constructor of a time section
       * \param[in] p: pointer to the time section
       */
      TimeSection(TimeSection* const = nullptr);
      //! \brief move constructor
      TimeSection(TimeSection&&) = default;
      //! \brief print the results
      void print(std::ostream&, const std::string&, const std::string&) const;
      //! \brief close the section
      void close(const uint64_t);
      /*!
       * \return a new timer for the given subsection
       * \param[in] n: name of the subsection
       */
      Timer getTimer(std::string_view);
      //! \brief paren section
      TimeSection* const parent;
      //! \brief total time spend in this region
      std::atomic<uint64_t> measure;
      //! \brief child sections
      std::map<std::string, std::unique_ptr<TimeSection>, std::less<>>
          subsections;
    };  // end of struct TimeSection

    //! \return the unique instance of this class
    static Profiler& getProfiler();
    //! \brief display the results
    void print(std::ostream&) const;
    /*!
     * \brief start measuring the time spend in a given portion of the code.
     * \param[in] name of the section
     *
     * \note the measurement
     */
    Timer getTimer(std::string_view);

   private:

    //! \brief default constructor
    Profiler();
    // explicitly deleted constructors and assignement operators
    Profiler(Profiler&&) = delete;
    Profiler(const Profiler&) = delete;
    Profiler& operator=(Profiler&&) = delete;
    Profiler& operator=(const Profiler&) = delete;
    //! \brief top level section
    TimeSection main;
    //! \brief current section
    TimeSection* current;
  };  // end of struct Profiler

  /*!
   * \brief return a timer for the given code section
   * \param[in] n: name
   */
  MFEM_MGIS_EXPORT Profiler::Timer getTimer(std::string_view);

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PROFILER_HXX */
