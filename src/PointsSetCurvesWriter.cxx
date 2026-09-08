/*!
 * \file   src/PointsSetCurvesWriter.cxx
 * \brief  This file implements the `PointsSetCurvesWriter` classs
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/PostProcessing/PointsSetCurvesWriter.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  PointsSetCurvesWriter::getParametersDescription() noexcept {
    auto d = PointsSetCurves::getParametersDescription();
    d.insert(
        {{"File", "name of the output file"},
         {"Precision",
          "number of significant digits in the output file. A precision of 0 "
          "means that only the integer part of the number is represented "
          "(see https://en.cppreference.com/w/cpp/io/manip/setprecision for "
          "details)"}});
    return d;
  }  // end of getParametersDescription

  PointsSetCurvesWriter::PointsSetCurvesWriter(const PhysicalSystem &ps,
                                               const Parameters &parameters)
      : PointsSetCurvesWriter(ps.getMeshDiscretization(), parameters) {
  }  // end of PointsSetCurvesWriter

  PointsSetCurvesWriter::PointsSetCurvesWriter(const MeshDiscretization &m,
                                               const Parameters &parameters)
      : mesh(m),
        curves(m,
               extract(throwing,
                       parameters,
                       PointsSetCurves::getParametersDescription())) {
    checkParameters(throwing, parameters,
                    PointsSetCurvesWriter::getParametersDescription());
    if (!contains(parameters, "File")) {
      raise("no parameter 'file' specified");
    }
    if (!is<std::string>(throwing, parameters, "File")) {
      raise("the parameter 'file' must be a string");
    }
    const auto &fname = get<std::string>(throwing, parameters, "File");
    auto success = true;
    if (isMainProcess(this->mesh)) {
      this->out.open(fname);
      // force the C locale. Otherwise, the LC_NUMERICS environment
      // variable is taken into account and that could lead to
      // inconsistent outputs between users. More over, the french
      // locale uses the comma as a separator, which is incompatible
      // with the `csv` output
      this->out.imbue(std::locale("C"));
      success = this->out.is_open();
    }
    if (!isTrueOnAllProcesses(m, success)) {
      raise("opening of file '" + fname + "' failed");
    }
    //
    if (contains(parameters, "Precision")) {
      const auto p = get<int>(throwing, parameters, "Precision");
      if (p < 0) {
        raise("invalid value for the 'Precision' parameter (" +
              std::to_string(p) + ")");
      }
      if (isMainProcess(this->mesh)) {
        this->out.precision(p);
      }
    }
  }  // end of PointsSetCurvesWriter

#ifdef MFEM_USE_MPI

  bool PointsSetCurvesWriter::add(Context &ctx,
                                  std::string_view n,
                                  const GridFunction<true> &f) noexcept {
    if (!this->allowNewPointsSetCurves) {
      return ctx.registerErrorMessage("no new points set curves allowed");
    }
    return this->curves.add(ctx, n, f);
  }  // end of add

#endif /* MFEM_USE_MPI */

  bool PointsSetCurvesWriter::add(Context &ctx,
                                  std::string_view n,
                                  const GridFunction<false> &f) noexcept {
    if (!this->allowNewPointsSetCurves) {
      return ctx.registerErrorMessage("no new points set curves allowed");
    }
    return this->curves.add(ctx, n, f);
  }  // end of add

  bool PointsSetCurvesWriter::writeFileHeader(Context &ctx) {
    this->allowNewPointsSetCurves = false;
    auto success = true;
    if (isMainProcess(this->mesh)) {
      try {
        auto nl = size_type{1};
        if (this->curves.exportCurvilinearAbscissa()) {
          this->out << "# first line: values of the curvilinear abscissa\n";
          ++nl;
        }
        if (this->curves.exportCoordinates()) {
          const auto d = this->curves.getSpaceDimension();
          for (size_type i = 0; i != d; ++i) {
            this->out << "# " << nl << "th line: values of " << i + 1
                      << " coordinate\n";
            ++nl;
          }
        }
        for (const auto &d : this->curves.getValuesDescription()) {
          this->out << "# " << nl << "th line: " << d << "\n";
          ++nl;
        }
      } catch (...) {
        std::ignore = registerExceptionInErrorBacktrace(ctx);
        success = false;
      }
    }
    if (!isTrueOnAllProcesses(this->mesh, success)) {
      return false;
    }
    return true;
  }  // end of writeFileHeader

  bool PointsSetCurvesWriter::writeValues(Context &ctx,
                                          const TimeStep &ts,
                                          const TimeStepStage &tss) {
    this->allowNewPointsSetCurves = false;
    const auto t = (tss == bts) ? ts.begin : ts.end;
    if (isMainProcess(this->mesh)) {
      this->out << "\n#Time " << t << '\n';
    }
    auto write_values = [this](const std::vector<real> &values) {
      if (isMainProcess(this->mesh)) {
        bool first = true;
        for (const auto &v : values) {
          if (!first) {
            this->out << " ";
          }
          this->out << v;
          first = false;
        }
        this->out << '\n';
      }
    };
    if (this->curves.exportCurvilinearAbscissa()) {
      auto success = true;
      auto ovalues = this->curves.getCurvilinearAbscissa(ctx);
      if (isInvalid(ovalues)) {
        success = false;
      } else {
        write_values(*ovalues);
      }
      if (!isTrueOnAllProcesses(this->mesh, success)) {
        return false;
      }
    }
    if (this->curves.exportCoordinates()) {
      auto success = true;
      auto ovalues = this->curves.getCoordinates(ctx);
      if (isInvalid(ovalues)) {
        success = false;
      } else {
        for (const auto &row : *ovalues) {
          write_values(row);
        }
      }
      if (!isTrueOnAllProcesses(this->mesh, success)) {
        return false;
      }
    }
    auto success = true;
    auto ovalues = this->curves.getValues(ctx, tss);
    if (isInvalid(ovalues)) {
      success = false;
    } else {
      for (const auto &row : *ovalues) {
        write_values(row);
      }
    }
    if (!isTrueOnAllProcesses(this->mesh, success)) {
      return false;
    }
    if (isMainProcess(this->mesh)) {
      this->out.flush();
    }
    return true;
  }  // end of writeValues

}  // namespace mfem_mgis
