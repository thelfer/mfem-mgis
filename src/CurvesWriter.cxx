/*!
 * \file   src/CurvesWriter.cxx
 * \brief  This file implements the `CurvesWriter` classs
 * \author Thomas Helfer
 * \date   07/09/2026
 */

#include "MFEMMGIS/MPI.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PhysicalSystem.hxx"
#include "MFEMMGIS/MeshDiscretization.hxx"
#include "MFEMMGIS/PostProcessing/CurvesWriter.hxx"

namespace mfem_mgis {

  std::map<std::string, std::string>
  CurvesWriter::getParametersDescription() noexcept {
    return {
        {"File", "name of the output file"},
        {"Precision",
         "number of significant digits in the output file. A precision of 0 "
         "means that only the integer part of the number is represented "
         "(see https://en.cppreference.com/w/cpp/io/manip/setprecision for "
         "details)"},
        {"FileFormat", "data file format. Expected values are 'txt' or 'csv'"}};
  }  // end of getParametersDescription

  CurvesWriter::CurvesWriter(Context &ctx,
                             const PhysicalSystem &ps,
                             const Parameters &parameters)
      : CurvesWriter(ctx, ps.getMeshDiscretization(), parameters) {
  }  // end of CurvesWriter

  CurvesWriter::CurvesWriter(Context &ctx,
                             const MeshDiscretization &m,
                             const Parameters &parameters)
      : isMainProcess(::mfem_mgis::isMainProcess(m)) {
    checkParameters(throwing, parameters,
                    CurvesWriter::getParametersDescription());
    if (!contains(parameters, "File")) {
      raise("no parameter 'File' specified");
    }
    if (!is<std::string>(throwing, parameters, "File")) {
      raise("the parameter 'File' must be a string");
    }
    const auto &fname = get<std::string>(throwing, parameters, "File");
    auto success = true;
    if (this->isMainProcess) {
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
        raise("invalid value for the 'precision' parameter (" +
              std::to_string(p) + ")");
      }
      if (this->isMainProcess) {
        this->out.precision(p);
      }
    }
    if (contains(parameters, "FileFormat")) {
      const auto &ff = get<std::string>(throwing, parameters, "FileFormat");
      const auto of = getDataFileFormat(ctx, ff);
      if (isInvalid(of)) {
        raise("invalid file format '" + ff + "'");
      }
      this->fileFormat = *of;
    } else {
      // try to deduce the file format from the file name
      const auto of = getDataFileFormatFromFileExtension(ctx, fname);
      if (isInvalid(of)) {
        raise("can't deduced file format");
      }
      this->fileFormat = *of;
    }
  }  // end of CurvesWriter

  bool CurvesWriter::addCurve(Context &ctx,
                              std::shared_ptr<const AbstractCurve> c) {
    if (!this->allowNewCurves) {
      return ctx.registerErrorMessage("no new curve allowed");
    }
    return this->curves.addCurve(ctx, c);
  }  // end of addCurve

  bool CurvesWriter::writeFileHeader(Context &ctx) {
    this->allowNewCurves = false;
    const auto cd = this->curves.getDescriptions();
    auto d = std::vector<std::string>{};
    d.push_back("time");
    d.insert(d.end(), cd.begin(), cd.end());
    return writeDataFileHeader(ctx, this->out, this->fileFormat, d);
  }  // end of writeFileHeader

  bool CurvesWriter::writeValues(Context &ctx,
                                 const TimeStep &ts,
                                 const TimeStepStage &tss) {
    this->allowNewCurves = false;
    const auto t = (tss == bts) ? ts.begin : ts.end;
    const auto os = getValueSeparator(ctx, this->fileFormat);
    if (isInvalid(os)) {
      return false;
    }
    const auto ovalues = this->curves.getValues(ctx, tss);
    if (isInvalid(ovalues)) {
      return false;
    }
    if (this->isMainProcess) {
      this->out << t;
      for (const auto &v : *ovalues) {
        this->out << *os << v;
      }
      this->out << '\n';
      this->out.flush();
    }
    return true;
  }  // end of writeValues

}  // namespace mfem_mgis
