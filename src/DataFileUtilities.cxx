/*!
 * \file   manta-core/io/data_file_utilities.cpp
 * \brief  This file implements various functions related to data files
 * \author Thomas Helfer
 * \date   27/12/2023
 */

#include <string>
#include <ostream>
#include "MFEMMGIS/Utilities/DataFileUtilities.hxx"

namespace mfem_mgis {

  std::optional<DataFileFormat> getDataFileFormat(Context &ctx,
                                                  std::string_view f) noexcept {
    if (f.empty()) {
      return ctx.registerErrorMessage("invalid file format");
    }
    if (f == "txt") {
      return DataFileFormat::TXT;
    }
    if (f == "csv") {
      return DataFileFormat::CSV;
    }
    return ctx.registerErrorMessage("unsupported file format '" +
                                    std::string{f} + "'");
  }

  std::optional<DataFileFormat> getDataFileFormatFromFileExtension(
      Context &ctx, std::string_view f) noexcept {
    if (f.empty()) {
      return ctx.registerErrorMessage("empty file name");
    }
    const auto p = f.find_last_of(".");
    if (p == std::string::npos) {
      return ctx.registerErrorMessage("no file extension");
    }
    return getDataFileFormat(ctx, f.substr(p + 1));
  }  // end of getDataFileFormatFromFileExtension

  std::optional<std::string_view> getValueSeparator(
      Context &ctx, const DataFileFormat f) noexcept {
    if (f == DataFileFormat::TXT) {
      return " ";
    }
    if (f == DataFileFormat::CSV) {
      return ",";
    }
    return ctx.registerErrorMessage("unsupported data file format");
  }  // end of getValueSeparator

  static std::string getColumnNumber(const size_type n) noexcept {
    if (n == 1) {
      return "first";
    }
    if (n == 2) {
      return "second";
    }
    if (n == 3) {
      return "third";
    }
    if (n == 4) {
      return "fourth";
    }
    if (n == 5) {
      return "fifth";
    }
    return std::to_string(n) + "th";
  }  // end of getColumnNumber

  static void writeTxtDataFileHeader(
      std::ostream &os, const std::vector<std::string> &cnames) noexcept {
    auto nc = size_type{1};
    for (const auto &n : cnames) {
      os << "# " << getColumnNumber(nc) << " column: " << n << '\n';
      ++nc;
    }
  }

  static void writeCsvDataFileHeader(
      std::ostream &os, const std::vector<std::string> &cnames) noexcept {
    auto p = cnames.begin();
    auto pe = cnames.end();
    while (p != pe) {
      os << '"' << *p << '"';
      if (++p != pe) {
        os << ",";
      }
    }
    os << '\n';
  }

  bool writeDataFileHeader(Context &ctx,
                           std::ostream &os,
                           const DataFileFormat f,
                           const std::vector<std::string> &cnames) noexcept {
    if (f == DataFileFormat::TXT) {
      writeTxtDataFileHeader(os, cnames);
    } else if (f == DataFileFormat::CSV) {
      writeCsvDataFileHeader(os, cnames);
    } else {
      return ctx.registerErrorMessage("unsupported data file format");
    }
    return true;
  }  // end of writeDataFileHandler

}  // namespace mfem_mgis