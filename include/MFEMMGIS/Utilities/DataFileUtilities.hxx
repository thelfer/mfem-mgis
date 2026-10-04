/*!
 * \file   MFEMMGIS/Utilities/DataFileUtilities.hxx
 * \brief  This file declares various functions related to data files
 * \author Thomas Helfer
 * \date   27/12/2023
 */

#ifndef LIB_MFEMMGIS_UTILITIES_DATAFILEUTILITIES_HXX
#define LIB_MFEMMGIS_UTILITIES_DATAFILEUTILITIES_HXX 1

#include <iosfwd>
#include <string>
#include <vector>
#include <optional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Context;

  /*!
   * \brief list of supported data file formats
   */
  enum struct DataFileFormat {
    TXT,  //!< space separated values
    CSV   //!< comma separated values
  };

  /*!
   * \brief get the data file format from the extension of a file name
   * \return the data file format for the given file name, using the file
   * extension. Empty on failure.
   * \param[in, out] ctx: execution context
   * \param[in] f: file name
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<DataFileFormat>
  getDataFileFormatFromFileExtension(Context& ctx, std::string_view f) noexcept;
  /*!
   * \brief get the data file format from its name
   * \return the data file format from a string, empty on failure
   * \param[in, out] ctx: execution context
   * \param[in] f: name of the data file format, `txt` or `csv`
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<DataFileFormat>
  getDataFileFormat(Context& ctx, std::string_view f) noexcept;
  /*!
   * \brief get the value separator of a data file format
   * \return the value separator of the given data file format, empty on
   * failure
   * \param[in, out] ctx: execution context
   * \param[in] f: data file format
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::string_view>
  getValueSeparator(Context& ctx, const DataFileFormat f) noexcept;
  /*!
   * \brief write the header of the data file
   * \param[in, out] ctx: execution context
   * \param[in, out] os: output file stream
   * \param[in] f: data file format
   * \param[in] cnames: column names
   * \return true on success
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool writeDataFileHeader(
      Context& ctx,
      std::ostream& os,
      const DataFileFormat f,
      const std::vector<std::string>& cnames) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_UTILITIES_DATAFILEUTILITIES_HXX */
