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
   * \brief list of supported data file format
   */
  enum struct DataFileFormat {
    TXT,  // space separated values
    CSV   // comma separated values
  };

  /*!
   * \return the data file format for the given file name, using the file
   * extension \param[in] ctx: exectution context \param[in] f: file name
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<DataFileFormat>
  getDataFileFormatFromFileExtension(Context& ctx, std::string_view f) noexcept;
  /*!
   * \return the data file format from a string
   * \param[in] ctx: exectution context
   * \param[in] f: data file format
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<DataFileFormat>
  getDataFileFormat(Context& ctx, std::string_view f) noexcept;
  /*!
   * \return the data file format for the given file name, using the file
   * extention \param[in] ctx: exectution context \param[in] f: data file format
   */
  MFEM_MGIS_EXPORT [[nodiscard]] std::optional<std::string_view>
  getValueSeparator(Context& ctx, const DataFileFormat f) noexcept;
  /*!
   * \brief write the header of the data file
   * \param[in] ctx: exectution context
   * \param[in] os: output file stream
   * \param[in] f: data file format
   * \param[in] cnames: column names
   */
  MFEM_MGIS_EXPORT [[nodiscard]] bool writeDataFileHeader(
      Context& ctx,
      std::ostream& os,
      const DataFileFormat f,
      const std::vector<std::string>& cnames) noexcept;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_UTILITIES_DATAFILEUTILITIES_HXX */
