/*!
 * \file   include/MFEMMGIS/BoundaryUtilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   28/03/2021
 */

#ifndef LIB_MFEM_MGIS_BOUNDARYUTILITIES_HXX
#define LIB_MFEM_MGIS_BOUNDARYUTILITIES_HXX

#include <vector>
#include <utility>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  /*!
   * \return a description of the boundary by a vector of pair
   * associating for each face its identifier and the identifier of the
   * adjacent element.
   * \tparam parallel: boolean stating if the computation is done in parallel.
   * \param[in] p: non linear evolution problem
   * \param[in] bid: boundary identifier
   */
  template <bool parallel>
  std::vector<std::pair<size_type, size_type>> buildFacesDescription(
      NonLinearEvolutionProblemImplementation<parallel>&, const size_type);

  /*!
   * \brief return a structure which associates the global number of the
   * selected elements to the local indexes of its degrees of freedom sorted by
   * components.
   * \tparam parallel: boolean stating if the computation is done in parallel.
   * \param[in] p: non linear evolution problem
   * \param[in] bid: boundary identifier
   */
  template <bool parallel>
  std::vector<std::pair<size_type,                  //< element number
                        std::vector<                //< storage per components
                            std::vector<size_type>  //< local index of
                                                    // the degree of
                                                    //  freedoms
                            >>>
  getElementsDegreesOfFreedomOnBoundary(
      NonLinearEvolutionProblemImplementation<parallel>&, const size_type);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/BoundaryUtilities.ixx"

#endif /* LIB_MFEM_MGIS_BOUNDARYUTILITIES_HXX */
