/*!
 * \file   include/MFEMMGIS/BoundaryUtilities.hxx
 * \brief  This file declares functions describing the elements and the degrees
 * of freedom on a boundary
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
   * \brief describe a boundary by its faces
   * \return a vector of pairs associating for each boundary element its index
   * and the index of the adjacent element.
   * \tparam parallel: boolean stating if the computation is done in parallel.
   * \param[in] p: non linear evolution problem
   * \param[in] bid: boundary identifier
   */
  template <bool parallel>
  std::vector<std::pair<size_type, size_type>> buildFacesDescription(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const size_type bid);

  /*!
   * \brief list the elements having degrees of freedom on a boundary
   * \return a vector of pairs associating the index of each element having
   * at least one degree of freedom on the boundary with the local indexes of
   * these degrees of freedom, grouped by component.
   * \tparam parallel: boolean stating if the computation is done in parallel.
   * \param[in] p: non linear evolution problem
   * \param[in] bid: boundary identifier
   * \note all the elements having a degree of freedom on the boundary are
   * needed, even those touching it only by an edge or a vertex: the inner
   * force at a node is the sum of the contributions of all the elements
   * sharing this node.
   * \note in parallel, MFEM gives a boundary element only to the process
   * owning the adjacent element. A process may thus own an element touching
   * the boundary without owning any boundary element. The degrees of freedom
   * of the boundary are therefore marked by `GetEssentialVDofs`, which
   * synchronizes the markers between the processes sharing them, as MFEM does
   * for essential boundary conditions. This function is thus collective: all
   * processes must call it.
   * \note the local indexes follow `GetElementVDofs`, which groups the degrees
   * of freedom of an element by component whatever the ordering of the space
   * (`byNODES` or `byVDIM`).
   */
  template <bool parallel>
  std::vector<std::pair<size_type,                  //< element number
                        std::vector<                //< storage per components
                            std::vector<size_type>  //< local index of
                                                    // the degrees of
                                                    //  freedom
                            >>>
  getElementsDegreesOfFreedomOnBoundary(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const size_type bid);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/BoundaryUtilities.ixx"

#endif /* LIB_MFEM_MGIS_BOUNDARYUTILITIES_HXX */
