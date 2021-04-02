/*!
 * \file   include/MFEMMGIS/MFEMForward.hxx
 * \brief
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_MFEM_FORWARD_HXX
#define LIB_MFEM_MGIS_MFEM_FORWARD_HXX


#include <type_traits>
#include "mfem/config/config.hpp"

#ifndef MFEM_USE_MPI
#define MPI_COMM_WORLD 0
#define MPI_Finalize(args...) {}
#define MPI_Init(args...) {}
#define MPI_Comm_rank(comm,rank) {*rank=0;}
#endif


namespace mfem {

  class Vector;
  class GridFunction;
  class DenseMatrix;
  class Mesh;
  class FiniteElementSpace;
  class FiniteElementCollection;
  class FiniteElement;
  class ElementTransformation;
  class IntegrationRule;
  class NonlinearForm;
  class NonlinearFormIntegrator;
  class ParGridFunction;
  class ParMesh;
  class ParFiniteElementSpace;
  class ParNonlinearForm;
  class NewtonSolver;
  class Solver;

}  // end of namespace mfem

namespace mfem_mgis {

  //! brief a simple alias
  using FiniteElementCollection = mfem::FiniteElementCollection;

  /*!
   * \brief a simple alias used to select the `MFEM` class handling the mesh
   * depending if a parallel computation is considered or not.
   * \tparam parallel: flag stating if a parallel computation is considered.
   */
  template <bool parallel>
  using Mesh = std::conditional_t<parallel, mfem::ParMesh, mfem::Mesh>;
  /*!
   * \brief a simple alias used to select the `MFEM` class handling the
   * finite element space depending if a parallel computation is considered or
   * not.
   * \tparam parallel: flag stating if a parallel computation is considered.
   */
  template <bool parallel>
  using FiniteElementSpace = std::conditional_t<parallel,
                                                mfem::ParFiniteElementSpace,
                                                mfem::FiniteElementSpace>;
  /*!
   * \brief a simple alias used to select the `MFEM` class representing a
   * non linear form depending if a parallel computation is considered or
   * not.
   * \tparam parallel: flag stating if a parallel computation is considered.
   */
  template <bool parallel>
  using NonlinearForm =
      std::conditional_t<parallel, mfem::ParNonlinearForm, mfem::NonlinearForm>;
  /*!
   * \brief a simple alias used to select the `MFEM` class representing a
   * non linear form integrator depending if a parallel computation is
   * considered or not.
   * \tparam parallel: flag stating if a parallel computation is considered.
   */
  using NonlinearFormIntegrator = mfem::NonlinearFormIntegrator;
  /*!
   * \brief a simple alias used to select the `MFEM` class representing a
   * grid function depending if a parallel computation is considered or not.
   * \tparam parallel: flag stating if a parallel computation is considered.
   */
  template <bool parallel>
  using GridFunction =
      std::conditional_t<parallel, mfem::ParGridFunction, mfem::GridFunction>;
  //! a simple alias
  using LinearSolver = mfem::Solver;
  //! a simple alias
  using NewtonSolver = mfem::NewtonSolver;

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MFEM_FORWARD_HXX */
