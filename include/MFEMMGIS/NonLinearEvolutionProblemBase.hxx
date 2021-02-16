/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_EVOLUTIONPROBLEMBASE_HXX
#define LIB_MFEM_MGIS_EVOLUTIONPROBLEMBASE_HXX

#include <memory>
#include "mfem/config/config.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/nonlinearform.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pnonlinearform.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemCommon.hxx"

namespace mfem_mgis {

  /*!
   * \brief class for solving non linear evolution problems
   */
  template <bool parallel>
  struct NonLinearEvolutionProblemBase;

#ifdef MFEM_USE_MPI

  /*!
   * \brief partial specialisation for parallel computations
   */
  template <>
  struct NonLinearEvolutionProblemBase<true>
      : public NonLinearEvolutionProblemCommon, public NonlinearForm<true> {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    NonLinearEvolutionProblemBase(std::shared_ptr<FiniteElementDiscretization>);
    //! \return the finite element space
    FiniteElementSpace<true>& getFiniteElementSpace();
    //! \return the finite element space
    const FiniteElementSpace<true>& getFiniteElementSpace() const;
    //! \return the Newton solver
    NewtonSolver& getSolver();
    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in] dt: time increment
     */
    virtual void solve(const real);
    //! \brief destructor
    ~NonLinearEvolutionProblemBase();

   protected:
    //! \brief newton solver
    NewtonSolver solver;
  };  // end of struct NonLinearEvolutionProblemBase<true>

#endif /* MFEM_USE_MPI */

  /*!
   * \brief partial specialisation for sequential computations
   */
  template <>
  struct NonLinearEvolutionProblemBase<false>
      : public NonLinearEvolutionProblemCommon, public NonlinearForm<false> {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    NonLinearEvolutionProblemBase(std::shared_ptr<FiniteElementDiscretization>);
    //! \return the finite element space
    FiniteElementSpace<false>& getFiniteElementSpace();
    //! \return the finite element space
    const FiniteElementSpace<false>& getFiniteElementSpace() const;
    //! \return the Newton solver
    NewtonSolver& getSolver();
    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in] dt: time increment
     */
    virtual void solve(const real);
    //! \brief destructor
    ~NonLinearEvolutionProblemBase();

   protected:
    //! \brief newton solver
    NewtonSolver solver;
  };  // end of struct NonLinearEvolutionProblemBase<false>

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_EVOLUTIONPROBLEMBASE */
