/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_EVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_EVOLUTIONPROBLEM_HXX

#include <memory>
#include "mfem/linalg/vector.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/nonlinearform.hpp"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"

namespace mfem_mgis {

  // forward declaration
  struct MGISIntegrator;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblem
      : public mfem::NonlinearForm {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fs: finite element space
     * \param[in] h: modelling hypothesis
     */
    NonLinearEvolutionProblem(std::shared_ptr<mfem::FiniteElementSpace>,
                              const Hypothesis);
    //! \return the Newton solver
    mfem::NewtonSolver& getSolver();
    /*! 
     * \brief solve the non linear problem over the given time step
     * \param[in] dt: time increment
     */
    void solve(const real);
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   private:
    //! \brief pointer to the underlying domain integrator
    MGISIntegrator* const mgis_integrator;
    //! \brief underlying finit element space
    const std::shared_ptr<const mfem::FiniteElementSpace> fe_space;
    //! \brief modelling hypothesis
    const Hypothesis hypothesis;
    //! \brief newton solver
    mfem::NewtonSolver solver;
    //! \brief unknowns at the beginning of the time step
    mfem::Vector u0;
    //! \brief unknowns at the end of the time step
    mfem::Vector u1;
  };  // end of struct NonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_EVOLUTIONPROBLEM */
