/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblemBase.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_EVOLUTIONPROBLEMBASE_HXX
#define LIB_MFEM_MGIS_EVOLUTIONPROBLEMBASE_HXX

#include <memory>
#include "mfem/linalg/vector.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/nonlinearform.hpp"
#include "MGIS/Behaviour/Hypothesis.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemBase
      : public mfem::NonlinearForm {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     */
    NonLinearEvolutionProblemBase(std::shared_ptr<FiniteElementDiscretization>);
    //! \return the finite element space
    mfem::FiniteElementSpace& getFiniteElementSpace();
    //! \return the finite element space
    const mfem::FiniteElementSpace& getFiniteElementSpace() const;
    //! \return the Newton solver
    mfem::NewtonSolver& getSolver();
    //! \return the unknowns at the beginning of the time step
    mfem::Vector& getUnknownsAtBeginningOfTheTimeStep();
    //! \return the unknowns at the beginning of the time step
    const mfem::Vector& getUnknownsAtBeginningOfTheTimeStep() const;
    //! \return the unknowns at the end of the time step
    mfem::Vector& getUnknownsAtEndOfTheTimeStep();
    //! \return the unknowns at the end of the time step
    const mfem::Vector& getUnknownsAtEndOfTheTimeStep() const;
    /*!
     * \brief revert the state to the beginning of the time step.
     */
    virtual void revert();
    /*!
     * \brief update the state to the end of the time step.
     */
    virtual void update();

    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in] dt: time increment
     */
    virtual void solve(const real);

    //! \brief destructor
    ~NonLinearEvolutionProblemBase() override;

   protected:
    /*!
     * \brief set the time time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real);
    //! \brief underlying finite element discretization
    const std::shared_ptr<FiniteElementDiscretization> fe_discretization;
    //! \brief newton solver
    mfem::NewtonSolver solver;
    //! \brief unknowns at the beginning of the time step
    mfem::Vector u0;
    //! \brief unknowns at the end of the time step
    mfem::Vector u1;
  };  // end of struct NonLinearEvolutionProblemBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_EVOLUTIONPROBLEMBASE */
