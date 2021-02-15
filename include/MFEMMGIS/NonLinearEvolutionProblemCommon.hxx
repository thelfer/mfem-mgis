/*!
 * \file   NonLinearEvolutionProblemCommon.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   15/02/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMCOMMON_HXX
#define LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMCOMMON_HXX

#include <memory>
#include "mfem/linalg/vector.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

namespace mfem_mgis {

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblemCommon {
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    NonLinearEvolutionProblemCommon(std::shared_ptr<FiniteElementDiscretization>);
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

    //! \brief destructor
    virtual ~NonLinearEvolutionProblemCommon();

   protected:
    /*!
     * \brief set the time time increment
     * \param[in] dt: time increment
     */
    virtual void setTimeIncrement(const real);
    /*!
     * \brief method called before each resolution
     */
    virtual void setup();
    //! \brief underlying finite element discretization
    const std::shared_ptr<FiniteElementDiscretization> fe_discretization;
    //! \brief unknowns at the beginning of the time step
    mfem::Vector u0;
    //! \brief unknowns at the end of the time step
    mfem::Vector u1;
  };  // end of struct NonLinearEvolutionProblemCommon

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_NONLINEAREVOLUTIONPROBLEMCOMMON_HXX */
