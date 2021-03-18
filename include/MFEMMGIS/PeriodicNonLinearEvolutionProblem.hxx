/*!
 * \file   include/MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date 11/12/2020
 */

#ifndef LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM_HXX

#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  /*!
   * \brief a base class handling the evolution of the macroscopic gradients
   */
  struct PeriodicNonLinearEvolutionProblemBase {
    /*!
     * \brief set the evolution of the macroscopic gradients
     * \param[in] e : function
     */
    virtual void setMacroscopicGradientsEvolution(
        const std::function<std::vector<real>(const real)>&);
    /*!
     * \return the value of the macroscopic gradients at the end of the time
     * step.
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual std::vector<real> getMacroscopicGradients(
        const real, const real) const;
    //! \brief destructor
    virtual ~PeriodicNonLinearEvolutionProblemBase();

   protected:
    //! \brief a function describing the evolution of the macroscopic gradients
    std::function<std::vector<real>(const real)>
        macroscopic_gradients_evolution;
  };  // end of PeriodicNonLinearEvolutionProblemBase

  /*!
   * \brief class for solving non linear evolution problems
   */
  template <bool parallel>
  struct PeriodicNonLinearEvolutionProblem;

#ifdef MFEM_USE_MPI

  template <>
  struct MFEM_MGIS_EXPORT PeriodicNonLinearEvolutionProblem<true>
      : public NonLinearEvolutionProblem<true>,
        public PeriodicNonLinearEvolutionProblemBase {
    /*!
     * \brief an helper function fix the degree of freedom of a point
     * \param[in] p: non linear evolution problem
     */
    static void setBoundaryConditions(
        mfem_mgis::NonLinearEvolutionProblemBase<true>&);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    PeriodicNonLinearEvolutionProblem(
        std::shared_ptr<FiniteElementDiscretization>);
    // disable adding boundary conditions
    [[noreturn]] void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) override;
    //! \brief destructor
    ~PeriodicNonLinearEvolutionProblem() override;

   protected:
    void setup(const real, const real) override;
  };  // end of struct PeriodicNonLinearEvolutionProblem

#endif /* MFEM_USE_MPI */

  template <>
  struct MFEM_MGIS_EXPORT PeriodicNonLinearEvolutionProblem<false>
      : public NonLinearEvolutionProblem<false>,
        public PeriodicNonLinearEvolutionProblemBase {
    /*!
     * \brief an helper function fix the degree of freedom of a point
     * \param[in] p: non linear evolution problem
     */
    static void setBoundaryConditions(
        mfem_mgis::NonLinearEvolutionProblemBase<false>&);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     */
    PeriodicNonLinearEvolutionProblem(
        std::shared_ptr<FiniteElementDiscretization>);
    // disable adding boundary conditions
    [[noreturn]] void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) override;
    //! \brief destructor
    ~PeriodicNonLinearEvolutionProblem() override;

   protected:
    void setup(const real, const real) override;
  };  // end of struct PeriodicNonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PERIODICNONLINEAREVOLUTIONPROBLEM */
