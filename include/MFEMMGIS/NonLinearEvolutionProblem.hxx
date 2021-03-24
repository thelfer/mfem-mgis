/*!
 * \file   include/MFEMMGIS/NonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#ifndef LIB_MFEMMGIS_NONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEMMGIS_NONLINEAREVOLUTIONPROBLEM_HXX

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT NonLinearEvolutionProblem
      : AbstractNonLinearEvolutionProblem {
    //! \brief a simple alias
    using Hypothesis = mgis::behaviour::Hypothesis;
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] h: modelling hypothesis
     * \param[in] p: parameters
     */
    NonLinearEvolutionProblem(std::shared_ptr<FiniteElementDiscretization>,
                              const Hypothesis,
                              const Parameters& = Parameters());
    //! \return the implementation
    template <bool parallel>
    NonLinearEvolutionProblemImplementation<parallel>& getImplementation();
    //! \return the implementation
    template <bool parallel>
    const NonLinearEvolutionProblemImplementation<parallel> &getImplementation()
        const;
    //
    FiniteElementDiscretization &getFiniteElementDiscretization() override;
    std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() override;
    void setSolverParameters(const Parameters &) override;
    void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) override;
    void addPostProcessing(
        const std::function<void(const real, const real)>&) override;
    void addPostProcessing(std::string_view, const Parameters&) override;
    void executePostProcessings(const real, const real) override;
    void addBehaviourIntegrator(const std::string &,
                                const size_type,
                                const std::string &,
                                const std::string &)override;
    std::vector<size_type> getMaterialIdentifiers() const override;
    const Material &getMaterial(const size_type) const override;
    Material &getMaterial(const size_type)override;
    const BehaviourIntegrator &getBehaviourIntegrator(
        const size_type) const override;
    BehaviourIntegrator &getBehaviourIntegrator(const size_type) override;
    void revert() override;
    void update() override;
    void solve(const real, const real) override;
    //! \brief destructor
    ~NonLinearEvolutionProblem() override;

   protected:
    /*!
     * \brief method called before each resolution
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void setup(const real, const real);
    //! \brief implementation of the non linear problem
    std::unique_ptr<AbstractNonLinearEvolutionProblem> pimpl;
  };  // end of struct NonLinearEvolutionProblem

#ifdef MFEM_USE_MPI

  template <>
  NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation();

  template <>
  const NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation() const;

#else /* MFEM_USE_MPI */

  template <>
  [[noreturn]] NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation();

  template <>
  [[noreturn]] const NonLinearEvolutionProblemImplementation<true>
      &NonLinearEvolutionProblem::getImplementation() const;

#endif /* MFEM_USE_MPI */

  template <>
  NonLinearEvolutionProblemImplementation<false>
      &NonLinearEvolutionProblem::getImplementation();

  template <>
  const NonLinearEvolutionProblemImplementation<false>
      &NonLinearEvolutionProblem::getImplementation() const;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_NONLINEAREVOLUTIONPROBLEM_HXX */
