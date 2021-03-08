/*!
 * \file   src/NonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   11/12/2020
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/PostProcessing.hxx"
#include "MFEMMGIS/MultiMaterialNonLinearIntegrator.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

  /*!
   * \brief post-processing defined by an std::function
   * \tparam parallel: boolean stating if the computations are performed in
   * parallel.
   */
  template <bool parallel>
  struct StdFunctionPostProcessing final : PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] fct: function executing the postprocessing
     */
    StdFunctionPostProcessing(
        const std::function<void(
            NonLinearEvolutionProblem<parallel>&, const real, const real)>& fct)
        : f(fct) {}  // end of StdFunctionPostProcessing
    //
    void execute(NonLinearEvolutionProblem<parallel>& p,
                 const real t,
                 const real dt) override {
      this->f(p, t, dt);
    }  // end of execute
    //! \brief destructor
    ~StdFunctionPostProcessing() override = default;

   private:
    //! \brief function executing the post-processing
    std::function<void(
        NonLinearEvolutionProblem<parallel>&, const real, const real)>
        f;
  };  // end of struct StdFunctionPostProcessing

#ifdef MFEM_USE_MPI

  NonLinearEvolutionProblem<true>::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      : NonLinearEvolutionProblemBase(fed),
        MultiMaterialEvolutionProblemBase(fed, h) {
    if (this->fe_discretization->getMesh<true>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblemBase::NonLinearEvolutionProblemBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    this->AddDomainIntegrator(this->mgis_integrator);
  }  // end of NonLinearEvolutionProblem

  void NonLinearEvolutionProblem<true>::addPostProcessing(
      std::unique_ptr<PostProcessing<true>> p) {
    this->postprocessings.push_back(std::move(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblem<true>::addPostProcessing(
      const std::function<
          void(NonLinearEvolutionProblem<true>&, const real, const real)>& p) {
    this->addPostProcessing(
        std::make_unique<StdFunctionPostProcessing<true>>(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblem<true>::executePostProcessings(const real t,
                                                               const real dt) {
    for (auto& p : this->postprocessings) {
      p->execute(*this, t, dt);
    }
  }  // end of executePostProcessings

  void NonLinearEvolutionProblem<true>::setup() {
    MultiMaterialEvolutionProblemBase::setup();
  }  // end of setup

  void NonLinearEvolutionProblem<true>::revert() {
    NonLinearEvolutionProblemBase<true>::revert();
    MultiMaterialEvolutionProblemBase::revert();
  }  // end of revert

  void NonLinearEvolutionProblem<true>::update() {
    NonLinearEvolutionProblemBase<true>::update();
    MultiMaterialEvolutionProblemBase::update();
  }  // end of update

  void NonLinearEvolutionProblem<true>::setTimeIncrement(const real dt) {
    MultiMaterialEvolutionProblemBase::setTimeIncrement(dt);
  }  // end of setTimeIncrement

  NonLinearEvolutionProblem<true>::~NonLinearEvolutionProblem() = default;

#endif /* MFEM_USE_MPI */

  NonLinearEvolutionProblem<false>::NonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed, const Hypothesis h)
      : NonLinearEvolutionProblemBase(fed),
        MultiMaterialEvolutionProblemBase(fed, h) {
    if (this->fe_discretization->getMesh<false>().Dimension() !=
        mgis::behaviour::getSpaceDimension(h)) {
      mgis::raise(
          "NonLinearEvolutionProblemBase::NonLinearEvolutionProblemBase: "
          "modelling hypothesis is not consistent with the spatial dimension "
          "of the mesh");
    }
    this->AddDomainIntegrator(this->mgis_integrator);
  }  // end of NonLinearEvolutionProblem

  void NonLinearEvolutionProblem<false>::addPostProcessing(
      std::unique_ptr<PostProcessing<false>> p) {
    this->postprocessings.push_back(std::move(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblem<false>::addPostProcessing(
      const std::function<
          void(NonLinearEvolutionProblem<false>&, const real, const real)>& p) {
    this->addPostProcessing(
        std::make_unique<StdFunctionPostProcessing<false>>(p));
  }  // end of addPostProcessing

  void NonLinearEvolutionProblem<false>::executePostProcessings(const real t,
                                                                const real dt) {
    for (auto& p : this->postprocessings) {
      p->execute(*this, t, dt);
    }
  }  // end of executePostProcessings

  void NonLinearEvolutionProblem<false>::setup() {
    MultiMaterialEvolutionProblemBase::setup();
  }  // end of setup

  void NonLinearEvolutionProblem<false>::revert() {
    NonLinearEvolutionProblemBase<false>::revert();
    MultiMaterialEvolutionProblemBase::revert();
  }  // end of revert

  void NonLinearEvolutionProblem<false>::update() {
    NonLinearEvolutionProblemBase<false>::update();
    MultiMaterialEvolutionProblemBase::update();
  }  // end of update

  void NonLinearEvolutionProblem<false>::setTimeIncrement(const real dt) {
    MultiMaterialEvolutionProblemBase::setTimeIncrement(dt);
  }  // end of setTimeIncrement

  NonLinearEvolutionProblem<false>::~NonLinearEvolutionProblem() = default;

}  // end of namespace mfem_mgis
