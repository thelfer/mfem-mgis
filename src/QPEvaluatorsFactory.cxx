/*!
 * \file   src/QPEvaluatorsFactory.cxx
 * \brief  This file implements the QPEvaluatorsFactory class.
 * \author Thomas Helfer
 * \date   07/10/2022
 */

#include "MFEMMGIS/QPEvaluator/AbstractQPEvaluator.hxx"
#include "MFEMMGIS/QPEvaluator/AbstractQPEvaluatorGenerator.hxx"
#include "MFEMMGIS/QPEvaluator/QPEvaluatorDescription.hxx"
// #include "MFEMMGIS/QPEvaluator/QPEvaluatorsFactory.hxx"

namespace mfem_mgis::internal {

  //   /*!
  //    * \brief a generator for evaluators based on a field on integration
  //    points.
  //    */
  //   class NonUniformQPEvaluatorGenerator final
  //       : public AbstractQPEvaluatorGenerator {
  //    public:
  //     /*!
  //      * \brief constructor
  //      * \param[in] qid: quadrature
  //      * \param[in] f: reference to the integration point field.
  //      */
  //     NonUniformQPEvaluatorGenerator(const QuadId qid, const ConstIPFieldView
  //     &f)
  //         : ipfield_(f), quadId_(qid) {}  // end of
  //         NonUniformQPEvaluatorGenerator
  //     //
  //     [[nodiscard]] size_type getNumRows() const override {
  //       return this->ipfield_.getNumRows();
  //     }
  //     [[nodiscard]] size_type getNumCols() const override {
  //       return this->ipfield_.getNumCols();
  //     }
  //     [[nodiscard]] Vector<QPEvaluatorDescription> getDependencies()
  //         const override {
  //       return {};
  //     }
  //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> operator()(
  //         Context &ctx,
  //         const QPEvaluatorsFactory &,
  //         const MeshSet &m,
  //         const QuadId qid,
  //         const TimeStepStage) const override {
  //       if ((&m != &this->ipfield_.getMeshSet()) || (qid != this->quadId_)) {
  //         return ctx.registerErrorMessage(
  //             "inconsistent mesh set or quadrature id");
  //       }
  //       // we can't use MANTA_MAKE_SHARED_AS here because constructor is
  //       private if ((this->ipfield_.getNumRows() == 1) &&
  //           (this->ipfield_.getNumCols() == 1)) {
  //         try {
  //           return std::shared_ptr<AbstractQPEvaluator>(
  //               new NonUniformScalarQPEvaluator(this->ipfield_));
  //         } catch (...) {
  //           registerExceptionInErrorBacktrace(ctx);
  //         }
  //         return {};
  //       }
  //       try {
  //         return std::shared_ptr<AbstractQPEvaluator>(
  //             new NonUniformQPEvaluator(this->ipfield_));
  //       } catch (...) {
  //         registerExceptionInErrorBacktrace(ctx);
  //       }
  //       return {};
  //     }
  //     //! \brief destructor
  //     ~NonUniformQPEvaluatorGenerator() override = default;
  //
  //    private:
  //     //! \brief integration point field
  //     ConstIPFieldView ipfield_;
  //     /*!
  //      * \brief quadrature id
  //      *
  //      * \note This data member is required since the `FieldRef` method
  //      * does not have a `getQuadId` method. This will be treated later
  //      * in a dedicated issue.
  //      */
  //     const QuadId quadId_;
  //   };  // end of NonUniformQPEvaluatorGenerator
  //
  //   class ScalarQPEvaluatorsUnaryFunctionGenerator final
  //       : public AbstractQPEvaluatorGenerator {
  //    public:
  //     /*!
  //      * \brief constructor
  //      * \param[in] m: mesh set
  //      * \param[in] qid: quadrature identification
  //      * \param[in] f: function
  //      * \param[in] a: name of evaluator used as the first argument of the
  //      * function
  //      */
  //     ScalarQPEvaluatorsUnaryFunctionGenerator(
  //         const std::function<Real(const Real)> &f, const std::string &a)
  //         : function_(f),
  //           argument_(a) {}  // end of
  //           ScalarQPEvaluatorsUnaryFunctionGenerator
  //     //
  //     [[nodiscard]] size_type getNumRows() const override { return 1u; }
  //     [[nodiscard]] size_type getNumCols() const override { return 1u; }
  //     [[nodiscard]] Vector<QPEvaluatorDescription> getDependencies()
  //         const override {
  //       return {{this->argument_, 1, 1}};
  //     }
  //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> operator()(
  //         Context &ctx,
  //         const QPEvaluatorsFactory &g,
  //         const MeshSet &m,
  //         const QuadId qid,
  //         const TimeStepStage s) const override {
  //       const auto a = g.generate(ctx, m, qid, s, this->argument_);
  //       if (a.get() == nullptr) {
  //         return {};
  //       }
  //       if ((a->getNumRows() != 1u) || (a->getNumCols() != 1u)) {
  //         return ctx.registerErrorMessage("argument is not scalar");
  //       }
  //       return MANTA_MAKE_SHARED_AS(AbstractQPEvaluator,
  //                                   ScalarQPEvaluatorsUnaryFunction, ctx, m,
  //                                   qid, this->function_, a);
  //     }
  //     //! \brief destructor
  //     ~ScalarQPEvaluatorsUnaryFunctionGenerator() override = default;
  //
  //    private:
  //     //! \brief function
  //     std::function<Real(const Real)> function_;
  //     //! \brief argument
  //     std::string argument_;
  //   };  // end of ScalarQPEvaluatorsUnaryFunctionGenerator
  //
  //   class ScalarQPEvaluatorsBinaryFunctionGenerator final
  //       : public AbstractQPEvaluatorGenerator {
  //    public:
  //     /*!
  //      * \brief constructor
  //      * \param[in] m: mesh set
  //      * \param[in] qid: quadrature identification
  //      * \param[in] f: function
  //      * \param[in] a1: name of evaluator used as the first argument of the
  //      * function \param[in] a1: name of evaluator used as the second
  //      argument of
  //      * the function
  //      */
  //     ScalarQPEvaluatorsBinaryFunctionGenerator(
  //         const std::function<Real(const Real, const Real)> &f,
  //         const std::string &a1,
  //         const std::string &a2)
  //         : function_(f),
  //           firstArgument_(a1),
  //           secondArgument_(a2) {
  //     }  // end of ScalarQPEvaluatorsBinaryFunctionGenerator
  //     //
  //     [[nodiscard]] size_type getNumRows() const override { return 1u; }
  //     [[nodiscard]] size_type getNumCols() const override { return 1u; }
  //     [[nodiscard]] Vector<QPEvaluatorDescription> getDependencies()
  //         const override {
  //       return {{this->firstArgument_, 1, 1}, {this->secondArgument_, 1, 1}};
  //     }
  //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> operator()(
  //         Context &ctx,
  //         const QPEvaluatorsFactory &g,
  //         const MeshSet &m,
  //         const QuadId qid,
  //         const TimeStepStage s) const override {
  //       const auto a1 = g.generate(ctx, m, qid, s, this->firstArgument_);
  //       if (a1.get() == nullptr) {
  //         return {};
  //       }
  //       if ((a1->getNumRows() != 1) || (a1->getNumCols() != 1)) {
  //         return ctx.registerErrorMessage("first argument is not scalar");
  //       }
  //       const auto a2 = g.generate(ctx, m, qid, s, this->secondArgument_);
  //       if (a2.get() == nullptr) {
  //         return {};
  //       }
  //       if ((a2->getNumRows() != 1) || (a2->getNumCols() != 1)) {
  //         return ctx.registerErrorMessage("second argument is not scalar");
  //       }
  //       return MANTA_MAKE_SHARED_AS(AbstractQPEvaluator,
  //                                   ScalarQPEvaluatorsBinaryFunction, ctx, m,
  //                                   qid, this->function_, a1, a2);
  //     }
  //     //! \brief destructor
  //     ~ScalarQPEvaluatorsBinaryFunctionGenerator() override = default;
  //
  //    private:
  //     //! \brief function
  //     std::function<Real(const Real, const Real)> function_;
  //     //! \brief first argument
  //     std::string firstArgument_;
  //     //! \brief second argument
  //     std::string secondArgument_;
  //   };  // end of ScalarQPEvaluatorsBinaryFunctionGenerator
  //
  //   //! \brief class generating an uniform evaluator
  //   class UniformQPEvaluatorGenerator final
  //       : public AbstractQPEvaluatorGenerator {
  //    public:
  //     /*!
  //      * \brief constructor
  //      * \param[in] v: value of the evaluator.
  //      */
  //     UniformQPEvaluatorGenerator(const ConstMatrixRef<dynamic, dynamic> &v)
  //         : value_(v) {}  // end of UniformQPEvaluatorGenerator
  //     //
  //     [[nodiscard]] size_type getNumRows() const override { return
  //     rows(this->value_); }
  //     [[nodiscard]] size_type getNumCols() const override { return
  //     cols(this->value_); }
  //     [[nodiscard]] Vector<QPEvaluatorDescription> getDependencies()
  //         const override {
  //       return {};
  //     }
  //     [[nodiscard]] std::shared_ptr<AbstractQPEvaluator> operator()(
  //         Context &ctx,
  //         const QPEvaluatorsFactory &,
  //         const MeshSet &m,
  //         const QuadId qid,
  //         const TimeStepStage s) const override {
  //       return MANTA_MAKE_SHARED_AS(AbstractQPEvaluator, UniformQPEvaluator,
  //       ctx,
  //                                   m, qid, this->value_);
  //     }
  //     //! \brief destructor
  //     ~UniformQPEvaluatorGenerator() override = default;
  //
  //    private:
  //     //! \brief values of the evaluator
  //     const Matrix<dynamic, dynamic> value_;
  //   };  // end of UniformQPEvaluatorGenerator

}  // end of namespace mfem_mgis::internal

namespace mfem_mgis {
  //
  //   static ExitStatus merge(
  //       Context &ctx,
  //       QPEvaluatorsFactory::LocalDependenciesAnalysisOutput &r,
  //       const QPEvaluatorsFactory::LocalDependenciesAnalysisOutput &r2)
  //       noexcept {
  //     if (r2.shallContinue()) {
  //       return ExitStatus::success;
  //     }
  //     r.setStatus(r2);
  //     if (r2 == AdvancedExitStatus::unrecoverableError) {
  //       r.missingIPDependencies.clear();
  //       return ExitStatus::failure;
  //     }
  //     for (const auto &d2 : r2.missingIPDependencies) {
  //       auto p2 = r.missingIPDependencies.begin();
  //       auto p2e = r.missingIPDependencies.end();
  //       const auto p3 = std::find_if(p2, p2e, [&d2](const auto &dep) {
  //         return dep.getName() == d2.getName();
  //       });
  //       if (p3 == p2e) {
  //         r.missingIPDependencies.push_back(d2);
  //         continue;
  //       }
  //       if (p3->getNumRows() != d2.getNumRows()) {
  //         r.setStatus(AdvancedExitStatus::unrecoverableError);
  //         r.missingIPDependencies.clear();
  //         return ctx.registerErrorMessage(
  //             "dependency '" + d2.getName() +
  //             "' is multiply defined with inconsistent size");
  //       }
  //       if (p3->getNumCols() != d2.getNumCols()) {
  //         r.setStatus(AdvancedExitStatus::unrecoverableError);
  //         r.missingIPDependencies.clear();
  //         return ctx.registerErrorMessage(
  //             "dependency '" + d2.getName() +
  //             "' is multiply defined with inconsistent size");
  //       }
  //     }
  //     return ExitStatus::success;
  //   }  // end of merge
  //
  //   static InvalidResult reportMissingGenerator(Context &ctx,
  //                                               const MeshSet &m,
  //                                               const TimeStepStage s,
  //                                               const std::string &n)
  //                                               noexcept {
  //     return ctx.registerErrorMessage(
  //         "no generator named '" + n + "' declared on " +
  //         QPEvaluatorsFactory::getLocationDescription(m, s));
  //   }  // end of reportMissingGenerator
  //
  //   static InvalidResult reportMissingGenerator(Context &ctx,
  //                                               const MeshSet &m,
  //                                               const QuadId qid,
  //                                               const TimeStepStage s,
  //                                               const std::string &n)
  //                                               noexcept {
  //     return ctx.registerErrorMessage(
  //         "no generator named '" + n + "' declared on " +
  //         QPEvaluatorsFactory::getLocationDescription(m, qid, s));
  //   }  // end of reportMissingGenerator
  //
  //   std::string QPEvaluatorsFactory::getLocationDescription(
  //       const MeshSet &m, const QuadId qid, const TimeStepStage s) noexcept {
  //     const auto tss = [&s]() -> std::string {
  //       if (s == TimeStepStage::endOfTimeStep) {
  //         return "end";
  //       }
  //       return "beginning";
  //     };
  //     return "mesh set '" + m.getName() + "' for quadrature '" +
  //     toString(qid) +
  //            "' at the " + tss() + " of the time step";
  //   }  // end of getLocationDescription
  //
  //   std::string QPEvaluatorsFactory::getLocationDescription(
  //       const MeshSet &m, const TimeStepStage s) noexcept {
  //     const auto tss = [&s]() -> std::string {
  //       if (s == TimeStepStage::endOfTimeStep) {
  //         return "end";
  //       }
  //       return "beginning";
  //     };
  //     return "mesh set '" + m.getName() + "' for unspecified quadrature at
  //     the " +
  //            tss() + " of the time step";
  //   }  // end of getLocationDescription
  //
  //   QPEvaluatorsFactory::QPEvaluatorsFactory(ResourcesManager &m) noexcept
  //       : resourcesManager_(m) {}  // end of QPEvaluatorsFactory
  //
  //   ResourcesManager &QPEvaluatorsFactory::getResourcesManager() noexcept {
  //     return this->resourcesManager_;
  //   }  // end of getResourcesManager
  //
  //   const ResourcesManager &QPEvaluatorsFactory::getResourcesManager()
  //       const noexcept {
  //     return this->resourcesManager_;
  //   }  // end of getResourcesManager
  //
  //   QPEvaluatorsFactory::GeneratorsContainer &
  //   QPEvaluatorsFactory::getGeneratorsContainer_(const TimeStepStage s)
  //   noexcept {
  //     return (s == TimeStepStage::beginningOfTimeStep)
  //                ? this->generatorsAtTheBeginningOfTheTimeStep_
  //                : this->generatorsAtTheEndOfTheTimeStep_;
  //   }  // end of getGeneratorsContainer
  //
  //   const QPEvaluatorsFactory::GeneratorsContainer &
  //   QPEvaluatorsFactory::getGeneratorsContainer_(
  //       const TimeStepStage s) const noexcept {
  //     return (s == TimeStepStage::beginningOfTimeStep)
  //                ? this->generatorsAtTheBeginningOfTheTimeStep_
  //                : this->generatorsAtTheEndOfTheTimeStep_;
  //   }  // end of getGeneratorsContainer
  //
  //   ExitStatus QPEvaluatorsFactory::registerGenerator(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n,
  //       const Generator &g) noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       return ExitStatus::failure;
  //     }
  //     if (g.get() == nullptr) {
  //       ctx.registerErrorMessage("invalid generator (null pointer)");
  //       return ExitStatus::failure;
  //     }
  //     auto &generators = this->getGeneratorsContainer_(s);
  //     auto &map = generators[&m][qid];
  //     if (map.contains(n)) {
  //       return ctx.registerErrorMessage("generator '" + n + "' already
  //       declared");
  //     }
  //     map.insert({n, g});
  //     return ExitStatus::success;
  //   }  // end of registerGenerator
  //
  //   ExitStatus QPEvaluatorsFactory::registerGenerator(Context &ctx,
  //                                                     const MeshSet &m,
  //                                                     const QuadId qid,
  //                                                     const TimeStepStage s,
  //                                                     const std::string &n,
  //                                                     const Real v) noexcept
  //                                                     {
  //     auto g = MANTA_MAKE_SHARED_AS(AbstractQPEvaluatorGenerator,
  //                                   UniformScalarQPEvaluatorGenerator, ctx,
  //                                   v);
  //     return this->registerGenerator(ctx, m, qid, s, n, g);
  //   }  // end of registerGenerator
  //
  //   ExitStatus QPEvaluatorsFactory::registerGenerator(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n,
  //       const ConstMatrixRef<> &v) noexcept {
  //     using namespace ::manta::internal;
  //     auto g = MANTA_MAKE_SHARED_AS(AbstractQPEvaluatorGenerator,
  //                                   UniformQPEvaluatorGenerator, ctx, v);
  //     return this->registerGenerator(ctx, m, qid, s, n, g);
  //   }  // end of registerGenerator
  //
  //   ExitStatus QPEvaluatorsFactory::registerGenerator(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n,
  //       const std::function<Real(const Real)> &f,
  //       const std::string &a) noexcept {
  //     auto g = MANTA_MAKE_SHARED_AS(
  //         AbstractQPEvaluatorGenerator,
  //         ::manta::internal::ScalarQPEvaluatorsUnaryFunctionGenerator, ctx,
  //         f, a);
  //     return this->registerGenerator(ctx, m, qid, s, n, g);
  //   }  // end of registerGenerator
  //
  //   ExitStatus QPEvaluatorsFactory::registerGenerator(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n,
  //       const std::function<Real(const Real, const Real)> &f,
  //       const std::string &a1,
  //       const std::string &a2) noexcept {
  //     auto g = MANTA_MAKE_SHARED_AS(
  //         AbstractQPEvaluatorGenerator,
  //         ::manta::internal::ScalarQPEvaluatorsBinaryFunctionGenerator, ctx,
  //         f, a1, a2);
  //     return this->registerGenerator(ctx, m, qid, s, n, g);
  //   }  // end of registerGenerator
  //
  //   ExitStatus QPEvaluatorsFactory::registerGenerator(
  //       Context &ctx,
  //       const TimeStepStage s,
  //       const std::string &n,
  //       const ConstIPFieldView &f) noexcept {
  //     const auto qid = f.getQuadId();
  //     auto g = MANTA_MAKE_SHARED_AS(
  //         AbstractQPEvaluatorGenerator,
  //         ::manta::internal::NonUniformQPEvaluatorGenerator, ctx, qid, f);
  //     return this->registerGenerator(ctx, f.getMeshSet(), qid, s, n, g);
  //   }
  //
  //   Bool QPEvaluatorsFactory::containsGenerator(
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n) const noexcept {
  //     const auto &generators = this->getGeneratorsContainer_(s);
  //     if (!generators.contains(&m)) {
  //       return false;
  //     }
  //     const auto &map = generators.at(&m);
  //     if (!map.contains(qid)) {
  //       return false;
  //     }
  //     return map.at(qid).contains(n);
  //   }  // end of contains
  //
  //   std::shared_ptr<AbstractQPEvaluator> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const QPEvaluatorDescription &d) const noexcept {
  //     auto report_inconsistency_check = [&ctx, &m, &qid, &s,
  //                                        &d](const char *const reason) {
  //       ctx.registerErrorMessage(
  //           "invalid evaluator '" + d.getName() + "' generated on " +
  //           getLocationDescription(m, qid, s) + ": " + std::string{reason});
  //     };
  //     if (d.shallNotBeUsedInEvaluatorsGeneration()) {
  //       report_inconsistency_check("dummy description given");
  //     }
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       return {};
  //     }
  //     auto p = this->generate(ctx, m, qid, s, d.getName());
  //     if (p.get() == nullptr) {
  //       return p;
  //     }
  //     if (p->getNumRows() != d.getNumRows()) {
  //       report_inconsistency_check("invalid number of rows");
  //       return {};
  //     }
  //     if (p->getNumCols() != d.getNumCols()) {
  //       report_inconsistency_check("invalid number of columns");
  //       return {};
  //     }
  //     return p;
  //   }  // end of generate
  //
  //   std::shared_ptr<AbstractQPEvaluator> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n) const noexcept {
  //     auto report_inconsistency_check = [&ctx, &m, &qid, &s,
  //                                        &n](const char *const reason) {
  //       ctx.registerErrorMessage("invalid evaluator '" + n + "' generated on
  //       " +
  //                                getLocationDescription(m, qid, s) + ": " +
  //                                std::string{reason});
  //     };
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       return {};
  //     }
  //     const auto &generators = this->getGeneratorsContainer_(s);
  //     if (!generators.contains(&m)) {
  //       return reportMissingGenerator(ctx, m, qid, s, n);
  //     }
  //     const auto &map = generators.at(&m);
  //     if (!map.contains(qid)) {
  //       return reportMissingGenerator(ctx, m, qid, s, n);
  //     }
  //     const auto &map2 = map.at(qid);
  //     if (!map2.contains(n)) {
  //       return reportMissingGenerator(ctx, m, qid, s, n);
  //     }
  //     const auto r = this->analyseDependencies_(ctx, map2, n, {});
  //     if (r != AdvancedExitStatus::success) {
  //       if (!r.missingIPDependencies.empty()) {
  //         auto msg = "generation of evaluator '" + n + "' on " +
  //                    getLocationDescription(m, qid, s) + " failed. ";
  //         msg += "The following dependencies are missing:";
  //         for (const auto &d : r.missingIPDependencies) {
  //           msg += "\n- " + d.getName() + " of ";
  //           msg += "size (" + std::to_string(d.getNumRows()) + ", " +
  //                  std::to_string(d.getNumCols()) + ")";
  //         }
  //         ctx.registerErrorMessage(msg);
  //       }
  //       return {};
  //     }
  //     const auto &g = *(map2.at(n));
  //     const auto p = g(ctx, *this, m, qid, s);
  //     if (p.get() == nullptr) {
  //       ctx.registerErrorMessage("generation of evaluator '" + n + "' on " +
  //                                getLocationDescription(m, qid, s) + "
  //                                failed");
  //       return p;
  //     }
  //     if (&(p->getMeshSet()) != &m) {
  //       report_inconsistency_check("inconsistent mesh set");
  //       return {};
  //     }
  //     if (p->getQuadId() != qid) {
  //       report_inconsistency_check("inconsistent quadrature id");
  //       return {};
  //     }
  //     if (p->getNumRows() != g.getNumRows()) {
  //       report_inconsistency_check("invalid number of rows");
  //       return {};
  //     }
  //     if (p->getNumCols() != g.getNumCols()) {
  //       report_inconsistency_check("invalid number of columns");
  //       return {};
  //     }
  //     return p;
  //   }  // end of generate
  //
  //   Optional<Vector<std::pair<size_type,
  //   std::shared_ptr<AbstractQPEvaluator>>>> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const Vector<QPEvaluatorDescription> &descriptions) const noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       return {};
  //     }
  //     auto r = Vector<std::pair<size_type,
  //     std::shared_ptr<AbstractQPEvaluator>>>{}; auto o = size_type{}; for
  //     (const auto &d : descriptions) {
  //       if (d.shallNotBeUsedInEvaluatorsGeneration()) {
  //         o += d.getNumRows() * d.getNumCols();
  //         continue;
  //       }
  //       const auto g = this->generate(ctx, m, qid, s, d);
  //       if (g.get() == nullptr) {
  //         return {};
  //       }
  //       r.push_back({o, g});
  //       o += g->getNumRows() * g->getNumCols();
  //     }
  //     return r;
  //   }  // end of generate
  //
  //   Optional<Vector<std::pair<size_type,
  //   std::shared_ptr<AbstractQPEvaluator>>>> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const Vector<std::string> &names) const noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       return {};
  //     }
  //     auto r = Vector<std::pair<size_type,
  //     std::shared_ptr<AbstractQPEvaluator>>>{}; auto o = size_type{}; for
  //     (const auto &n : names) {
  //       const auto g = this->generate(ctx, m, qid, s, n);
  //       if (g.get() == nullptr) {
  //         return {};
  //       }
  //       r.push_back({o, g});
  //       o += g->getNumRows() * g->getNumCols();
  //     }
  //     return r;
  //   }  // end of generate
  //
  //   std::shared_ptr<AbstractQPEvaluator> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const TimeStepStage s,
  //       const QPEvaluatorDescription &d) const noexcept {
  //     auto report_inconsistency_check =
  //         [&ctx, &m, &s, &d](const QuadId qid, const char *const reason) {
  //           ctx.registerErrorMessage(
  //               "invalid evaluator '" + d.getName() + "' generated on " +
  //               getLocationDescription(m, qid, s) + ": " +
  //               std::string{reason});
  //         };
  //     if (d.shallNotBeUsedInEvaluatorsGeneration()) {
  //       return ctx.registerErrorMessage(
  //           "invalid evaluator '" + d.getName() + "' generated on " +
  //           getLocationDescription(m, s) + ": dummy evaluator description
  //           given");
  //     }
  //     if (!this->checkMeshSetConsistency_(ctx, m)) {
  //       return {};
  //     }
  //     auto p = this->generate(ctx, m, s, d.getName());
  //     if (p.get() == nullptr) {
  //       return p;
  //     }
  //     if (p->getNumRows() != d.getNumRows()) {
  //       report_inconsistency_check(p->getQuadId(), "invalid number of rows");
  //       return {};
  //     }
  //     if (p->getNumCols() != d.getNumCols()) {
  //       report_inconsistency_check(p->getQuadId(), "invalid number of
  //       columns"); return {};
  //     }
  //     return p;
  //   }  // end of generate
  //
  //   std::shared_ptr<AbstractQPEvaluator> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const TimeStepStage s,
  //       const std::string &n) const noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m)) {
  //       return {};
  //     }
  //     const auto &generators = this->getGeneratorsContainer_(s);
  //     if (!generators.contains(&m)) {
  //       return reportMissingGenerator(ctx, m, s, n);
  //     }
  //     const auto &map = generators.at(&m);
  //     auto oqid = Optional<QuadId>{};
  //     for (const auto &[qid, mg] : map) {
  //       if (mg.contains(n)) {
  //         if (isValid(oqid)) {
  //           // ambiguous case
  //           return ctx.registerErrorMessage(
  //               "the ip evaluator '" + n +
  //               "' is available for several quadratures.");
  //         }
  //         oqid = qid;
  //       }
  //     }
  //     if (isInvalid(oqid)) {
  //       return reportMissingGenerator(ctx, m, s, n);
  //     }
  //     return this->generate(ctx, m, *oqid, s, n);
  //   }  // end of generate
  //
  //   Optional<Vector<std::pair<size_type,
  //   std::shared_ptr<AbstractQPEvaluator>>>> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const TimeStepStage s,
  //       const Vector<QPEvaluatorDescription> &descriptions) const noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m)) {
  //       return {};
  //     }
  //     auto r = Vector<std::pair<size_type,
  //     std::shared_ptr<AbstractQPEvaluator>>>{}; auto o = size_type{}; for
  //     (const auto &d : descriptions) {
  //       if (d.shallNotBeUsedInEvaluatorsGeneration()) {
  //         o += d.getNumRows() * d.getNumCols();
  //         continue;
  //       }
  //       const auto g = this->generate(ctx, m, s, d);
  //       if (g.get() == nullptr) {
  //         return {};
  //       }
  //       r.push_back({o, g});
  //       o += g->getNumRows() * g->getNumCols();
  //     }
  //     return r;
  //   }  // end of generate
  //
  //   Optional<Vector<std::pair<size_type,
  //   std::shared_ptr<AbstractQPEvaluator>>>> QPEvaluatorsFactory::generate(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const TimeStepStage s,
  //       const Vector<std::string> &names) const noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m)) {
  //       return {};
  //     }
  //     auto r = Vector<std::pair<size_type,
  //     std::shared_ptr<AbstractQPEvaluator>>>{}; auto o = size_type{}; for
  //     (const auto &n : names) {
  //       const auto g = this->generate(ctx, m, s, n);
  //       if (g.get() == nullptr) {
  //         return {};
  //       }
  //       r.push_back({o, g});
  //       o += g->getNumRows() * g->getNumCols();
  //     }
  //     return r;
  //   }  // end of generate
  //
  //   QPEvaluatorsFactory::LocalDependenciesAnalysisOutput
  //   QPEvaluatorsFactory::analyseDependencies(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const QPEvaluatorDescription &d) const noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       auto r = LocalDependenciesAnalysisOutput{};
  //       r.setStatus(AdvancedExitStatus::unrecoverableError);
  //       return r;
  //     }
  //     const auto &generators = this->getGeneratorsContainer_(s);
  //     if (!generators.contains(&m) == 0) {
  //       auto r = LocalDependenciesAnalysisOutput{};
  //       r.setStatus(AdvancedExitStatus::recoverableError);
  //       r.missingIPDependencies.push_back(d);
  //       return r;
  //     }
  //     const auto &map = generators.at(&m);
  //     if (!map.contains(qid)) {
  //       auto r = LocalDependenciesAnalysisOutput{};
  //       r.setStatus(AdvancedExitStatus::recoverableError);
  //       r.missingIPDependencies.push_back(d);
  //       return r;
  //     }
  //     const auto &map2 = map.at(qid);
  //     const auto r = this->analyseDependencies_(ctx, map2, d.getName(), {});
  //     if (r == AdvancedExitStatus::unrecoverableError) {
  //       ctx.registerErrorMessage("dependency resolution for evaluator '" +
  //                                d.getName() + "' on " +
  //                                getLocationDescription(m, qid, s) + "
  //                                failed");
  //     }
  //     return r;
  //   }  // end of analyseDependencies
  //
  //   QPEvaluatorsFactory::LocalDependenciesAnalysisOutput
  //   QPEvaluatorsFactory::analyseDependencies(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const Vector<QPEvaluatorDescription> &deps) const noexcept {
  //     auto r = LocalDependenciesAnalysisOutput{};
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       r.setStatus(AdvancedExitStatus::unrecoverableError);
  //       return r;
  //     }
  //     if (deps.empty()) {
  //       r.setStatus(AdvancedExitStatus::success);
  //       return r;
  //     }
  //     auto reportMissingGenerators = [&deps, &r] {
  //       r.setStatus(AdvancedExitStatus::recoverableError);
  //       for (const auto &d : deps) {
  //         r.missingIPDependencies.push_back(d);
  //       }
  //       return r;
  //     };
  //     const auto &generators = this->getGeneratorsContainer_(s);
  //     if (!generators.contains(&m)) {
  //       return reportMissingGenerators();
  //     }
  //     const auto &map = generators.at(&m);
  //     if (!map.contains(qid)) {
  //       return reportMissingGenerators();
  //     }
  //     const auto &map2 = map.at(qid);
  //     for (const auto &d : deps) {
  //       const auto &n = d.getName();
  //       if (map2.count(n) == 0) {
  //         r.setStatus(AdvancedExitStatus::recoverableError);
  //         r.missingIPDependencies.push_back(d);
  //       } else {
  //         if (!merge(ctx, r, this->analyseDependencies_(ctx, map2, n, {}))) {
  //           r.missingIPDependencies.clear();
  //           r.setStatus(AdvancedExitStatus::unrecoverableError);
  //           return r;
  //         }
  //       }
  //     }
  //     return r;
  //   }  // end of analyseDependencies
  //
  //   QPEvaluatorsFactory::LocalDependenciesAnalysisOutput
  //   QPEvaluatorsFactory::analyseDependencies(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n) const noexcept {
  //     auto r = LocalDependenciesAnalysisOutput{};
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       r.setStatus(AdvancedExitStatus::unrecoverableError);
  //       return r;
  //     }
  //     const auto &generators = this->getGeneratorsContainer_(s);
  //     if (!generators.contains(&m)) {
  //       return static_cast<ExitStatus>(reportMissingGenerator(ctx, m, qid, s,
  //       n));
  //     }
  //     const auto &map = generators.at(&m);
  //     if (!map.contains(qid)) {
  //       return static_cast<ExitStatus>(reportMissingGenerator(ctx, m, qid, s,
  //       n));
  //     }
  //     const auto &map2 = map.at(qid);
  //     r = this->analyseDependencies_(ctx, map2, n, {});
  //     if (r == AdvancedExitStatus::unrecoverableError) {
  //       ctx.registerErrorMessage("dependency resolution for evaluator '" + n
  //       +
  //                                "' on " + getLocationDescription(m, qid, s)
  //                                + " failed");
  //     }
  //     return r;
  //   }  // end of analyseDependencies
  //
  //   QPEvaluatorsFactory::LocalDependenciesAnalysisOutput
  //   QPEvaluatorsFactory::analyseDependencies(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const Vector<std::string> &names) const noexcept {
  //     auto r = LocalDependenciesAnalysisOutput{};
  //     if (!this->checkMeshSetConsistency_(ctx, m, qid)) {
  //       r.setStatus(AdvancedExitStatus::unrecoverableError);
  //       return r;
  //     }
  //     if (names.empty()) {
  //       r.setStatus(AdvancedExitStatus::success);
  //       return r;
  //     }
  //     const auto &generators = this->getGeneratorsContainer_(s);
  //     if (!generators.contains(&m)) {
  //       return static_cast<ExitStatus>(
  //           reportMissingGenerator(ctx, m, qid, s, names[0]));
  //     }
  //     const auto &map = generators.at(&m);
  //     if (!map.contains(qid)) {
  //       return static_cast<ExitStatus>(
  //           reportMissingGenerator(ctx, m, qid, s, names[0]));
  //     }
  //     const auto &map2 = map.at(qid);
  //     for (const auto &n : names) {
  //       if (!merge(ctx, r, this->analyseDependencies_(ctx, map2, n, {}))) {
  //         r.missingIPDependencies.clear();
  //         r.setStatus(AdvancedExitStatus::unrecoverableError);
  //         return r;
  //       }
  //     }
  //     return r;
  //   }  // end of analyseDependencies
  //
  //   QPEvaluatorsFactory::LocalDependenciesAnalysisOutput
  //   QPEvaluatorsFactory::analyseDependencies_(
  //       Context &ctx,
  //       const HMap<std::string, Generator> &m,
  //       const std::string &n,
  //       const Vector<QPEvaluatorDescription> &previous_dependencies)
  //       const noexcept {
  //     auto r = LocalDependenciesAnalysisOutput{};
  //     if (!m.contains(n)) {
  //       ctx.registerErrorMessage("no generator named '" + n + "' declared");
  //       r.setStatus(AdvancedExitStatus::unrecoverableError);
  //       return r;
  //     }
  //     const auto &g = *(m.at(n));
  //     const auto current_dependencies = [&previous_dependencies, &n, &g] {
  //       auto tmp = previous_dependencies;
  //       tmp.emplace_back(n, g.getNumRows(), g.getNumCols());
  //       return tmp;
  //     }();
  //     for (const auto &d : g.getDependencies()) {
  //       const auto dn = d.getName();
  //       const auto p = current_dependencies.begin();
  //       const auto pe = current_dependencies.end();
  //       if (std::find_if(p, pe, [&dn](const auto &d2) {
  //             return d2.getName() == dn;
  //           }) != pe) {
  //         ctx.registerErrorMessage("cyclic dependency detected for evalutor
  //         '" +
  //                                  n + "'");
  //         r.setStatus(AdvancedExitStatus::unrecoverableError);
  //         return r;
  //       }
  //       if (!m.contains(dn)) {
  //         r.setStatus(AdvancedExitStatus::recoverableError);
  //         r.missingIPDependencies.push_back(d);
  //         continue;
  //       }
  //       const auto &g2 = *(m.at(dn));
  //       if (g2.getNumRows() != d.getNumRows()) {
  //         ctx.registerErrorMessage("dependency '" + dn + "' for evaluator '"
  //         + n +
  //                                  "' has unexpected number of rows (" +
  //                                  std::to_string(g2.getNumRows()) + " vs " +
  //                                  std::to_string(d.getNumRows()) + ")");
  //         r.setStatus(AdvancedExitStatus::unrecoverableError);
  //         return r;
  //       }
  //       if (g2.getNumCols() != d.getNumCols()) {
  //         ctx.registerErrorMessage("dependency '" + dn + "'for evaluator '" +
  //         n +
  //                                  "' has unexpected number of columns (" +
  //                                  std::to_string(g2.getNumCols()) + " vs " +
  //                                  std::to_string(d.getNumCols()) + ")");
  //         r.setStatus(AdvancedExitStatus::unrecoverableError);
  //         return r;
  //       }
  //       if (!merge(
  //               ctx, r,
  //               this->analyseDependencies_(ctx, m, dn,
  //               current_dependencies))) {
  //         r.setStatus(AdvancedExitStatus::unrecoverableError);
  //         return r;
  //       }
  //     }
  //     return r;
  //   }  // end of analyseDependencies_
  //
  //   ExitStatus QPEvaluatorsFactory::checkMeshSetConsistency_(
  //       Context &ctx, const MeshSet &m) const noexcept {
  //     const auto &m3 =
  //         this->getGeneratorsContainer_(TimeStepStage::beginningOfTimeStep);
  //     const auto &m4 =
  //         this->getGeneratorsContainer_(TimeStepStage::endOfTimeStep);
  //     if ((m3.contains(&m)) || (m4.contains(&m))) {
  //       // check already performend
  //       return ExitStatus::success;
  //     }
  //     if (&(this->resourcesManager_) != &(m.getResourcesManager())) {
  //       return ctx.registerErrorMessage("the mesh set '" + m.getName() +
  //                                       "' is defined on an inconsistent
  //                                       model");
  //     }
  //     // check if the mesh only contains one geometric type
  //     const auto otype = getGeometricType(ctx, m);
  //     if (!otype.hasValue()) {
  //       return ExitStatus::failure;
  //     }
  //     return ExitStatus::success;
  //   }  // end of checkMeshSetConsistency_
  //
  //   ExitStatus QPEvaluatorsFactory::checkMeshSetConsistency_(
  //       Context &ctx, const MeshSet &m, const QuadId qid) const noexcept {
  //     if (!this->checkMeshSetConsistency_(ctx, m)) {
  //       return ExitStatus::failure;
  //     }
  //     // check if the mesh only contains one geometric type
  //     const auto otype = getGeometricType(ctx, m);
  //     if (!otype.hasValue()) {
  //       return ExitStatus::failure;
  //     }
  //     const auto gtype = getGeometricSupport(*otype);
  //     if (gtype != Quadrature::getGeomSupport(qid)) {
  //       return ctx.registerErrorMessage(
  //           "the quadrature '" + toString(qid) +
  //           "'is not compatible with the geometric type of the mesh set '" +
  //           m.getName() + "'");
  //     }
  //     return ExitStatus::success;
  //   }  // end of checkMeshSetConsistency_
  //
  //   ExitStatus QPEvaluatorsFactory::checkIPFieldViewSpecifications_(
  //       Context &ctx,
  //       const MeshSet &m,
  //       const QuadId qid,
  //       const TimeStepStage s,
  //       const std::string &n,
  //       const std::string &fn,
  //       const size_type r,
  //       const size_type c,
  //       const IPFieldViewStorageSpecifications &bs) noexcept {
  //     auto report = [&ctx, &m, qid, s, &n, &fn](const char *const details) {
  //       auto msg = "error while generating view '" + n + "' on field '" + fn
  //       +
  //                  "' on " +
  //                  QPEvaluatorsFactory::getLocationDescription(m, qid, s) +
  //                  ": ";
  //       msg += details;
  //       ctx.registerErrorMessage(msg);
  //       return ExitStatus::failure;
  //     };
  //     if ((!(bs.block[0] > 0)) || (!(bs.block[1] > 0))) {
  //       return report("invalid view size");
  //     }
  //     if ((!(bs.blockTL[0] >= 0)) || (!(bs.blockTL[1] >= 0))) {
  //       return report("invalid offset");
  //     }
  //     if (bs.blockTL[0] + bs.block[0] > r) {
  //       return report("invalid number of rows");
  //     }
  //     if (bs.blockTL[1] + bs.block[1] > c) {
  //       return report("invalid number of columns");
  //     }
  //     return ExitStatus::success;
  //   }  // end of checkIPFieldViewSpecifications_
  //
  //   std::string QPEvaluatorsFactory::getRegisteredGeneratorsList()
  //       const noexcept {
  //     const auto l1 =
  //         this->getRegisteredGeneratorsList_(TimeStepStage::endOfTimeStep);
  //     const auto l2 =
  //         this->getRegisteredGeneratorsList_(TimeStepStage::endOfTimeStep);
  //     if (l2.empty()) {
  //       return l1;
  //     }
  //     if (l1.empty()) {
  //       return l2;
  //     }
  //     return l1 + '\n' + l2;
  //   }  // end of getRegisteredGeneratorsList
  //
  //   std::string QPEvaluatorsFactory::getRegisteredGeneratorsList_(
  //       const TimeStepStage ts) const noexcept {
  //     auto r = std::string{};
  //     for (const auto &[m, map] : this->getGeneratorsContainer_(ts)) {
  //       for (const auto &[qid, generators] : map) {
  //         if (!r.empty()) {
  //           r += '\n';
  //         }
  //         r += "- generators defined on " + getLocationDescription(*m, qid,
  //         ts) +
  //              ":";
  //         for (const auto &[n, g] : generators) {
  //           r += "\n  - " + n;
  //         }
  //       }
  //     }
  //     return r;
  //   }

  //  QPEvaluatorsFactory::~QPEvaluatorsFactory() noexcept = default;

}  // namespace mfem_mgis
