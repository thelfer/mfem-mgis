/*!
 * \file   include/MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEM_MGIS_NONLINEARSOLVERS_NONLINEARSOLVERFACTORY_HXX
#define LIB_MFEM_MGIS_NONLINEARSOLVERS_NONLINEARSOLVERFACTORY_HXX

#include <map>
#include <memory>
#include <functional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;
  struct AbstractNonLinearSolver;

  //! \brief interface class for generators of nonlinear solvers
  struct AbstractNonLinearSolverGenerator {
#ifdef MFEM_USE_MPI
    /*!
     * \brief generate a nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] p: nonlinear problem
     * \param[in] parameters: parameters used to initialize the nonlinear
     * solver
     */
    [[nodiscard]] virtual std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context&,
        NonLinearEvolutionProblemImplementation<true>&,
        const Parameters&) noexcept = 0;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief generate a nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] p: nonlinear problem
     * \param[in] parameters: parameters used to initialize the nonlinear
     * solver
     */
    [[nodiscard]] virtual std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context&,
        NonLinearEvolutionProblemImplementation<false>&,
        const Parameters&) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearSolverGenerator() noexcept;
  };  // end of AbstractNonLinearSolverGenerator

  //! \brief generator suitable for most nonlinear solvers
  template <std::derived_from<AbstractNonLinearSolver> SolverType>
  struct StandardNonLinearSolverGenerator : AbstractNonLinearSolverGenerator {
#ifdef MFEM_USE_MPI
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context&,
        NonLinearEvolutionProblemImplementation<true>&,
        const Parameters&) noexcept override;
#endif /* MFEM_USE_MPI */
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context&,
        NonLinearEvolutionProblemImplementation<false>&,
        const Parameters&) noexcept override;
    //! \brief destructor
    ~StandardNonLinearSolverGenerator() noexcept override;
  };  // end of StandardNonLinearSolverGenerator

  //! \brief abstract factory for nonlinear solvers
  struct MFEM_MGIS_EXPORT NonLinearSolverFactory {
    //! \return the unique instance of the class
    [[nodiscard]] static NonLinearSolverFactory& get() noexcept;
    //
    NonLinearSolverFactory(NonLinearSolverFactory&&) = delete;
    NonLinearSolverFactory(const NonLinearSolverFactory&) = delete;
    NonLinearSolverFactory& operator=(NonLinearSolverFactory&&) = delete;
    NonLinearSolverFactory& operator=(const NonLinearSolverFactory&) = delete;
    /*!
     * \brief register a new nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the nonlinear solver
     * \param[in] g: generator of the nonlinear solver
     */
    [[nodiscard]] bool add(
        Context&,
        std::string_view,
        std::unique_ptr<AbstractNonLinearSolverGenerator>) noexcept;
#ifdef MFEM_USE_MPI
    /*!
     * \return the requested nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the nonlinear solver
     * \param[in] p: non linear evolution postprocessing
     * \param[in] params: parameters passed to the nonlinear solver
     */
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> generate(
        Context&,
        std::string_view,
        NonLinearEvolutionProblemImplementation<true>&,
        const Parameters&) const noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \return the requested nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the nonlinear solver
     * \param[in] p: non linear evolution postprocessing
     * \param[in] params: parameters passed to the nonlinear solver
     */
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> generate(
        Context&,
        std::string_view,
        NonLinearEvolutionProblemImplementation<false>&,
        const Parameters&) const noexcept;

   private:
    //! \brief default destructor
    NonLinearSolverFactory() noexcept;
    //! \brief destructor
    ~NonLinearSolverFactory();
    //! \brief registred factories
    std::map<std::string,
             std::unique_ptr<AbstractNonLinearSolverGenerator>,
             std::less<>>
        generators;
  };  // end of struct NonLinearSolverFactory

}  // end of namespace mfem_mgis

#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.ixx"

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_NONLINEARSOLVERFACTORY_HXX*/
