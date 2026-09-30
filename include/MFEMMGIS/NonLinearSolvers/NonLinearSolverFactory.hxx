/*!
 * \file   include/MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.hxx
 * \brief  This file declares the `NonLinearSolverFactory` class
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
     * \return the nonlinear solver
     */
    [[nodiscard]] virtual std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<true>& p,
        const Parameters& parameters) noexcept = 0;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief generate a nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] p: nonlinear problem
     * \param[in] parameters: parameters used to initialize the nonlinear
     * solver
     * \return the nonlinear solver
     */
    [[nodiscard]] virtual std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<false>& p,
        const Parameters& parameters) noexcept = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearSolverGenerator() noexcept;
  };  // end of AbstractNonLinearSolverGenerator

  //! \brief generator suitable for most nonlinear solvers
  template <std::derived_from<AbstractNonLinearSolver> SolverType>
  struct StandardNonLinearSolverGenerator : AbstractNonLinearSolverGenerator {
#ifdef MFEM_USE_MPI
    /*!
     * \brief generate a nonlinear solver of type `SolverType` and set its
     * parameters
     * \param[in, out] ctx: execution context
     * \param[in] p: nonlinear problem
     * \param[in] parameters: parameters used to initialize the nonlinear
     * solver
     * \return the nonlinear solver
     */
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<true>& p,
        const Parameters& parameters) noexcept override;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief generate a nonlinear solver of type `SolverType` and set its
     * parameters
     * \param[in, out] ctx: execution context
     * \param[in] p: nonlinear problem
     * \param[in] parameters: parameters used to initialize the nonlinear
     * solver
     * \return the nonlinear solver
     */
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> operator()(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<false>& p,
        const Parameters& parameters) noexcept override;
    //! \brief destructor
    ~StandardNonLinearSolverGenerator() noexcept override;
  };  // end of StandardNonLinearSolverGenerator

  //! \brief abstract factory for nonlinear solvers
  struct MFEM_MGIS_EXPORT NonLinearSolverFactory {
    //! \return the unique instance of the class
    [[nodiscard]] static NonLinearSolverFactory& get() noexcept;
    //
    //! \brief deleted move constructor
    NonLinearSolverFactory(NonLinearSolverFactory&&) = delete;
    //! \brief deleted copy constructor
    NonLinearSolverFactory(const NonLinearSolverFactory&) = delete;
    //! \brief deleted move assignment
    NonLinearSolverFactory& operator=(NonLinearSolverFactory&&) = delete;
    //! \brief deleted copy assignment
    NonLinearSolverFactory& operator=(const NonLinearSolverFactory&) = delete;
    /*!
     * \brief register a new nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the nonlinear solver
     * \param[in] g: generator of the nonlinear solver
     * \return true on success
     */
    [[nodiscard]] bool add(
        Context& ctx,
        std::string_view n,
        std::unique_ptr<AbstractNonLinearSolverGenerator> g) noexcept;
#ifdef MFEM_USE_MPI
    /*!
     * \brief generate a nonlinear solver
     * \return the requested nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the nonlinear solver
     * \param[in] p: non linear evolution problem
     * \param[in] parameters: parameters passed to the nonlinear solver
     */
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> generate(
        Context& ctx,
        std::string_view n,
        NonLinearEvolutionProblemImplementation<true>& p,
        const Parameters& parameters) const noexcept;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief generate a nonlinear solver
     * \return the requested nonlinear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the nonlinear solver
     * \param[in] p: non linear evolution problem
     * \param[in] parameters: parameters passed to the nonlinear solver
     */
    [[nodiscard]] std::unique_ptr<AbstractNonLinearSolver> generate(
        Context& ctx,
        std::string_view n,
        NonLinearEvolutionProblemImplementation<false>& p,
        const Parameters& parameters) const noexcept;

   private:
    //! \brief default constructor
    NonLinearSolverFactory() noexcept;
    //! \brief destructor
    ~NonLinearSolverFactory();
    //! \brief registered generators
    std::map<std::string,
             std::unique_ptr<AbstractNonLinearSolverGenerator>,
             std::less<>>
        generators;
  };  // end of struct NonLinearSolverFactory

}  // end of namespace mfem_mgis

#include "MFEMMGIS/NonLinearSolvers/NonLinearSolverFactory.ixx"

#endif /* LIB_MFEM_MGIS_NONLINEARSOLVERS_NONLINEARSOLVERFACTORY_HXX*/
