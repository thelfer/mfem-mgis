/*!
 * \file   include/MFEMMGIS/LinearSolverFactory.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEM_MGIS_LINEARSOLVERFACTORY_HXX
#define LIB_MFEM_MGIS_LINEARSOLVERFACTORY_HXX

#include <map>
#include <memory>
#include <functional>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/LinearSolverHandler.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief an abstract factory for linear solvers
   * \tparam parallel: boolean stating if parallel linear solvers are
   * considered
   *
   * \note if a linear solver is added, the `hasConverged` function in
   * `SolverUtilities.hxx` shall also be modified
   */
  template <bool parallel>
  struct LinearSolverFactory;

#ifdef MFEM_USE_MPI

  //! \brief specialisation in parallel
  template <>
  struct MFEM_MGIS_EXPORT LinearSolverFactory<true> {
    //! \brief a simple alias
    using Generator = std::function<LinearSolverHandler(
        Context&, FiniteElementSpace<true>&, const Parameters&)>;
    //! \return the unique instance of the class
    static LinearSolverFactory& getFactory();
    /*!
     * \brief register a new linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] g: generator of the linear solver
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           Generator g) noexcept;
    /*!
     * \return the requested linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] fespace: finite element space
     * \param[in] params: parameters passed to the linear solver
     */
    LinearSolverHandler generate(Context& ctx,
                                 std::string_view n,
                                 FiniteElementSpace<true>& fespace,
                                 const Parameters& params) const;

   private:
    //! \brief default constructor
    LinearSolverFactory();
    //! \brief destructor
    ~LinearSolverFactory();
    //! \brief registered generators
    std::map<std::string, Generator, std::less<>> generators;
  };  // end of struct LinearSolverFactory

#endif /* MFEM_USE_MPI */

  //! \brief specialisation in sequential
  template <>
  struct MFEM_MGIS_EXPORT LinearSolverFactory<false> {
    //! \brief a simple alias
    using Generator = std::function<LinearSolverHandler(
        Context&, FiniteElementSpace<false>&, const Parameters&)>;
    //! \return the unique instance of the class
    static LinearSolverFactory& getFactory();
    /*!
     * \brief register a new linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] g: generator of the linear solver
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           Generator g) noexcept;
    /*!
     * \return the requested linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] fespace: finite element space
     * \param[in] params: parameters passed to the linear solver
     */
    LinearSolverHandler generate(Context& ctx,
                                 std::string_view n,
                                 FiniteElementSpace<false>& fespace,
                                 const Parameters& params) const;

   private:
    //! \brief default constructor
    LinearSolverFactory();
    //! \brief destructor
    ~LinearSolverFactory();
    //! \brief registered generators
    std::map<std::string, Generator, std::less<>> generators;
  };  // end of struct LinearSolverFactory

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LINEARSOLVERFACTORY_HXX */
