/*!
 * \file   include/MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx
 * \brief  This file declares the `UniformImposedPressureBoundaryCondition`
 * class
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#ifndef LIB_MFEM_MGIS_UNIFORMIMPOSEDPRESSUREBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_UNIFORMIMPOSEDPRESSUREBOUNDARYCONDITION_HXX

#include <memory>
#include <vector>
#include <mfem/linalg/densemat.hpp>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractBoundaryCondition.hxx"

namespace mfem_mgis {

  // forward declarations
  struct Parameters;
  struct FiniteElementDiscretization;
  struct AbstractNonLinearEvolutionProblem;

  /*!
   * \brief a uniform pressure imposed on a set of boundaries. Its value
   * evolves in time.
   */
  struct MFEM_MGIS_EXPORT UniformImposedPressureBoundaryCondition final
      : public AbstractBoundaryCondition {
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     * \param[in] params: parameters defining the boundary condition:
     * `Boundary` or `Boundaries`, and `LoadingEvolution`
     */
    UniformImposedPressureBoundaryCondition(
        AbstractNonLinearEvolutionProblem& p, const Parameters& params);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] bid: id of the boundary
     * \param[in] prvalues: function returning the imposed values
     */
    UniformImposedPressureBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const size_type bid,
        std::function<real(const real)> prvalues);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretization
     * \param[in] bid: regular expression selecting the boundaries by name
     * \param[in] prvalues: function returning the imposed values
     */
    UniformImposedPressureBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::string_view bid,
        std::function<real(const real)> prvalues);
#ifdef MFEM_USE_MPI
    /*!
     * \brief add the nonlinear form integrator describing the imposed
     * pressure
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] f: form
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     * \return true on success
     */
    [[nodiscard]] bool addNonlinearFormIntegrator(
        Context& ctx,
        NonlinearForm<true>& f,
        const mfem::Vector& u) noexcept override;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add the nonlinear form integrator describing the imposed
     * pressure
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] f: form
     * \param[in] u: current estimate of the solution at the end of the time
     * step
     * \return true on success
     */
    [[nodiscard]] bool addNonlinearFormIntegrator(
        Context& ctx,
        NonlinearForm<false>& f,
        const mfem::Vector& u) noexcept override;
#ifdef MFEM_USE_MPI
    /*!
     * \brief add the linear form integrator describing the imposed pressure
     * at the end of the time step
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] a: form computing the jacobian matrix
     * \param[in, out] b: form computing the right hand side
     * \param[in] u: solution at the beginning of the time step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool addLinearFormIntegrators(
        Context& ctx,
        BilinearForm<true>& a,
        LinearForm<true>& b,
        const mfem::Vector& u,
        const real t,
        const real dt) noexcept override;
#endif /* MFEM_USE_MPI */
    /*!
     * \brief add the linear form integrator describing the imposed pressure
     * at the end of the time step
     *
     * \param[in, out] ctx: execution context
     * \param[in, out] a: form computing the jacobian matrix
     * \param[in, out] b: form computing the right hand side
     * \param[in] u: solution at the beginning of the time step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool addLinearFormIntegrators(
        Context& ctx,
        BilinearForm<false>& a,
        LinearForm<false>& b,
        const mfem::Vector& u,
        const real t,
        const real dt) noexcept override;
    /*!
     * \brief set the pressure to its value at the end of the time step
     * \param[in, out] ctx: execution context
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool setup(Context& ctx,
                             const real t,
                             const real dt) noexcept override;
    //! \brief destructor
    ~UniformImposedPressureBoundaryCondition() override;

   protected:
    //! \brief internal structure
    struct UniformImposedPressureFormIntegratorBase;
    //! \brief linear form integrator imposing the pressure
    struct UniformImposedPressureLinearFormIntegrator;
    //! \brief nonlinear form integrator imposing the pressure
    struct UniformImposedPressureNonlinearFormIntegrator;
    //! \brief finite element discretization
    std::shared_ptr<FiniteElementDiscretization> finiteElementDiscretization;
    //! \brief list of boundary identifiers
    std::vector<size_type> bids;
    //! \brief markers of the boundaries
    mfem::Array<mfem_mgis::size_type> boundaries_markers;
    //! \brief function returning the value of the imposed pressure
    std::function<real(const real)> prfct;
    //! \brief underlying integrator
    UniformImposedPressureNonlinearFormIntegrator* const nfi = nullptr;
    //! \brief if the integrator must be freed by the destructor
    bool shallFreeIntegrator = true;
  };  // end of UniformImposedPressureBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_UNIFORMIMPOSEDPRESSUREBOUNDARYCONDITION_HXX */
