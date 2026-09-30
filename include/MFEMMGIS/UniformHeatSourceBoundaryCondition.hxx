/*!
 * \file   include/MFEMMGIS/UniformHeatSourceBoundaryCondition.hxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#ifndef LIB_MFEM_MGIS_UNIFORMHEATSOURCEBOUNDARYCONDITION_HXX
#define LIB_MFEM_MGIS_UNIFORMHEATSOURCEBOUNDARYCONDITION_HXX

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
   * \brief base class for non linear integrators based on an MGIS' behaviours.
   * This class manages an mapping associating a material and its identifier
   */
  struct MFEM_MGIS_EXPORT UniformHeatSourceBoundaryCondition final
      : public AbstractBoundaryCondition {
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     * \param[in] params: parameters defining the boundary condition
     */
    UniformHeatSourceBoundaryCondition(AbstractNonLinearEvolutionProblem& p,
                                       const Parameters& params);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] mid: material identifier
     * \param[in] prvalues: function returning the imposed values
     */
    UniformHeatSourceBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const size_type mid,
        std::function<real(const real)> prvalues);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] mid: material identifier
     * \param[in] prvalues: function returning the imposed values
     */
    UniformHeatSourceBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization> fed,
        const std::string_view mid,
        std::function<real(const real)> prvalues);
    //
#ifdef MFEM_USE_MPI
    [[nodiscard]] bool addNonlinearFormIntegrator(
        Context& ctx,
        NonlinearForm<true>& f,
        const mfem::Vector& u) noexcept override;
#endif /* MFEM_USE_MPI */
    [[nodiscard]] bool addNonlinearFormIntegrator(
        Context& ctx,
        NonlinearForm<false>& f,
        const mfem::Vector& u) noexcept override;
#ifdef MFEM_USE_MPI
    [[nodiscard]] bool addLinearFormIntegrators(
        Context& ctx,
        BilinearForm<true>& a,
        LinearForm<true>& b,
        const mfem::Vector& u,
        const real t,
        const real dt) noexcept override;
#endif /* MFEM_USE_MPI */
    [[nodiscard]] bool addLinearFormIntegrators(
        Context& ctx,
        BilinearForm<false>& a,
        LinearForm<false>& b,
        const mfem::Vector& u,
        const real t,
        const real dt) noexcept override;
    [[nodiscard]] bool setup(Context& ctx,
                             const real t,
                             const real dt) noexcept override;
    //! \brief destructor
    ~UniformHeatSourceBoundaryCondition() override;

   protected:
    //! \brief internal structure
    struct UniformHeatSourceFormIntegratorBase;
    //! \brief internal structure
    struct UniformHeatSourceLinearFormIntegrator;
    //! \brief internal structure
    struct UniformHeatSourceNonlinearFormIntegrator;
    //! \brief finite element discretization
    std::shared_ptr<FiniteElementDiscretization> finiteElementDiscretization;
    //! \brief list of material identifiers
    std::vector<size_type> mids;
    //
    mfem::Array<mfem_mgis::size_type> materials_markers;
    //! \brief function returning the value of the heat source
    std::function<real(const real)> qfct;
    //! \brief underlying integrator
    UniformHeatSourceNonlinearFormIntegrator* const nfi = nullptr;
    //
    bool shallFreeIntegrator = true;
  };  // end of UniformHeatSourceBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_UNIFORMHEATSOURCEBOUNDARYCONDITION_HXX */
