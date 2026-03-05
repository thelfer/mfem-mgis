/*!
 * \file   include/MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx
 * \brief
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
   * \brief base class for non linear integrators based on an MGIS' behaviours.
   * This class manages an mapping associating a material and its identifier
   */
  struct MFEM_MGIS_EXPORT UniformImposedPressureBoundaryCondition final
      : public AbstractBoundaryCondition {
    /*!
     * \brief constructor
     * \param[in] p: non linear evolution problem
     * \param[in] params: parameters defining the boundary condition
     */
    UniformImposedPressureBoundaryCondition(AbstractNonLinearEvolutionProblem&,
                                            const Parameters&);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] prvalues: function returning the imposed values
     */
    UniformImposedPressureBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization>,
        const size_type,
        std::function<real(const real)>);
    /*!
     * \brief constructor
     * \param[in] fed: finite element discretiszation
     * \param[in] bid: id of the boundary
     * \param[in] prvalues: function returning the imposed values
     */
    UniformImposedPressureBoundaryCondition(
        std::shared_ptr<FiniteElementDiscretization>,
        const std::string_view,
        std::function<real(const real)>);
    //
#ifdef MFEM_USE_MPI
    [[nodiscard]] bool addNonlinearFormIntegrator(
        Context&, NonlinearForm<true>&, const mfem::Vector&) noexcept override;
#endif /* MFEM_USE_MPI */
    [[nodiscard]] bool addNonlinearFormIntegrator(
        Context&, NonlinearForm<false>&, const mfem::Vector&) noexcept override;
#ifdef MFEM_USE_MPI
    [[nodiscard]] bool addLinearFormIntegrators(Context&,
                                                BilinearForm<true>&,
                                                LinearForm<true>&,
                                                const mfem::Vector&,
                                                const real,
                                                const real) noexcept override;
#endif /* MFEM_USE_MPI */
    [[nodiscard]] bool addLinearFormIntegrators(Context&,
                                                BilinearForm<false>&,
                                                LinearForm<false>&,
                                                const mfem::Vector&,
                                                const real,
                                                const real) noexcept override;
    void setup(const real, const real) override;
    //! \brief destructor
    virtual ~UniformImposedPressureBoundaryCondition();

   protected:
    //! \brief internal structure
    struct UniformImposedPressureFormIntegratorBase;
    //! \brief linear form to impose the internal structure
    struct UniformImposedPressureLinearFormIntegrator;
    //! \brief non linear form to impose the internal structure
    struct UniformImposedPressureNonlinearFormIntegrator;
    //! \brief finite element discretization
    std::shared_ptr<FiniteElementDiscretization> finiteElementDiscretization;
    //! \brief list of boundary identifiers
    std::vector<size_type> bids;
    //
    mfem::Array<mfem_mgis::size_type> boundaries_markers;
    //! \brief function returning the value of the imposed pressure
    std::function<real(const real)> prfct;
    //! \brief underlying integrator
    UniformImposedPressureNonlinearFormIntegrator* const nfi = nullptr;
    //
    bool shallFreeIntegrator = true;
  };  // end of UniformImposedPressureBoundaryCondition

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_UNIFORMIMPOSEDPRESSUREBOUNDARYCONDITION_HXX */
