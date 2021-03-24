/*!
 * \file   include/MFEMMGIS/AbstractNonLinearEvolutionProblem.hxx
 * \brief
 * \author Thomas Helfer
 * \date   23/03/2021
 */

#ifndef LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEM_HXX
#define LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEM_HXX

#include <memory>
#include <functional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct FiniteElementDiscretization;
  // forward declaration
  struct DirichletBoundaryCondition;
  // forward declaration
  struct BehaviourIntegrator;
  // forward declaration
  struct Material;

  /*!
   * \brief class for solving non linear evolution problems
   */
  struct MFEM_MGIS_EXPORT AbstractNonLinearEvolutionProblem {
    //! \return the underlying finite element discretization
    virtual FiniteElementDiscretization &getFiniteElementDiscretization() = 0;
    //! \return the underlying finite element discretization
    virtual std::shared_ptr<FiniteElementDiscretization>
    getFiniteElementDiscretizationPointer() = 0;
    //! \return the unknowns at the beginning of the time step
    virtual mfem::Vector &getUnknownsAtBeginningOfTheTimeStep() = 0;
    //! \return the unknowns at the beginning of the time step
    virtual const mfem::Vector &getUnknownsAtBeginningOfTheTimeStep() const = 0;
    //! \return the unknowns at the end of the time step
    virtual mfem::Vector &getUnknownsAtEndOfTheTimeStep() = 0;
    //! \return the unknowns at the end of the time step
    virtual const mfem::Vector &getUnknownsAtEndOfTheTimeStep() const = 0;
    //! \return the Newton solver
    virtual NewtonSolver &getSolver() = 0;
    /*!
     * \brief solve the non linear problem over the given time step
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void solve(const real, const real) = 0;
    /*!
     * \brief add a new behaviour integrator
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material id
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    virtual void addBehaviourIntegrator(const std::string &,
                                        const size_type,
                                        const std::string &,
                                        const std::string &) = 0;
    /*!
     * \return the list of material identifiers for which a behaviour
     * integrator has been defined.
     */
    virtual std::vector<size_type> getMaterialIdentifiers() const = 0;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual const Material &getMaterial(const size_type) const = 0;
    /*!
     * \return the material with the given id
     * \param[in] m: material id
     */
    virtual Material &getMaterial(const size_type) = 0;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in] m: material id
     */
    virtual const BehaviourIntegrator &getBehaviourIntegrator(
        const size_type) const = 0;
    /*!
     * \return the behaviour integrator with the given material id
     * \param[in] m: material id
     */
    virtual BehaviourIntegrator &getBehaviourIntegrator(const size_type) = 0;
    /*!
     * \brief add a Dirichlet boundary condition
     * \param[in] bc: boundary condition
     */
    virtual void addBoundaryCondition(
        std::unique_ptr<DirichletBoundaryCondition>) = 0;
    /*!
     * \brief add a new post-processing
     * \param[in] p: post-processing
     */
    virtual void addPostProcessing(
        const std::function<void(const real, const real)>&) = 0;
    /*!
     * \brief execute the registred postprocessings
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     */
    virtual void executePostProcessings(const real, const real) = 0;
    //! \brief revert the state to the beginning of the time step.
    virtual void revert() = 0;
    //! \brief update the state to the end of the time step.
    virtual void update() = 0;
    //! \brief destructor
    virtual ~AbstractNonLinearEvolutionProblem();
  };  // end of struct NonLinearEvolutionProblem

}  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_ABSTRACTNONLINEAREVOLUTIONPROBLEM_HXX */
