/*!
 * \file   include/MFEMMGIS/MGISIntegrator.hxx
 * \brief
 * \author Thomas Helfer
 * \date   8/06/2020
 */

#ifndef LIB_MFEM_MGIS_MGISINTEGRATOR_HXX
#define LIB_MFEM_MGIS_MGISINTEGRATOR_HXX

#include <memory>
#include <unordered_map>

#include "mfem/fem/nonlininteg.hpp"

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/MGISForward.hxx"
#include "MFEMMGIS/Material.hxx"

namespace mfem_mgis {

  // forward declaration
  struct BehaviourIntegrator;

  /*!
   * \brief Base class for non linear integrators based on an MGIS' behaviours.
   * This class manages an mapping associating a material and its identifier
   */
  struct MFEM_MGIS_EXPORT MGISIntegrator
      : public mfem::NonlinearFormIntegrator {
    //! \brief a simple alias
    using Behaviour = mgis::behaviour::Behaviour;

    /*!
     * \brief constructor
     * \param[in] fs: finite element space
     * \param[in] h: modelling hypothesis
     */
    MGISIntegrator(std::shared_ptr<const mfem::FiniteElementSpace>,
                   const Hypothesis);

    void AssembleElementVector(const mfem::FiniteElement &,
                               mfem::ElementTransformation &,
                               const mfem::Vector &,
                               mfem::Vector &) override;

    void AssembleElementGrad(const mfem::FiniteElement &,
                             mfem::ElementTransformation &,
                             const mfem::Vector &,
                             mfem::DenseMatrix &) override;

    /*!
     * \brief add a new material
     * \param[in] n: name of the behaviour integrator
     * \param[in] m: material id
     * \param[in] l: library name
     * \param[in] b: behaviour name
     */
    virtual void addBehaviourIntegrator(const std::string &,
                                        const size_type,
                                        const std::string &,
                                        const std::string &);

    //! \brief destructor
    ~MGISIntegrator() override;

   private:
    //! \brief underlying finit element space
    const std::shared_ptr<const mfem::FiniteElementSpace> fe_space;
    /*!
     * \brief mapping between the material integrator and the behaviour
     * integrator.
     */
    std::unordered_map<size_type,  // material id
                       std::shared_ptr<BehaviourIntegrator>>
        behaviour_integrators;
    //! \brief modelling hypothesis
    Hypothesis hypothesis;
  };  // end of MGISIntegratorBase

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MGISINTEGRATOR_HXX */
