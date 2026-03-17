/*!
 * \file   MFEMMGIS/CouplingSchemeConvergenceCriterionBase.hxx
 * \brief  This file declares the `CouplingSchemeConvergenceCriterionBase` class
 * \date   05/12/2022
 */

#ifndef LIB_MFEM_MGIS_COUPLING_SCHEME_CONVERGENCE_CRITERION_BASE_HXX
#define LIB_MFEM_MGIS_COUPLING_SCHEME_CONVERGENCE_CRITERION_BASE_HXX

#include <vector>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractCouplingSchemeConvergenceCriterion.hxx"

namespace mfem_mgis {

  //! \brief a base class for most coupling schemes
  struct MFEM_MGIS_EXPORT CouplingSchemeConvergenceCriterionBase
      : AbstractCouplingSchemeConvergenceCriterion {
    //! \brief constructor
    CouplingSchemeConvergenceCriterionBase();
    //! \brief destructor
    ~CouplingSchemeConvergenceCriterionBase() noexcept override;
  };  // end of CouplingSchemeConvergenceCriterionBase

}  // namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_COUPLING_SCHEME_CONVERGENCE_CRITERION_BASE_HXX */
