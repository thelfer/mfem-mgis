/*!
 * \file   MFEMMGIS/IntegrationType.hxx
 * \brief
 * \author Thomas Helfer
 * \date   31/03/2021
 */

#ifndef LIB_MFEM_MGIS_INTEGRATIONTYPE_HXX
#define LIB_MFEM_MGIS_INTEGRATIONTYPE_HXX

namespace mfem_mgis {

  //! \brief type of integration to be performed
  enum struct IntegrationType {
    PREDICTION_TANGENT_OPERATOR = -3,
    PREDICTION_SECANT_OPERATOR = -2,
    PREDICTION_ELASTIC_OPERATOR = -1,
    INTEGRATION_NO_TANGENT_OPERATOR = 0,
    INTEGRATION_ELASTIC_OPERATOR = 1,
    INTEGRATION_SECANT_OPERATOR = 2,
    INTEGRATION_TANGENT_OPERATOR = 3,
    INTEGRATION_CONSISTENT_TANGENT_OPERATOR = 4
  };  // end of enum IntegrationType

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_INTEGRATIONTYPE_HXX */
