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
    //! \brief compute the tangent prediction operator, no integration
    PREDICTION_TANGENT_OPERATOR = -3,
    //! \brief compute the secant prediction operator, no integration
    PREDICTION_SECANT_OPERATOR = -2,
    //! \brief compute the elastic prediction operator, no integration
    PREDICTION_ELASTIC_OPERATOR = -1,
    //! \brief integrate the behaviour without computing a tangent operator
    INTEGRATION_NO_TANGENT_OPERATOR = 0,
    //! \brief integrate the behaviour and compute the elastic operator
    INTEGRATION_ELASTIC_OPERATOR = 1,
    //! \brief integrate the behaviour and compute the secant operator
    INTEGRATION_SECANT_OPERATOR = 2,
    //! \brief integrate the behaviour and compute the tangent operator
    INTEGRATION_TANGENT_OPERATOR = 3,
    //! \brief integrate and compute the consistent tangent operator
    INTEGRATION_CONSISTENT_TANGENT_OPERATOR = 4
  };  // end of enum IntegrationType

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_INTEGRATIONTYPE_HXX */
