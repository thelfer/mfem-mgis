/*!
 * \file   include/MFEMMGIS/MFEMForward.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   11/06/2020
 */

#ifndef LIB_MFEM_MGIS_MFEM_FORWARD_HXX
#define LIB_MFEM_MGIS_MFEM_FORWARD_HXX

namespace mfem {

  class Vector;
  class DenseMatrix;
  class Mesh;
  class FiniteElementSpace;
  class FiniteElementCollection;
  class FiniteElement;
  class ElementTransformation;
  class IntegrationRule;
#ifdef MFEM_USE_MPI
  class ParMesh;
  class ParFiniteElementSpace;
#endif

}  // end of namespace mfem

#endif /* LIB_MFEM_MGIS_MFEM_FORWARD_HXX */
