/*!
 * \file   PartialQuadratureFunctionUtilities.cxx
 * \brief
 * \author Thomas Helfer
 * \date   30/04/2025
 */

#include "MFEMMGIS/Algorithms.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionEvaluator.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionUtilities.hxx"

namespace mfem_mgis {

  bool rotateThermodynamicsForcesToGlobalFrame(
      Context &ctx,
      PartialQuadratureFunction &f,
      const Material &m,
      const Material::StateSelection s) {
    RotatedThermodynamicForcesMatrixPartialQuadratureFunctionEvalutor e(m, s);
    return assign(ctx, f, e);
  }  // end of rotateThermodynamicsForces

}  // end of namespace mfem_mgis
