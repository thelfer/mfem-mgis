/*!
 * \file   include/MFEMMGIS/Material.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_MATERIAL_HXX
#define LIB_MFEM_MGIS_MATERIAL_HXX

#include <span>
#include <array>
#include <vector>
#include <memory>
#include <optional>
#include <string_view>
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/RotationMatrix.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"

namespace mfem_mgis {

  // forward declarations
  struct PartialQuadratureSpace;
  struct BehaviourIntegrator;

  /*!
   * \brief a simple structure describing a material an the associated data.
   */
  struct MFEM_MGIS_EXPORT Material : mgis::behaviour::MaterialDataManager {
    /*!
     * \brief a small enumeration
     */
    enum StateSelection {
      BEGINNING_OF_TIME_STEP,
      END_OF_TIME_STEP,
    };
    /*!
     * \brief constructor
     * \param[in] s: quadrature space
     * \param[in] b_ptr: behaviour
     */
    Material(std::shared_ptr<const PartialQuadratureSpace>,
             std::unique_ptr<const Behaviour>);
    /*!
     * \brief set the macroscropic gradients
     * \param[in] g: macroscopic gradients
     */
    void setMacroscopicGradients(std::span<const real>);
    /*!
     * \brief set the rotation matrix
     * \param[in] r: rotation matrix
     * \note this call is only meaningfull in 2D hypotheses for orthotropic
     * behaviours.
     */
    void setRotationMatrix(const RotationMatrix2D &);
    /*!
     * \brief set the rotation matrix
     * \param[in] r: rotation matrix
     * \note this call is only meaningfull in 3D for orthotropic behaviours
     */
    void setRotationMatrix(const RotationMatrix3D &);
    //! \return the quadrature space
    const PartialQuadratureSpace &getPartialQuadratureSpace() const;
    //! \return the quadrature space
    std::shared_ptr<const PartialQuadratureSpace>
    getPartialQuadratureSpacePointer() const;
    /*!
     * \return the rotation matrix for the given integration point
     * \param[in] o: offset of the integration point
     * \note this method is only valid for orthotropic behaviours
     */
    std::array<real, 9u> getRotationMatrixAtIntegrationPoint(const size_type) const;
    //! \brief destructor
    ~Material();

   protected:
    /*!
     * \brief underlying quadrature space
     */
    const std::shared_ptr<const PartialQuadratureSpace> quadrature_space;
    /*!
     * \brief macroscopic gradients
     */
    std::vector<real> macroscopic_gradients;

   protected:
    //! \brief the rotation matrix in 3D
    RotationMatrix2D r2D;
    //! \brief the rotation matrix in 3D
    RotationMatrix3D r3D;
    //! \brief pointer to a function returning the rotation matrix
    std::array<real, 9u> (*get_rotation_fct_ptr)(const RotationMatrix2D &,
                                                 const RotationMatrix3D &,
                                                 const size_type);

   private:
    //! \brief copy constructor (disabled)
    Material(const Material &) = delete;
    //! \brief move constructor (disabled)
    Material(Material &&) = delete;
    //! \brief standard assignement (disabled)
    Material &operator=(const Material &) = delete;
    //! \brief move assignement (disabled)
    Material &operator=(Material &&) = delete;
    /*!
     * \brief underlying behaviour. Only stored for memory management.
     * \note The behaviour can be accessed through the `b` member which is
     * inherited from the `mgis::behaviour::MaterialDataManager` class.
     */
    const std::unique_ptr<const Behaviour> behaviour_ptr;

  };  // end of struct Material

  /*!
   * \brief rotate the thermodynamic forces in the global frame
   * \param[out] f: quadrature function containing the thermodynamic forces in the global frame
   * \param[in] m: material
   * \param[in] s: state considered
   */
  [[nodiscard]] bool rotateThermodynamicsForces(
      PartialQuadratureFunction &,
      Material &,
      const Material::StateSelection = Material::END_OF_TIME_STEP);

  /*!
   * \return a partial quadrature function for the given gradient
   * \param[in] m: material
   * \param[in] n: name of the gradient
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT PartialQuadratureFunction
  getGradient(Material &,
              const std::string_view,
              const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function for the given gradient
   * \param[in] m: material
   * \param[in] n: name of the gradient
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT ImmutablePartialQuadratureFunctionView
  getGradient(const Material &,
              const std::string_view,
              const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function for the given thermodynamic force
   * \param[in] m: material
   * \param[in] n: name of the thermodynamic force
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT PartialQuadratureFunction getThermodynamicForce(
      Material &,
      const std::string_view,
      const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function for the given thermodynamic force
   * \param[in] m: material
   * \param[in] n: name of the thermodynamic force
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT ImmutablePartialQuadratureFunctionView getThermodynamicForce(
      const Material &,
      const std::string_view,
      const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function for the given state variable
   * \param[in] m: material
   * \param[in] n: name of the state variable
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT PartialQuadratureFunction getInternalStateVariable(
      Material &,
      const std::string_view,
      const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function for the given state variable
   * \param[in] m: material
   * \param[in] n: name of the state variable
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT ImmutablePartialQuadratureFunctionView
  getInternalStateVariable(
      const Material &,
      const std::string_view,
      const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function holding the stored energy
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction> getStoredEnergy(
      Material &, const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function holding the stored energy
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT std::optional<ImmutablePartialQuadratureFunctionView>
  getStoredEnergy(const Material &,
                  const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function holding the dissipated energy
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT std::optional<PartialQuadratureFunction> getDissipatedEnergy(
      Material &, const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return a partial quadrature function holding the dissipated energy
   * \param[in] m: material
   * \param[in] s: state considered
   */
  MFEM_MGIS_EXPORT std::optional<ImmutablePartialQuadratureFunctionView>
  getDissipatedEnergy(
      const Material &,
      const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return the stored energy by the whole material if the behaviour computes
   * it, zero otherwise
   * \param[in] bi: behaviour integrator
   * \param[in] s: selection of the state
   */
  MFEM_MGIS_EXPORT real computeStoredEnergy(
      const BehaviourIntegrator &,
      const Material::StateSelection = Material::END_OF_TIME_STEP);
  /*!
   * \return the stored energy if the behaviour computes it, zero otherwise
   * \param[in] bi: behaviour integrator
   * \param[in] s: selection of the state
   */
  MFEM_MGIS_EXPORT real computeDissipatedEnergy(
      const BehaviourIntegrator &,
      const Material::StateSelection = Material::END_OF_TIME_STEP);

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Material.ixx"

#endif /* LIB_MFEM_MGIS_MATERIAL_HXX */
