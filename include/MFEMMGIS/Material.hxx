/*!
 * \file   include/MFEMMGIS/Material.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_MATERIAL_HXX
#define LIB_MFEM_MGIS_MATERIAL_HXX

#include <array>
#include <vector>
#include <memory>
#include "MGIS/Span.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"
#include "MFEMMGIS/RotationMatrix.hxx"

namespace mfem_mgis {

  /*!
   * \brief PartialQuadratureSpace
   */
  struct PartialQuadratureSpace;

  /*!
   * \brief a simple structure describing a material an the associated data.
   */
  struct MFEM_MGIS_EXPORT Material : mgis::behaviour::MaterialDataManager {
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
    void setMacroscopicGradients(mgis::span<const real>);
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

}  // end of namespace mfem_mgis

#include "MFEMMGIS/Material.ixx"

#endif /* LIB_MFEM_MGIS_MATERIAL_HXX */
