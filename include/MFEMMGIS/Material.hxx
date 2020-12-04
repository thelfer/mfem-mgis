/*!
 * \file   include/MFEMMGIS/Material.hxx
 * \brief
 * \author Thomas Helfer
 * \date   26/08/2020
 */

#ifndef LIB_MFEM_MGIS_MATERIAL_HXX
#define LIB_MFEM_MGIS_MATERIAL_HXX

#include <memory>
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Behaviour.hxx"

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

    //! \brief destructor
    ~Material();

   protected:

    /*!
     * \brief underlying quadrature space
     */
    const std::shared_ptr<const PartialQuadratureSpace> quadrature_space;

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

};  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_MATERIAL_HXX */
