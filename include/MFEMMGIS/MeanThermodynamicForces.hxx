/*!
 * \file   include/MFEMMGIS/MeanThermodynamicForces.hxx
 * \brief
 * \author Thomas Helfer, Hugo Copin
 * \date   08/04/2021
 */

#ifndef LIB_MFEM_MGIS_MEANTHERMODYNAMICFORCES_HXX
#define LIB_MFEM_MGIS_MEANTHERMODYNAMICFORCES_HXX

#include <string>
#include <vector>
#include <fstream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing which computes the mean values of each component
   * of the thermodynamic forces and prints them in a file.
   */
  template <bool parallel>
  struct MeanThermodynamicForces final
      : public AbstractNonLinearEvolutionProblemPostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    MeanThermodynamicForces(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const Parameters& params);
    /*!
     * \brief do nothing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessing(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const real t) noexcept override;
    /*!
     * \brief execute the post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool execute(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const real t,
        const real dt) noexcept override;
    //! \brief destructor
    ~MeanThermodynamicForces() override;

   private:
    /*!
     * \brief open the output file and write the header
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear problem
     * \param[in] f: file name
     * \return true on success
     */
    [[nodiscard]] bool openFile(
        Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const std::string& f) noexcept;
    /*!
     * \brief write the mean of the value of the thermodynamic forces of a
     * material to the output file.
     * \param[in] tf_integral: integral of the thermodynamic forces over the
     * material.
     * \param[in] v: volume of the material
     */
    void writeResults(const std::vector<real>& tf_integral, const real v);

    //! \brief output file
    std::ofstream out;
  };  // end of struct MeanThermodynamicForces

}  // end of namespace mfem_mgis

#include "MFEMMGIS/MeanThermodynamicForces.ixx"

#endif /* LIB_MFEM_MGIS_MEANTHERMODYNAMICFORCES_HXX */
