/*!
 * \file   include/MFEMMGIS/EnergyPostProcessings.hxx
 * \brief  This file declares the post-processings exporting the stored and
 * dissipated energies
 * \author Thomas Helfer
 * \date   14/12/2021
 */

#ifndef LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_HXX
#define LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_HXX

#include <string>
#include <vector>
#include <fstream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/AbstractNonLinearEvolutionProblemPostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing which exports the stored or dissipated energies of
   * a set of materials in a file.
   */
  template <bool parallel>
  struct EnergyPostProcessingBase
      : public AbstractNonLinearEvolutionProblemPostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     * \param[in] etype: type of energy post-processed
     */
    EnergyPostProcessingBase(
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const Parameters& params,
        const std::string_view etype);
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
    ~EnergyPostProcessingBase() override;

   protected:
    /*!
     * \brief compute the energies of the materials
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \return the energies of the materials, empty on failure
     */
    [[nodiscard]] virtual std::optional<std::vector<real>> computeEnergies(
        Context& ctx,
        const AbstractNonLinearEvolutionProblem& p) const noexcept = 0;
    //! \brief materials
    std::vector<size_type> materials_identifiers;

   private:
    /*!
     * \brief open the output file and write the header
     * \param[in] f: file name
     * \param[in] etype: type of energy post-processed
     */
    void openFile(const std::string& f, const std::string_view etype);
    /*!
     * \brief write the energies of the materials in the output file
     * \param[in] energies: energies of the materials
     */
    void writeResults(const std::vector<real>& energies);
    //! \brief output file
    std::ofstream out;
  };  // end of struct EnergyPostProcessingBase

  /*!
   * \brief a post-processing which exports the stored energies in a file.
   */
  template <bool parallel>
  struct StoredEnergyPostProcessing final : EnergyPostProcessingBase<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    StoredEnergyPostProcessing(
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const Parameters& params);
    //! \brief destructor
    ~StoredEnergyPostProcessing() override;

   private:
    [[nodiscard]] std::optional<std::vector<real>> computeEnergies(
        Context& ctx,
        const AbstractNonLinearEvolutionProblem& p) const noexcept override;
  };  // end of struct StoredEnergyPostProcessing

  /*!
   * \brief a post-processing which exports the dissipated energies in a file.
   */
  template <bool parallel>
  struct DissipatedEnergyPostProcessing final
      : EnergyPostProcessingBase<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    DissipatedEnergyPostProcessing(
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const Parameters& params);
    //! \brief destructor
    ~DissipatedEnergyPostProcessing() override;

   private:
    [[nodiscard]] std::optional<std::vector<real>> computeEnergies(
        Context& ctx,
        const AbstractNonLinearEvolutionProblem& p) const noexcept override;
  };  // end of struct DissipatedEnergyPostProcessing

}  // end of namespace mfem_mgis

#include "MFEMMGIS/EnergyPostProcessings.ixx"

#endif /* LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_HXX */
