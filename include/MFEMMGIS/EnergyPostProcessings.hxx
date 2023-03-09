/*!
 * \file   include/MFEMMGIS/StoredEnergyPostProcessing.hxx
 * \brief
 * \author Thomas Helfer
 * \date   14/12/2021
 */

#ifndef LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_HXX
#define LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_HXX

#include <string>
#include <vector>
#include <fstream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing which export the stored or dissipated energies of
   * a set of materials in a file.
   */
  template <bool parallel>
  struct EnergyPostProcessingBase : public PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     * \param[in] etype: type of energy post-processed
     */
    EnergyPostProcessingBase(
        NonLinearEvolutionProblemImplementation<parallel> &,
        const Parameters &,
        const std::string_view);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel> &,
                 const real,
                 const real) override;
    //! \brief destructor
    ~EnergyPostProcessingBase() override;

   protected:
    //!
    virtual std::vector<real> computeEnergies(
        const AbstractNonLinearEvolutionProblem &) const = 0;
    //! \brief materials
    std::vector<size_type> materials_identifiers;

   private:
    /*!
     * \brief open the output file and write the header
     * \param[in] f: file name
     * \param[in] etype: type of energy post-processed
     */
    void openFile(const std::string &, const std::string_view);
    //!
    void writeResults(const std::vector<real> &);
    //! \brief output file
    std::ofstream out;
  };  // end of struct EnergyPostProcessingBase

  /*!
   * \brief a post-processing which exported the stored energies in a file.
   */
  template <bool parallel>
  struct StoredEnergyPostProcessing final : EnergyPostProcessingBase<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    StoredEnergyPostProcessing(
        NonLinearEvolutionProblemImplementation<parallel> &,
        const Parameters &);
    //! \brief destructor
    ~StoredEnergyPostProcessing() override;

   private:
    std::vector<real> computeEnergies(
        const AbstractNonLinearEvolutionProblem &) const override;
  };  // end of struct StoredEnergyPostProcessing

  /*!
   * \brief a post-processing which exported the stored energies in a file.
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
        NonLinearEvolutionProblemImplementation<parallel> &,
        const Parameters &);
    //! \brief destructor
    ~DissipatedEnergyPostProcessing() override;

   private:
    std::vector<real> computeEnergies(
        const AbstractNonLinearEvolutionProblem &) const override;
  };  // end of struct DissipatedEnergyPostProcessing

}  // end of namespace mfem_mgis

#include "MFEMMGIS/EnergyPostProcessings.ixx"

#endif /* LIB_MFEM_MGIS_ENERGYPOSTPROCESSINGS_HXX */
