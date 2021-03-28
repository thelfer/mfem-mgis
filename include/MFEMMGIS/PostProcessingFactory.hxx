/*!
 * \file   include/MFEMMGIS/PostProcessingFactory.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEM_MGIS_POSTPROCESSINGFACTORY_HXX
#define LIB_MFEM_MGIS_POSTPROCESSINGFACTORY_HXX

#include <map>
#include <memory>
#include <functional>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;
  // forward declaration
  template <bool parallel>
  struct NonLinearEvolutionProblemImplementation;
  // forward declaration
  template <bool parallel>
  struct PostProcessing;

  /*!
   * \brief an abstract factory for behaviour integrators
   * \tparam parallel: boolean stating if parallel post-processing are
   * considered
   */
  template <bool parallel>
  struct PostProcessingFactory;

#ifdef MFEM_USE_MPI

  //! \brief partial specialisation in parallel
  template <>
  struct MFEM_MGIS_EXPORT PostProcessingFactory<true> {
    //! a simple alias
    using Generator = std::function<std::unique_ptr<PostProcessing<true>>(
        NonLinearEvolutionProblemImplementation<true>&, const Parameters&)>;
    //! \return the unique instance of the class
    static PostProcessingFactory& getFactory();
    /*!
     * \brief register a new post-processing
     * \param[in] n: name of the post-processing
     * \param[in] g: generator of the post-processing
     */
    void add(std::string_view, Generator);
    /*!
     * \return the requested post-processing
     * \param[in] n: name of the post-processing
     * \param[in] p: non linear evolution postprocessing
     * \param[in] params: parameters passed to the post-processing
     */
    std::unique_ptr<PostProcessing<true>> generate(
        std::string_view,
        NonLinearEvolutionProblemImplementation<true>&,
        const Parameters&) const;

   private:
    //! \brief default destructor
    PostProcessingFactory();
    //! \brief destructor
    ~PostProcessingFactory();
    //! \brief registred factories
    std::map<std::string, Generator, std::less<>> generators;
  };  // end of struct PostProcessingFactory

#endif /* MFEM_USE_MPI */

  //! \brief partial specialisation in sequential
  template <>
  struct MFEM_MGIS_EXPORT PostProcessingFactory<false> {
    //! a simple alias
    using Generator = std::function<std::unique_ptr<PostProcessing<false>>(
        NonLinearEvolutionProblemImplementation<false>&, const Parameters&)>;
    //! \return the unique instance of the class
    static PostProcessingFactory& getFactory();
    /*!
     * \brief register a new post-processing
     * \param[in] n: name of the post-processing
     * \param[in] g: generator of the post-processing
     */
    void add(std::string_view, Generator);
    /*!
     * \return the requested post-processing
     * \param[in] n: name of the post-processing
     * \param[in] p: non linear evolution postprocessing
     * \param[in] params: parameters passed to the post-processing
     */
    std::unique_ptr<PostProcessing<false>> generate(
        std::string_view,
        NonLinearEvolutionProblemImplementation<false>&,
        const Parameters&) const;

   private:
    //! \brief default destructor
    PostProcessingFactory();
    //! \brief destructor
    ~PostProcessingFactory();
    //! \brief registred factories
    std::map<std::string, Generator, std::less<>> generators;
  };  // end of struct PostProcessingFactory

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_POSTPROCESSINGFACTORY_HXX */
