/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.hxx
 * \brief
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PostProcessing.hxx"

namespace mfem_mgis {

  /*!
   * \brief Utilitary classes and functions to post-process MFEM Gridfunctions
   */
  namespace {
    class DiagCoefficient : public mfem::Coefficient
    {
    protected:
      //OLD: Coefficient &lambda, &mu;
      mfem_mgis::Material *mat;
      mfem::GridFunction *u; // displacement
      int si, sj; // component to evaluate, 0 <= si,sj < dim
      
    public:
      DiagCoefficient(mfem_mgis::Material* mat_)
	: mat(mat_), u(NULL), si(0), sj(0) { }
      
      void SetDisplacement(mfem::GridFunction &u_) { u = &u_; }
      
      virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) = 0;
    };
    
    class StressCoefficient : public DiagCoefficient
    {
    public:
      using DiagCoefficient::DiagCoefficient;
      void SetComponent(int i, int j) { si = i; sj = j; }
      double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) {
	MFEM_ASSERT(u != NULL, "displacement field is not set");
	return 0.;
	//TODO:      double L = lambda.Eval(T, ip);
	//TODO:      double M = mu.Eval(T, ip);
	//TODO:      u->GetVectorGradient(T, grad);
	//TODO:      if (si == sj)
	//TODO:	{
	//TODO:	  double div_u = grad.Trace();
	//TODO:	  return L*div_u + 2*M*grad(si,si);
	//TODO:	}
	//TODO:      else
	//TODO:	{
	//TODO:	  return M*(grad(si,sj) + grad(sj,si));
	//TODO:	}
      }
    };
    
    class StrainCoefficient : public DiagCoefficient
    {
    public:
      using DiagCoefficient::DiagCoefficient;
      void SetComponent(int i, int j) { si = i; sj = j; }
      double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip)
      {
	MFEM_ASSERT(u != NULL, "displacement field is not set");
	return 0.;
	//TODO:      u->GetVectorGradient(T, grad);
	//TODO:      return (grad(si,sj)+grad(sj,si));
      }
    };
  }

  
  
  /*!
   * \brief a post-processing to export the results to paraview
   */
  template <bool parallel>
  struct ParaviewExportResults final : public PostProcessing<parallel> {
    /*!
     * \brief constructor
     * \param[in] p: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportResults(NonLinearEvolutionProblemImplementation<parallel>&,
                          const Parameters&);
    //
    void execute(NonLinearEvolutionProblemImplementation<parallel>&,
                 const real,
                 const real) override;
    //! \brief destructor
    ~ParaviewExportResults() override;

   private:
    const int dim, sdim;
    //! \brief Paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief Exported displacement grid function
    mfem_mgis::GridFunction<parallel> displacement;
    //!  \brief List of materials
    std::vector<mfem_mgis::Material*> mgis_materials;
    //!  \brief Number of materials
    int nb_materials;
    //!  \brief Temporary data structure for stress calculation
    std::vector<StressCoefficient*> stress_c;
    //!  \brief Temporary data structure for strain calculation
    std::vector<StrainCoefficient*> strain_c;
    //!  \brief Exported local stress grid functions
    std::vector<mfem_mgis::GridFunction<parallel>*> stress;
    //! \brief Exported local strain grid functions
    std::vector<mfem_mgis::GridFunction<parallel>*> strain;
    //!
    size_type cycle;
  };  // end of struct ParaviewExportResults

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportResults.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX */
