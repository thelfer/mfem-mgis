#include <memory>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "mpi.h"
#include "mfem/general/optparser.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/BehaviourIntegrator.hxx"
#include "MFEMMGIS/PostProcessing.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/MFEMForward.hxx"


struct MeanStressPostProcessing
  : public mfem_mgis::PostProcessing<true>// false=séquenciel true=parallèle
{

  MeanStressPostProcessing(const mfem_mgis::NonLinearEvolutionProblem<true>& p,
			   const std::string& f)
  {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0) {
          this->out.open(f);
          out << "# first column: time\n";
          auto c = mfem_mgis::size_type{2};
          for (const auto &mi : p.getMaterialIdentifiers()) {
            const auto &bi = p.getBehaviourIntegrator(mi);
            const auto &s1 = bi.getMaterial().s1;
            const auto thsize = static_cast<mfem_mgis::size_type>(
                s1.thermodynamic_forces_stride);
            for (mfem_mgis::size_type k = 0; k != thsize; ++k, ++c) {
              out << "# " << c << " column: " << k
                  << " component of the mean thermodynamic force for material "
                  << mi << '\n';
            }
          }
        }
  }
  
  void execute(mfem_mgis::NonLinearEvolutionProblem<true>& p, // *=pointeur, &=référence
	       const mfem_mgis::real t,
	       const mfem_mgis::real dt) override 
  {
    const auto nmax = p.getFiniteElementSpace().GetMesh()->attributes.Max() + 1;
    std::vector<std::vector<mfem_mgis::real>> stress_integrals(nmax);
    const auto& mis = p.getMaterialIdentifiers();
    for (const auto& mi : mis)
    {
      const auto& bi = p.getBehaviourIntegrator(mi);
      const auto& s1 = bi.getMaterial().s1;
      const auto thsize = s1.thermodynamic_forces_stride;
      stress_integrals[mi].resize(thsize, mfem_mgis::real(0));
    }
    std::vector<mfem_mgis::real> volume(stress_integrals.size(),
					mfem_mgis::real(0));
    const auto& fes = p.getFiniteElementSpace();

    for (mfem_mgis::size_type i = 0; i < fes.GetNE(); i++)
    {
      auto& e = *(fes.GetFE(i));
      auto& tr  = *(fes.GetElementTransformation(i));
      const auto m = tr.Attribute;
      const auto& bi = p.getBehaviourIntegrator(m);      
      const auto& ir = bi.getIntegrationRule(e, tr);
      const auto& s1 = bi.getMaterial().s1;
      const auto thsize =
	static_cast<mfem_mgis::size_type>(s1.thermodynamic_forces_stride);
//static_cast effectue une conversion de type explicite. 
      auto& s = stress_integrals[m];
      // offset associé à l'élément courant
      const auto eoffset = bi.getPartialQuadratureSpace().getOffset(tr.ElementNo);
      for (mfem_mgis::size_type j = 0; j < ir.GetNPoints(); j++)
      {
	const auto o = eoffset + j;
	const auto &ip = ir.IntPoint(j);
	tr.SetIntPoint(&ip);
	const auto thf = s1.thermodynamic_forces.subspan(o * thsize, thsize);
	const auto w = ip.weight * tr.Weight();
	for (mfem_mgis::size_type k = 0; k != thsize; ++k)
        {
	  s[k] += w * thf[k];
	//std::cerr << "s[k] = " << s[k] << std::endl;
	  
//changer expression pour passer au reduce pour multiplication et division
	}
        volume[m] += w;
	
      }
    }
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        if (rank == 0) {
          this->out << t + dt;
        }
    for (const auto& mi : p.getMaterialIdentifiers())
    {
	{
          
	  auto &s = stress_integrals[mi];
	  double res_mpi_v = 0;
	  const auto& bi = p.getBehaviourIntegrator(mi); 
	  const auto &s1 = bi.getMaterial().s1;
	  //const auto thsize = s1.thermodynamic_forces_stride;
	  const auto thsize =static_cast<mfem_mgis::size_type>(s1.thermodynamic_forces_stride);
          std::vector<double> stress_integral(s.size(), 0);
          MPI_Reduce(&volume[mi], &res_mpi_v, 1, MPI_DOUBLE, MPI_SUM, 0,
                     MPI_COMM_WORLD);
		std::cout << " MPI_V " << res_mpi_v << std::endl;
          MPI_Reduce(s.data(), stress_integral.data(), s.size(), MPI_DOUBLE,
                     MPI_SUM, 0, MPI_COMM_WORLD);
          if (rank == 0) {
            for (mfem_mgis::size_type k = 0; k != thsize; ++k) {
              this->out << " " << stress_integral[k] / res_mpi_v;
            }
          }
        }   
	//std::cerr << "volume[m] = " << volume[m] << std::endl;
        if (rank == 0) {
          this->out << '\n';
        }
    }
  }
  ~MeanStressPostProcessing() override;
//Destructeur de MeanStressPostProcessing
//override = default, pour le constructeur informer usage par défault 
//Toujours utiliser override sur les méthodes qui redéfinissent une méthode d’une classe mère.

protected:
  std::ofstream out;
  
//Destructeur de MeanStressPostProcessing
//override = default, pour le constructeur informer usage par défault 
//Toujours utiliser override sur les méthodes qui redéfinissent une méthode d’une classe mère.
};

MeanStressPostProcessing::~MeanStressPostProcessing() = default; // demander détails thomas

template <bool parallel>
void exit_on_failure (){
  if constexpr (parallel) MPI_Finalize();
  std::exit(EXIT_FAILURE);
}


int main(int argc, char *argv[]) 
{
	MPI_Init(&argc,&argv);	
	int code;
	//int var_rank;
	MPI_Comm_size(MPI_COMM_WORLD,&code);
	MPI_Comm_rank(MPI_COMM_WORLD,&code);
	constexpr const auto dim = mfem_mgis::size_type{3}; // déf en 3D
	const char* mesh_file = "cube.mesh"; //def nom du mesh et de la forme 
	const char* behaviour = "Plasticity"; // déformation plastique 
	const char* library = "src/libBehaviour.so"; // chercher la lib
	const auto order = 1; // déf de l'ordre en 1
	// loading the mesh
	std::shared_ptr<mfem_mgis::Mesh<true>> mesh;
	auto smesh = std::make_shared<mfem::Mesh>(mesh_file, 1, 1);
	if (dim != smesh->Dimension()) {
	  std::cerr << "Invalid mesh dimension \n";
	  exit_on_failure<true>();
	}
	mesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, *smesh);
	//	auto mesh = std::make_shared<mfem_mgis::Mesh<true>>(mesh_file, 1, 1);//make_shared poiteur pour une alloc
	if (mesh->Dimension() != dim) 
	{//To access members of a structure through a pointer, use "->"
	  exit_on_failure<true>();
	}
	// building the non linear problem
	//	std::cerr << __LINE__ << std::endl;
	auto feh = std::make_shared<mfem::H1_FECollection>(order,dim);
	auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(mesh, feh, 3);
	mfem_mgis::NonLinearEvolutionProblem<true> problem(fed,
							   mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
	//FEDmesh,FECollection,FESpace)
  	problem.addBehaviourIntegrator("Mechanics", 1, library, behaviour);
	//1 numero d'un mat qui compose l'ensemble
	auto& m1 = problem.getMaterial(1);
	mgis::behaviour::setExternalStateVariable(m1.s0, "Temperature", 293.15);
	//s0 etat début de pas s1 fin de pas
	mgis::behaviour::setExternalStateVariable(m1.s1, "Temperature", 293.15);
	auto fixed_dirichlet_dofs = mfem::Array<mfem_mgis::size_type>{};
  	auto append_dof = [&](const auto boundary, const auto normal) 
	{
    		auto tmp = mfem::Array<mfem_mgis::size_type>{};
    		auto boundaries_markers = mfem::Array<mfem_mgis::size_type>(mesh->bdr_attributes.Max());
    		boundaries_markers = 0;
    		boundaries_markers[boundary] = 1;
    		problem.getFiniteElementSpace().GetEssentialTrueDofs(boundaries_markers,tmp, normal);
    		fixed_dirichlet_dofs.Append(tmp);
    		return tmp;
  	};
	append_dof(0, 1);                     // xz0
	append_dof(1, 2);                     // xy0
	append_dof(4, 0);                     // yz0
	auto yz1_ux_dofs = append_dof(2, 0);  // yz1
	problem.SetEssentialTrueDofs(fixed_dirichlet_dofs);
	// solving the problem
#ifdef MFEM_USE_SUITESPARSE
	mfem::UMFPackSolver lsolver;
#else
	mfem::CGSolver lsolver(MPI_COMM_WORLD);
	lsolver.SetRelTol(1e-12);
	lsolver.SetMaxIter(300);
	lsolver.SetPrintLevel(1);
#endif
	auto& solver = problem.getSolver();
	solver.SetSolver(lsolver);//AFFECTE UN POINTEUR
	solver.SetPrintLevel(0);
	solver.SetRelTol(1e-12);
	solver.SetAbsTol(0);
	solver.SetMaxIter(10);

	const auto u_t = [](const auto t) 
	{
	    	if (t < 0.3) 
		{
			return 3e-2 * t;
	    	}
		else if (t < 0.6) 
		{
			return 0.009 - 0.1 * (t - 0.3);
	    	}
	    	return -0.021 + 0.1 * (t - 0.6);
  	};
// post-processing 1
	problem.addPostProcessing(std::make_unique<MeanStressPostProcessing>(problem,"mean_stresses.txt"));
// post-processing 2
	std::ofstream out("stress-strain.txt");
	problem.addPostProcessing([&problem, &out](const mfem_mgis::real t, const mfem_mgis::real dt)
	{
   		const auto& s1 = problem.getMaterial(1).s1;
        	const auto e = s1.gradients[0];
        	const auto s = s1.thermodynamic_forces[0];
        	out << t + dt << " " << e << " " << s << "\n";
      	}			);
// loop over time step
	const auto nsteps = mfem_mgis::size_type{100};
	auto t = mfem_mgis::real{0};
	const auto dt = mfem_mgis::real{1} / nsteps;
	for (mfem_mgis::size_type i = 0; i != nsteps; ++i) 
	{
// setting the boundary values
		auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
		const auto u = u_t(t + dt);
		for (mfem_mgis::size_type idx = 0; idx != yz1_ux_dofs.Size(); ++idx)
		{
			u1[yz1_ux_dofs[idx]] = u;
		}
// resolution
	problem.solve(t, dt);
	problem.executePostProcessings(t, dt);
	problem.update();
	t += dt;
	}
	MPI_Finalize();
	return EXIT_SUCCESS;	
	MPI_Finalize();
	return 0;
}
