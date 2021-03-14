/*!
 * \file   tests/PeriodicTest.cxx
 * \brief
 * This example code solves a simple linear elasticity problem
 * describing a multi-material square.
 * This problem has a 1D analytic solution along x1 dimension,
 * the solution is constant along x2 dimension which is also periodic.
 *
 * The geometry of the domain is assumed to be as
 * follows:
 *
 *             x2=1  +----------+----------+
 *                   | material | material |
 *                   |    1     |    2     |
 *             x2=0  +----------+----------+
 *                 x1=0                   x1=1
 *
 *
 * Specifically, we approximate the weak form of -div(sigma(u))=0
 * where sigma(u)=lambda*div(u)*I+mu*(grad*u+u*grad) is the stress
 * tensor corresponding to displacement field u, and lambda and mu
 * are the material Lame constants. The boundary conditions are
 * periodic.
 *
 * Mechanical strain:
 *                 eps = E + grad_s v
 *
 *           with  E the given macrocoscopic strain
 *                 v the periodic displacement fluctuation
 * Displacement:
 *                   u = U + v
 *
 *           with  U the given displacement associated to E
 *                   E = grad_s U
 *
 * The local microscopic strain is equal, on average, to the macroscopic strain:
 *           <eps> = <E>
 * \author Thomas Helfer, Guillaume Latu
 * \date   14/10/2020
 */

#include <memory>
#include <cstdlib>
#include <iostream>
#include "mfem/general/optparser.hpp"
#include "mfem/fem/intrules.hpp"
#include "mfem/linalg/solvers.hpp"
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"

#define USE_PROFILER 1
#define LIB_PROFILER_PRINTF MpiPrintf
#include "MFEMMGIS/libProfiler.h"

void (*getSolution(const std::size_t i))(const mfem::Vector&, mfem::Vector&) {
  constexpr const auto xthr = mfem_mgis::real(1) / 2;
  std::array<void (*)(const mfem::Vector&, mfem::Vector&), 6u> solutions = {
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(4) / 30;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(0) = gradx * x(0);
        } else {
          u(0) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(1) = gradx * x(0);
        } else {
          u(1) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector& x, mfem::Vector& u) {
        constexpr const auto gradx = mfem_mgis::real(1) / 3;
        u = mfem_mgis::real{};
        if (x(0) < xthr) {
          u(2) = gradx * x(0);
        } else {
          u(2) = gradx * xthr - gradx * (x(0) - xthr);
        }
      },
      +[](const mfem::Vector&, mfem::Vector& u) { u = mfem_mgis::real{}; }};
  return solutions[i];
}

std::shared_ptr<mfem::Solver> getLinearSolver(const std::size_t i) {
  using generator = std::shared_ptr<mfem::Solver> (*)();
  const generator generators[] = {
      []() -> std::shared_ptr<mfem::Solver> {
        std::shared_ptr<mfem::GMRESSolver> pgmres(nullptr);
	pgmres = std::make_shared<mfem::GMRESSolver>(MPI_COMM_WORLD);
        pgmres->iterative_mode = false;
        pgmres->SetRelTol(1e-9);
        pgmres->SetAbsTol(1e-9);
        pgmres->SetMaxIter(800);
        pgmres->SetPrintLevel(1);
        return pgmres;
      },
      []() -> std::shared_ptr<mfem::Solver> {
        std::shared_ptr<mfem::CGSolver> pcg(nullptr);                                
	pcg = std::make_shared<mfem::CGSolver>(MPI_COMM_WORLD);
        pcg->SetRelTol(1e-9);
        pcg->SetAbsTol(1e-9);
        pcg->SetMaxIter(800);
        pcg->SetPrintLevel(1);
        return pcg;
      }
  };
  return (generators[i])();
}

void setBoundaryConditions(mfem_mgis::NonLinearEvolutionProblemBase<true>& problem) {
  // Impose no displacement on the first node
  // which needs to be on x=xmin or x=xmax axis.
  // ux=0, uy=0, uz=0 on this point.
  auto *pmesh = problem.getFiniteElementSpace().GetMesh();
  const auto dim = pmesh->Dimension();
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  mfem::Array<int> ess_tdof_list;
  ess_tdof_list.SetSize(0);
  {
    mfem::GridFunction nodes(&problem.getFiniteElementSpace());
    int found = 0;
    bool reorder_space = true;
    pmesh->GetNodes(nodes);
    const auto size = nodes.Size()/dim;
    std::cerr << "Number of nodes: " << size << std::endl;

     // Traversal of all dofs to detect which one is (0,0,0)
     for (int i = 0; i < size; ++i) {
       double coord[dim]; // coordinates of a node
       double dist = 0.;
       for (int j = 0; j < dim; ++j) {
	 if (reorder_space)
	   coord[j] = (nodes)[j * size + i];
	 else
	   coord[j] = (nodes)[i * dim + j];
	 // because of periodic BC, 0. is also 1.
	 if (abs(coord[j] - 1.) < 1e-7)
	   coord[j] = 0.;
	 dist += coord[j] * coord[j];
       }
       // If distance is close to zero, we have our reference point
       if (dist < 1.e-16) {
	 std::cout << myid << ":coord: " <<coord[0] << " " << coord[1] << std::endl;
	 for (int j = 0; j < dim; ++j) {
	   int id_unk;
	   if (reorder_space)
	     {
	     //id_unk = (j * size + i);
	     id_unk = problem.getFiniteElementSpace().GetLocalTDofNumber(j * size + i);
	     }
	   else
	     {
	       //id_unk = (i * dim + j);
	       id_unk = problem.getFiniteElementSpace().GetLocalTDofNumber(i * dim + j);
	     }
	   if (id_unk >= 0)
	     {
	       found = 1;
	       ess_tdof_list.Append(id_unk);
	       std::cout <<  myid << ":bloqued unknown: " << id_unk << std::endl;
	     }
	 }
       }
     }
     MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
     MFEM_VERIFY(found, "Reference point at (0,0) was not found");
  }
  problem.SetEssentialTrueDofs(ess_tdof_list);
}

void setSolverParameters(mfem_mgis::NonLinearEvolutionProblemBase<true>& problem,
                         mfem::Solver& lsolver) {
  auto& solver = problem.getSolver();
  solver.iterative_mode = true;
  solver.SetSolver(lsolver);
  solver.SetPrintLevel(0);
  solver.SetRelTol(1e-9);
  solver.SetAbsTol(1e-9);
  solver.SetMaxIter(10);
}  // end of setSolverParmeters

bool checkSolution(mfem_mgis::NonLinearEvolutionProblemBase<true>& problem,
                   const std::size_t i) {
  constexpr const auto eps = mfem_mgis::real{1e-10};
  const auto dim = problem.getFiniteElementSpace().GetMesh()->Dimension();
  // recover the solution as a grid function
  auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
  mfem_mgis::GridFunction<true> x(&problem.getFiniteElementSpace());
  x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
  x.SetFromTrueVector();
  // comparison to analytical solution
  mfem::VectorFunctionCoefficient sol_coef(dim, getSolution(i));
  const auto error = x.ComputeL2Error(sol_coef);
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (error > eps) {
    if (myid == 0)
      std::cerr << "Error is greater than threshold (" << error << " > " << eps
              << ")\n";
    return false;
  } 
  if (myid == 0)
    std::cerr << "Error is lower than threshold (" << error << " < " << eps
	      << ")\n";
  return true;
}

void exportResults(mfem_mgis::NonLinearEvolutionProblemBase<true>& problem,
                   const std::size_t tcase) {
  auto* const mesh = problem.getFiniteElementSpace().GetMesh();
  auto& u1 = problem.getUnknownsAtEndOfTheTimeStep();
  mfem_mgis::GridFunction<true> x(&problem.getFiniteElementSpace());
  x.MakeTRef(&problem.getFiniteElementSpace(), u1, 0);
  x.SetFromTrueVector();

  mfem::ParaViewDataCollection paraview_dc(
      "PeriodicTestOutput-" + std::to_string(tcase), mesh);
  paraview_dc.RegisterField("u", &x);
  paraview_dc.SetCycle(0);
  paraview_dc.SetTime(0.0);
  paraview_dc.Save();
};

struct TestParameters {
  const char* mesh_file = nullptr;
  const char* library = nullptr;
  int algo  = 0;
  int order = 1;
  int tcase = 0;
  int refine = 2;
  int linearsolver = 0;
};

TestParameters parseCommandLineOptions(const int argc, char** const argv) {
  TestParameters p;
  // options treatment
  mfem::OptionsParser args(argc, argv);
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refine, "-r", "--refine",
                 "Level of refinement");
  args.AddOption(&p.algo, "-a", "--algo",
                 "Algorithm to be used");
  args.AddOption(&p.tcase, "-t", "--test-case",
                 "identifier of the case : Exx->0, Eyy->1, Ezz->2, Exy->3, "
                 "Exz->4, Eyz->5");
  args.AddOption(
      &p.linearsolver, "-ls", "--linearsolver",
      "identifier of the linear solver: 0 -> GMRES, 1 -> CG, 2 -> UMFPack");
  args.Parse();
  int myid;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if ((!args.Good()) || (p.mesh_file == nullptr) || (p.library == nullptr)) {
    if (myid == 0)
      args.PrintUsage(std::cout);
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  if (myid == 0)
    args.PrintOptions(std::cout);
  if ((p.tcase < 0) || (p.tcase > 5)) {
    if (myid == 0)
      std::cerr << "Invalid test case\n";
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  return p;
}

void executeMFEMMGISTest(const TestParameters& p) {
  PROFILER_START(2_read_mesh);
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace
  mfem::Mesh* smesh = new mfem::Mesh(p.mesh_file, 1, 1);
  if (dim != smesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  mfem::ParMesh *pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *smesh);
  delete smesh;
  MPI_Barrier(MPI_COMM_WORLD);
  std::shared_ptr<mfem::ParMesh> mesh(pmesh);

  if (dim != mesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  PROFILER_END(); PROFILER_START(3_refine_mesh);

  for (int i = 0 ; i < p.refine ; i++)
    mesh->UniformRefinement();

  PROFILER_END(); PROFILER_START(4_initialize_fem); 
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblem<true> problem(
      std::make_shared<mfem_mgis::FiniteElementDiscretization>(
          mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim), 3),
      mgis::behaviour::Hypothesis::TRIDIMENSIONAL);
  problem.addBehaviourIntegrator("Mechanics", 1, p.library, "Elasticity");
  problem.addBehaviourIntegrator("Mechanics", 2, p.library, "Elasticity");
  // materials
  auto& m1 = problem.getMaterial(1);
  auto& m2 = problem.getMaterial(2);
  // setting the material properties
  auto set_properties = [](auto& m, const double l, const double mu) {
    mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
    mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
  };
  set_properties(m1, 100, 75);
  set_properties(m2, 200, 150);
  //
  auto set_temperature = [](auto& m) {
    mgis::behaviour::setExternalStateVariable(m.s0, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(m.s1, "Temperature", 293.15);
  };
  set_temperature(m1);
  set_temperature(m2);
  // macroscopic strain
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{});
  if (p.tcase < 3) {
    e[p.tcase] = 1;
  } else {
    e[p.tcase] = 1.41421356237309504880 / 2;
  }
  m1.setMacroscopicGradients(e);
  m2.setMacroscopicGradients(e);
  //
  setBoundaryConditions(problem);
  //
  auto lsolver = getLinearSolver(p.linearsolver);
  setSolverParameters(problem, *(lsolver.get()));
  // solving the problem
  PROFILER_END(); PROFILER_START(5_solve);
  problem.solve(1);
  PROFILER_END(); PROFILER_START(6_postprocess);
  //
  if (!checkSolution(problem, p.tcase)) {
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  //
  exportResults(problem, p.tcase);
  PROFILER_END(); 
}

struct ElasticityNonLinearIntegrator final
    : public mfem::NonlinearFormIntegrator {
  ElasticityNonLinearIntegrator(mfem::Coefficient& l,
                                mfem::Coefficient& m,
                                int _tcase) {
    lambda = &l;
    mu = &m;
    tcase = _tcase;
#ifdef MFEM_THREAD_SAFE
    MFEM_VERIFY(0, "MPI_THREAD_SAFE should be off");
#endif
  }

  using size_type = mfem_mgis::size_type;
  using real = mfem_mgis::real;

  static constexpr const auto icste = real{0.70710678118654752440};

  inline void updateGradients(real* g,
                              const mfem::Vector& u,
                              const mfem::DenseMatrix& dN,
                              const size_type ni) noexcept {
    const auto nnodes = dN.NumRows();
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_2_2 = dN(ni, 2);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto Bi_4_0 = dN(ni, 2) * icste;
    const auto Bi_4_2 = dN(ni, 0) * icste;
    const auto Bi_5_1 = dN(ni, 2) * icste;
    const auto Bi_5_2 = dN(ni, 1) * icste;
    const auto u_0 = u[ni];
    const auto u_1 = u[ni + nnodes];
    const auto u_2 = u[ni + 2 * nnodes];
    g[0] += u_0 * Bi_0_0;
    g[1] += u_1 * Bi_1_1;
    g[2] += Bi_2_2 * u_2;
    g[3] += u_0 * Bi_3_0 + u_1 * Bi_3_1;
    g[4] += Bi_4_2 * u_2 + Bi_4_0 * u_0;
    g[5] += u_2 * Bi_5_2 + u_1 * Bi_5_1;
  }  // end of updateGradients

  inline void updateInnerForces(mfem::Vector& Fe,
                                const real* s,
                                const mfem::DenseMatrix& dN,
                                const real w,
                                const size_type ni) const noexcept {
    const auto nnodes = dN.NumRows();
    const auto Bi_0_0 = dN(ni, 0);
    const auto Bi_1_1 = dN(ni, 1);
    const auto Bi_2_2 = dN(ni, 2);
    const auto Bi_3_0 = dN(ni, 1) * icste;
    const auto Bi_3_1 = dN(ni, 0) * icste;
    const auto Bi_4_0 = dN(ni, 2) * icste;
    const auto Bi_4_2 = dN(ni, 0) * icste;
    const auto Bi_5_1 = dN(ni, 2) * icste;
    const auto Bi_5_2 = dN(ni, 1) * icste;
    const auto ni_0 = ni;
    const auto ni_1 = ni + nnodes;
    const auto ni_2 = ni + 2 * nnodes;
    Fe[ni_0] += w * (s[3] * Bi_3_0 + s[4] * Bi_4_0 + s[0] * Bi_0_0);
    Fe[ni_1] += w * (s[1] * Bi_1_1 + s[3] * Bi_3_1 + Bi_5_1 * s[5]);
    Fe[ni_2] += w * (Bi_4_2 * s[4] + Bi_2_2 * s[2] + s[5] * Bi_5_2);
  }  // end of updateInnerForces

  /*!
   * \brief Compute part of the inner forces using a single element.
   *
   * this function is called automatically by NonlinearForm::Mult
   * to perform the assembly of the RHS.
   */
  void AssembleElementVector(const mfem::FiniteElement& el,
                             mfem::ElementTransformation& Tr,
                             const mfem::Vector& U,
                             mfem::Vector& elvect) override {
    int dof = el.GetDof();
    int spaceDim = Tr.GetSpaceDim();

    dshape.SetSize(dof, spaceDim);
    Qvec.SetSize(3);
    Svec.SetSize(3);
    elvect.SetSize(dof * spaceDim);
    elvect = 0.0;

    const mfem::IntegrationRule* ir = mfem::NonlinearFormIntegrator::IntRule;
    if (ir == NULL) {
      int intorder = 2 * el.GetOrder();
      ir = &mfem::IntRules.Get(el.GetGeomType(), intorder);
    }
    const auto nnodes = el.GetDof();
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      Tr.SetIntPoint(&ip);
      el.CalcPhysDShape(Tr, dshape);
      // get the weights associated to point ip
      const auto w = ip.weight * Tr.Weight();
      // elastic properties
      double M = mu->Eval(Tr, ip);
      double L = lambda->Eval(Tr, ip);
      // macroscopic strains for each test cases
      real em[6][6] = {{1, 0, 0, 0, 0, 0},      //
                       {0, 1, 0, 0, 0, 0},      //
                       {0, 0, 1, 0, 0, 0},      //
                       {0, 0, 0, icste, 0, 0},  //
                       {0, 0, 0, 0, icste, 0},  //
                       {0, 0, 0, 0, 0, icste}};
      auto* e = &(em[tcase][0]);
      for (size_type ni = 0; ni != nnodes; ++ni) {
        this->updateGradients(e, U, dshape, ni);
      }
      const auto tr = e[0] + e[1] + e[2];
      const real sig[6] = {L * tr + 2 * M * e[0],  //
                           L * tr + 2 * M * e[1],  //
                           L * tr + 2 * M * e[2],  //
                           2 * M * e[3],           //
                           2 * M * e[4],           //
                           2 * M * e[5]};
      for (size_type ni = 0; ni != nnodes; ++ni) {
        this->updateInnerForces(elvect, sig, dshape, w, ni);
      }
    }
  }  // end of AssembleElementVector

  /*!
   * \brief Compute part of the stiffness matrix using a single element
   *
   * this function is called automatically by NonlinearForm::GetGradient
   * to perform the assembly of Stiffness matrix.
   */
  void AssembleElementGrad(const mfem::FiniteElement& el,
                           mfem::ElementTransformation& Trans,
                           const mfem::Vector&,
                           mfem::DenseMatrix& elmat) override {
    int dof = el.GetDof();
    int dim = el.GetDim();
    double w, L, M;

    MFEM_ASSERT(dim == Trans.GetSpaceDim(), "");

    cshape.SetSize(dof);
    gshape.SetSize(dof, dim);
    pelmat.SetSize(dof);  // size dof*dof
    divshape.SetSize(dim * dof);
    elmat.SetSize(dof * dim);

    const mfem::IntegrationRule* ir = mfem::NonlinearFormIntegrator::IntRule;
    if (ir == NULL) {
      int order = 2 * Trans.OrderGrad(&el);  // correct order?
      ir = &mfem::IntRules.Get(el.GetGeomType(), order);
    }

    elmat = 0.0;

    for (int nn = 0; nn < ir->GetNPoints(); nn++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(nn);

      Trans.SetIntPoint(&ip);
      // Each row of the result dshape contains
      // the derivatives of one shape function at the point ip.
      el.CalcPhysDShape(Trans, gshape);
      el.CalcPhysShape(Trans, cshape);
      // Get the transformation Trans for point ip
      // Get the weights associated to point ip
      w = ip.weight * Trans.Weight();
      // Multiply the derivatives by the inverse jacobian matrix
      // to get the derivatives along x, y and z
      // gshape contains these derivatives gshape(dof,dim)
      M = mu->Eval(Trans, ip);
      L = lambda->Eval(Trans, ip);

      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++) {
          for (int k = 0; k < dof; k++)
            for (int l = 0; l < dof; l++) {
              elmat(dof * i + k, dof * j + l) +=
                  (L * w) * gshape(k, i) * gshape(l, j) +
                  (M * w) * gshape(k, j) * gshape(l, i);
            }
          if (j == i) {
            for (int k = 0; k < dof; k++)
              for (int l = 0; l < dof; l++) {
                double ftmp = 0;
                for (int d = 0; d < dim; d++)
                  ftmp += gshape(k, d) * gshape(l, d);
                elmat(dof * i + k, dof * j + l) += (M * w) * ftmp;
              }
          }
        }
    }
  }  // end of AssembleElementGrad

 protected:
  // Coefficients related to materials
  mfem::Coefficient *lambda, *mu;
  int tcase;

 private:
  // Several temporary buffers
  mfem::Vector shape, Qvec, Svec;
  mfem::DenseMatrix dshape;
  mfem::Vector cshape;
  mfem::DenseMatrix gshape, pelmat;
  mfem::Vector divshape;
};

void executeMFEMTest(const TestParameters& p) {
  PROFILER_START(2_read_mesh);
  constexpr const auto dim = mfem_mgis::size_type{3};
  // creating the finite element workspace
  mfem::Mesh* smesh = new mfem::Mesh(p.mesh_file, 1, 1);
  if (dim != smesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  mfem::ParMesh *pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *smesh);
  delete smesh;
  MPI_Barrier(MPI_COMM_WORLD);
  std::shared_ptr<mfem::ParMesh> mesh(pmesh);

  if (dim != mesh->Dimension()) {
    std::cerr << "Invalid mesh dimension \n";
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  PROFILER_END(); PROFILER_START(3_refine_mesh);

  for (int i = 0 ; i < p.refine ; i++)
    mesh->UniformRefinement();

  PROFILER_END(); PROFILER_START(4_initialize_fem); 
  // building the non linear problem
  mfem_mgis::NonLinearEvolutionProblemBase<true> problem(
      std::make_shared<mfem_mgis::FiniteElementDiscretization>(
          mesh, std::make_shared<mfem::H1_FECollection>(p.order, dim), 3));
  mfem::Vector lambda(mesh->attributes.Max());
  lambda = 100.0;
  if (mesh->attributes.Max() > 1) lambda(1) = lambda(0) * 2;
  mfem::PWConstCoefficient lambda_func(lambda);
  // Question: pourquoi piece wise constant ?
  // lié aux attributs : constant par zone matériau
  mfem::Vector mu(mesh->attributes.Max());
  mu = 75.0;
  if (mesh->attributes.Max() > 1) mu(1) = mu(0) * 2;
  mfem::PWConstCoefficient mu_func(mu);
  problem.AddDomainIntegrator(
      new ElasticityNonLinearIntegrator(lambda_func, mu_func, p.tcase));
  //
  setBoundaryConditions(problem);
  //
  auto lsolver = getLinearSolver(p.linearsolver);
  setSolverParameters(problem, *(lsolver.get()));
  // solving the problem
  PROFILER_END(); PROFILER_START(5_solve);
  problem.solve(1);
  PROFILER_END(); PROFILER_START(6_postprocess);
  //
  if (!checkSolution(problem, p.tcase)) {
    MPI_Finalize();
    std::exit(EXIT_FAILURE);
  }
  //
  exportResults(problem, p.tcase);
  PROFILER_END(); 
}

int main(int argc, char* argv[]) {

  PROFILER_ENABLE;
   // 1. Initialize MPI.
  int num_procs, myid;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  PROFILER_START(0_total);
  PROFILER_START(1_initialize);
  const auto p = parseCommandLineOptions(argc, argv);
  PROFILER_END();
  if (p.algo == 0) {
    executeMFEMMGISTest(p);
  }
  else {
    executeMFEMTest(p);
  }
  PROFILER_END();
  if (myid == 0)
    LogProfiler();
  PROFILER_DISABLE;
  MPI_Finalize();
  return EXIT_SUCCESS;
}
