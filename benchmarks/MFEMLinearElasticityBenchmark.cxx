/// Modified version of the MFEM example 2  ///

#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[]) {
  // 1. Initialize MPI and HYPRE.
  Mpi::Init(argc, argv);
  int myid = Mpi::WorldRank();
  Hypre::Init();

  // 2. Parse command-line options.
  const char *mesh_file = "../beam-tet.mesh";
  int order = 1;
  int refinement = 3;
  int verbosity = 1;
  bool visualization = 0;
  bool reorder_space = false;

  OptionsParser args(argc, argv);
  args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&reorder_space, "-nodes", "--by-nodes", "-vdim", "--by-vdim",
                 "Use byNODES ordering of vector space instead of byVDIM");
  args.AddOption(&refinement, "-r", "--refinenemt",
                 "Define the uniform refinement level");
  args.AddOption(&verbosity, "-v", "--verbosity", "Linear solver verbosity");
  args.Parse();
  if (!args.Good()) {
    if (myid == 0) {
      args.PrintUsage(cout);
    }
    return 1;
  }
  if (myid == 0) {
    args.PrintOptions(cout);
  }

  double start, end;
  double assemble_t1, assemble_t2;
  double solve_t1, solve_t2;
  start = MPI_Wtime();

  /** Read Mesh */
  Mesh *mesh = new Mesh(mesh_file, 1, 1);
  int dim = mesh->Dimension();
  if (mesh->attributes.Max() < 2 || mesh->bdr_attributes.Max() < 2) {
    if (myid == 0)
      cerr << "\nInput mesh should have at least two materials and "
           << "two boundary attributes! (See schematic in ex2.cpp)\n"
           << endl;
    return 3;
  }

  ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
  delete mesh;
  for (int l = 0; l < refinement; l++) {
    pmesh->UniformRefinement();
  }

  FiniteElementCollection *fec;
  ParFiniteElementSpace *fespace;
  fec = new H1_FECollection(order, dim);
  if (reorder_space) {
    fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byNODES);
  } else {
    fespace = new ParFiniteElementSpace(pmesh, fec, dim, Ordering::byVDIM);
  }
  HYPRE_BigInt size = fespace->GlobalTrueVSize();
  if (myid == 0) {
    cout << "Number of finite element unknowns: " << size << endl
         << "Assembling: " << flush;
  }

  /** Boundary conditions */
  Array<int> ess_tdof_list, ess_bdr(pmesh->bdr_attributes.Max());
  ess_bdr = 0;
  ess_bdr[0] = 1;
  ess_bdr[1] = 1;
  fespace->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

  ParLinearForm *b = new ParLinearForm(fespace);
  if (myid == 0) {
    cout << "r.h.s. ... " << flush;
  }

  ParGridFunction x(fespace);
  x = 0.0;

  {
    ess_bdr[0] = 0;
    ess_bdr[1] = 1;
    VectorArrayCoefficient bdc_value(3);
    for (int i = 0; i < 2; i++) {
      bdc_value.Set(i, new ConstantCoefficient(0.0));
    }
    {
      Vector value(pmesh->bdr_attributes.Max());
      value = 0.0;
      value(1) = -1;  // -1.0e-2;
      bdc_value.Set(2, new PWConstCoefficient(value));
    }

    x.ProjectBdrCoefficient(bdc_value, ess_bdr);
  }

  // 12. Set up the parallel bilinear form a(.,.) on the finite element space
  //     corresponding to the linear elasticity integrator with piece-wise
  //     constants coefficient lambda and mu.
  Vector lambda(pmesh->attributes.Max());
  lambda = 1.0;
  lambda(0) = lambda(1) * 50;
  PWConstCoefficient lambda_func(lambda);
  Vector mu(pmesh->attributes.Max());
  mu = 1.0;
  mu(0) = mu(1) * 50;
  PWConstCoefficient mu_func(mu);

  ParBilinearForm *a = new ParBilinearForm(fespace);
  a->AddDomainIntegrator(new ElasticityIntegrator(lambda_func, mu_func));

  /** Assemble */
  if (myid == 0) {
    cout << "matrix ... " << flush;
  }
  assemble_t1 = MPI_Wtime();
  b->Assemble();
  a->Assemble();
  assemble_t2 = MPI_Wtime();

  HypreParMatrix A;
  Vector B, X;
  a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
  if (myid == 0) {
    cout << "done." << endl;
    cout << "Size of linear system: " << A.GetGlobalNumRows() << endl;
  }

  /** Define and Use Linear Solver */
  HypreDiagScale *ds = new HypreDiagScale(A);
  // HypreGMRES *pcg = new HypreGMRES(A);
  HyprePCG *pcg = new HyprePCG(A);
  pcg->SetTol(1e-14);
  pcg->SetMaxIter(10000);
  pcg->SetPrintLevel(verbosity);
  pcg->SetPreconditioner(*ds);

  solve_t1 = MPI_Wtime();
  pcg->Mult(B, X);
  solve_t2 = MPI_Wtime();

  a->RecoverFEMSolution(X, *b, x);

  /** Paraview output */
  if (visualization) {
    ParaViewDataCollection paraview_dc("Example2P", pmesh);
    paraview_dc.SetPrefixPath("ParaView");
    paraview_dc.RegisterField("Displacement", &x);
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.SetDataFormat(VTKFormat::BINARY);
    paraview_dc.SetHighOrderOutput(true);
    paraview_dc.SetCycle(0);
    paraview_dc.SetTime(0.0);
    paraview_dc.Save();
  }

  end = MPI_Wtime();

  if (myid == 0) {
    printf("Assembly:    %1.6f s \n", assemble_t2 - assemble_t1);
    fflush(stdout);
    printf("Solve(Mult): %1.6f s \n", solve_t2 - solve_t1);
    fflush(stdout);
    printf("Duration:    %1.6f s \n", end - start);
    fflush(stdout);
  }

  // 19. Free the used memory.
  delete pcg;
  delete ds;
  delete a;
  delete b;
  if (fec) {
    delete fespace;
    delete fec;
  }
  delete pmesh;

  return 0;
}
