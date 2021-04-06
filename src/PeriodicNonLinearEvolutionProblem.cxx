/*!
 * \file   src/PeriodicNonLinearEvolutionProblem.cxx
 * \brief
 * \author Thomas Helfer
 * \date   10/03/2021
 */

#include "MGIS/Raise.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"

namespace mfem_mgis {

#ifdef MFEM_USE_MPI

  void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<true>& p) {
    const auto* const mesh = p.getFiniteElementSpace().GetMesh();
    const auto dim = mesh->Dimension();
    mfem::Array<int> ess_tdof_list;
    ess_tdof_list.SetSize(0);
    mfem::GridFunction nodes(&p.getFiniteElementSpace());
    int found = 0;
    bool reorder_space = true;
    mesh->GetNodes(nodes);
    const auto size = nodes.Size() / dim;
    // Traversal of all dofs to detect which one is (0,0,0)
    for (int i = 0; i < size; ++i) {
      double coord[dim];  // coordinates of a node
      double dist = 0.;
      for (int j = 0; j < dim; ++j) {
        if (reorder_space)
          coord[j] = (nodes)[j * size + i];
        else
          coord[j] = (nodes)[i * dim + j];
        // because of periodic BC, 0. is also 1.
        if (abs(coord[j] - 1.) < 1e-7) coord[j] = 0.;
        dist += coord[j] * coord[j];
      }
      // If distance is close to zero, we have our reference point
      if (dist < 1.e-16) {
        for (int j = 0; j < dim; ++j) {
          int id_unk;
          if (reorder_space) {
            // id_unk = (j * size + i);
            id_unk = p.getFiniteElementSpace().GetLocalTDofNumber(j * size + i);
          } else {
            // id_unk = (i * dim + j);
            id_unk = p.getFiniteElementSpace().GetLocalTDofNumber(i * dim + j);
          }
          if (id_unk >= 0) {
            found = 1;
            ess_tdof_list.Append(id_unk);
          }
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MFEM_VERIFY(found, "Reference point at (0,0) was not found");
    p.SetEssentialTrueDofs(ess_tdof_list);
  }  // end of setPeriodicBoundaryConditions

#endif

  void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<false>& p) {
    const auto mesh = p.getFiniteElementSpace().GetMesh();
    const auto dim = mesh->Dimension();
    mfem::Array<int> ess_tdof_list;
    const auto nnodes = p.getFiniteElementSpace().GetTrueVSize() / dim;
    ess_tdof_list.SetSize(dim);
    for (int k = 0; k < dim; k++) {
      int tgdof = 0 + k * nnodes;
      ess_tdof_list[k] = tgdof;
    }
    p.SetEssentialTrueDofs(ess_tdof_list);
  }  // end of setPeriodicBoundaryConditions

  PeriodicNonLinearEvolutionProblem::PeriodicNonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed)
      : NonLinearEvolutionProblem(fed,
                                  mgis::behaviour::Hypothesis::TRIDIMENSIONAL) {
    if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      setPeriodicBoundaryConditions(this->getImplementation<true>());
#else
      mgis::raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "unsupported parallel computations");
#endif
    } else {
      setPeriodicBoundaryConditions(this->getImplementation<false>());
    }
  }  // end of PeriodicNonLinearEvolutionProblem

  void PeriodicNonLinearEvolutionProblem::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition>) {
    mgis::raise(
        "PeriodicNonLinearEvolutionProblem::addBoundaryCondition: "
        "invalid call");
  }  // end of addBoundaryCondition

  void PeriodicNonLinearEvolutionProblem::setMacroscopicGradientsEvolution(
      const std::function<std::vector<real>(const real)>& ev) {
    this->macroscopic_gradients_evolution = ev;
  }  // end of setMacroscopicGradientsEvolution

  std::vector<real> PeriodicNonLinearEvolutionProblem::getMacroscopicGradients(
      const real t, const real dt) const {
    if (!this->macroscopic_gradients_evolution) {
      mgis::raise(
          "PeriodicNonLinearEvolutionProblem::getMacroscopicGradients: "
          "the evolution of the macroscopic gradients has not been set");
    }
    return this->macroscopic_gradients_evolution(t + dt);
  }  // end of getMacroscopicGradients

  void PeriodicNonLinearEvolutionProblem::setup(const real t, const real dt) {
    auto& impl = dynamic_cast<NonLinearEvolutionProblemImplementationBase&>(
        *(this->pimpl));
    impl.setMacroscopicGradients(this->getMacroscopicGradients(t, dt));
  }  // end of setup

  PeriodicNonLinearEvolutionProblem::~PeriodicNonLinearEvolutionProblem() =
      default;

}  // end of namespace mfem_mgis