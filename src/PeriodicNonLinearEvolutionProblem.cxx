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

  real getNodesDistance(const mfem::GridFunction & nodes,
			const bool reorder_space,
			const size_t dim,
			const int index,
			const int size, 
			const mgis::span<const real>& corner1,
			const mgis::span<const real>& corner2) {
    real coord[dim];  // coordinates of a node
    real dist = 0.;
    for (int j = 0; j < dim; ++j) {
      if (reorder_space)
	coord[j] = (nodes)[j * size + index];
      else
	coord[j] = (nodes)[index * dim + j];
      real dist1 = (coord[j]-corner1[j]) * (coord[j]-corner1[j]);
      real dist2 = (coord[j]-corner2[j]) * (coord[j]-corner2[j]);
      dist += std::min(dist1,dist2);
    }
    return (dist);
  }
  
#ifdef MFEM_USE_MPI

  void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<true>& p,
      const mgis::span<const real>& corner1,
      const mgis::span<const real>& corner2) {
    const FiniteElementSpace<true> & fes = p.getFiniteElementSpace();
    const auto* const mesh = p.getFiniteElementSpace().GetMesh();
    const auto dim = mesh->Dimension();
    mfem::Array<int> ess_tdof_list;
    ess_tdof_list.SetSize(0);
    mfem::GridFunction nodes(&p.getFiniteElementSpace());
    int found = 0;
    bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
    mesh->GetNodes(nodes);
    const auto size = nodes.Size() / dim;
    
    // Traversal of all dofs to detect which one is the corner
    for (int i = 0; i < size; ++i) {
      real dist = getNodesDistance(nodes, bynodes, 
				   dim, i, size, corner1, corner2);
      // If distance is close to zero, we have our reference point
      if (dist < 1.e-12) {
        for (int j = 0; j < dim; ++j) {
          int id_unk;
          if (bynodes) {
            // id_unk = (j * size + i);
            id_unk = p.getFiniteElementSpace().GetLocalTDofNumber(j * size + i);
          } else {
            // id_unk = (i * dim + j);
            id_unk = p.getFiniteElementSpace().GetLocalTDofNumber(i * dim + j);
          }
          if (id_unk >= 0) {
            found = 1;
            ess_tdof_list.Append(id_unk);
	    std::cout << "unknown " << id_unk << std::endl;
	  }
        }
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    
    MFEM_VERIFY(found, "Corner point was not found");
    p.SetEssentialTrueDofs(ess_tdof_list);
  }  // end of setPeriodicBoundaryConditions

  void setPeriodicBoundaryConditions(
      NonLinearEvolutionProblemImplementation<true>& p,
      const mfem_mgis::BoundaryConditionType bct) {
    const FiniteElementSpace<true> & fes = p.getFiniteElementSpace();
    const auto* const mesh = p.getFiniteElementSpace().GetMesh();
    const auto dim = mesh->Dimension();
    mfem::Array<int> ess_tdof_list;
    ess_tdof_list.SetSize(0);
    mfem::GridFunction nodes(&p.getFiniteElementSpace());
    int found = 0;
    bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
    mesh->GetNodes(nodes);
    const auto size = nodes.Size() / dim;

    // Initialize reference values to largest possible numbers
    real refcoord[dim];
    int id_unk[dim];
    for (int j = 0; j < dim; ++j) {
      refcoord[j] = std::numeric_limits<double>::max();
      id_unk[j] = -1;
    }
    // Traversal of all dofs to detect which one is minimal in X, Y or Z direction
    // depending on the `bct` variable.
    for (int i = 0; i < size; ++i) {
      real curcoord[dim];
      for (int j = 0; j < dim; ++j) {
	if (bynodes)
	  curcoord[j] = (nodes)[j * size + i];
	else
	  curcoord[j] = (nodes)[i * dim + j];
      }
      if (((bct == FIX_XMIN) && (curcoord[0] < refcoord[0])) ||
	  ((bct == FIX_YMIN) && (curcoord[1] < refcoord[1])) ||
	  ((bct == FIX_ZMIN) && (curcoord[2] < refcoord[2]))) {
	for (int j = 0; j < dim; ++j) {
	  refcoord[j] = curcoord[j];
	  if (bynodes) {
	    // id_unk = (j * size + i);
	    id_unk[j] = p.getFiniteElementSpace().GetLocalTDofNumber(j * size + i);
	  } else {
	    // id_unk = (i * dim + j);
	    id_unk[j] = p.getFiniteElementSpace().GetLocalTDofNumber(i * dim + j);
	  }
	}
      }
    }
    // MPI communications to identify where is the minimum among all processes
    {
      int nbranks, myrank;
      MPI_Comm_size(MPI_COMM_WORLD, &nbranks);
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      std::vector<double> recv_buf(nbranks,-1);
      double mymin = -1;
      if (bct == FIX_XMIN) mymin = refcoord[0];
      if (bct == FIX_YMIN) mymin = refcoord[1];
      if (bct == FIX_ZMIN) mymin = refcoord[2];
      // gathering all minimum values overs all processes
      MPI_Allgather(&mymin, 1, MPI_DOUBLE, recv_buf.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
      // locate the minimum among the mimum values
      auto result = std::min_element(recv_buf.begin(),recv_buf.end());
      // locate on which process we have the minimum
      int target_pid = result-recv_buf.begin();
      // if the min value belongs to my process, I register the associated unknowns
      if (target_pid == myrank) {
	for (int j = 0; j < dim; ++j) {
	  ess_tdof_list.Append(id_unk[j]);
	}
	found = 1;
      }
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    MFEM_VERIFY(found == 1, "Not able to define proper periodic boundary conditions");
    p.SetEssentialTrueDofs(ess_tdof_list);
  }  // end of setPeriodicBoundaryConditions
#endif /* MFEM_USE_MPI */

  void setPeriodicBoundaryConditions(
    NonLinearEvolutionProblemImplementation<false>& p,
    const mgis::span<const real>& corner1,
    const mgis::span<const real>& corner2) {
    const FiniteElementSpace<false> & fes = p.getFiniteElementSpace();
    const auto* const mesh = p.getFiniteElementSpace().GetMesh();
    const auto dim = mesh->Dimension();
    mfem::Array<int> ess_tdof_list;
    ess_tdof_list.SetSize(0);
    mfem::GridFunction nodes(&p.getFiniteElementSpace());
    int found = 0;
    bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
    mesh->GetNodes(nodes);
    const auto size = nodes.Size() / dim;
    // Traversal of all dofs to detect which one is the corner
    for (int i = 0; i < size; ++i) {
      real dist = getNodesDistance(nodes, bynodes,
				   dim, i, size, corner1, corner2);
      // If distance is close to zero, we have our reference point
      if (dist < 1.e-12) {
        for (int j = 0; j < dim; ++j) {
          int id_unk;
          if (bynodes) {
	    id_unk = (j * size + i);
          } else {
	    id_unk = (i * dim + j);
          }
          if (id_unk >= 0) {
            found = 1;
            ess_tdof_list.Append(id_unk);
          }
        }
      }
    }
    
    MFEM_VERIFY(found, "Corner point was not found");
    p.SetEssentialTrueDofs(ess_tdof_list);
  }  // end of setPeriodicBoundaryConditions

  void setPeriodicBoundaryConditions(
    NonLinearEvolutionProblemImplementation<false>& p,
    const mfem_mgis::BoundaryConditionType bct) {
    const FiniteElementSpace<false> & fes = p.getFiniteElementSpace();
    const auto* const mesh = p.getFiniteElementSpace().GetMesh();
    const auto dim = mesh->Dimension();
    mfem::Array<int> ess_tdof_list;
    ess_tdof_list.SetSize(0);
    mfem::GridFunction nodes(&p.getFiniteElementSpace());
    int found = 0;
    bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
    mesh->GetNodes(nodes);
    const auto size = nodes.Size() / dim;

    // Initialize reference values to largest possible numbers
    real refcoord[dim];
    int id_unk[dim];
    for (int j = 0; j < dim; ++j) {
      refcoord[j] = std::numeric_limits<double>::max();
      id_unk[j] = -1;
    }

    // Traversal of all dofs to detect which one is minimal in X, Y or Z direction
    // depending on the `bct` variable.
    for (int i = 0; i < size; ++i) {
      real curcoord[dim];
      for (int j = 0; j < dim; ++j) {
	if (bynodes)
	  curcoord[j] = (nodes)[j * size + i];
	else
	  curcoord[j] = (nodes)[i * dim + j];
      }
      if ((bct == FIX_XMIN) && (curcoord[0] < refcoord[0]) ||
	  (bct == FIX_YMIN) && (curcoord[1] < refcoord[1]) ||
	  (bct == FIX_ZMIN) && (curcoord[2] < refcoord[2])) {
	for (int j = 0; j < dim; ++j) {
	  found = 1;
	  refcoord[j] = curcoord[j];
	  if (bynodes) {
	    id_unk[j] = (j * size + i);
	  } else {
	    id_unk[j] = (i * dim + j);
	  }
	}
      }
    }
    MFEM_VERIFY(found == 1, "Not able to define proper periodic boundary conditions");
    for (int j = 0; j < dim; ++j) {
      ess_tdof_list.Append(id_unk[j]);
    }
    p.SetEssentialTrueDofs(ess_tdof_list);
  }  // end of setPeriodicBoundaryConditions
  
  PeriodicNonLinearEvolutionProblem::PeriodicNonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const mgis::span<const real>& corner1,
      const mgis::span<const real>& corner2)
      : NonLinearEvolutionProblem(fed,
                                  mgis::behaviour::Hypothesis::TRIDIMENSIONAL) {
    if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      setPeriodicBoundaryConditions(this->getImplementation<true>(), corner1, corner2);
#else
      raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "unsupported parallel computations");
#endif
    } else {
      setPeriodicBoundaryConditions(this->getImplementation<false>(), corner1, corner2);
    }
  }  // end of PeriodicNonLinearEvolutionProblem

  PeriodicNonLinearEvolutionProblem::PeriodicNonLinearEvolutionProblem(
      std::shared_ptr<FiniteElementDiscretization> fed,
      const mfem_mgis::BoundaryConditionType bct)
      : NonLinearEvolutionProblem(fed,
                                  mgis::behaviour::Hypothesis::TRIDIMENSIONAL) {
    if (fed->describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      setPeriodicBoundaryConditions(this->getImplementation<true>(), bct);
#else
      raise(
          "NonLinearEvolutionProblem::NonLinearEvolutionProblem: "
          "unsupported parallel computations");
#endif
    } else {
      setPeriodicBoundaryConditions(this->getImplementation<false>(), bct);
    }
  }  // end of PeriodicNonLinearEvolutionProblem

  void PeriodicNonLinearEvolutionProblem::addBoundaryCondition(
      std::unique_ptr<DirichletBoundaryCondition>) {
    raise(
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
      raise(
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
