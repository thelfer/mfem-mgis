/*!
 * \file   src/ImposedDirichletBoundaryConditionAtClosestNode.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   31/03/2022
 */

#ifdef MFEM_MGIS_DEBUG
#include <iostream>
#endif /* MFEM_MGIS_DEBUG */

#include "mfem/mesh/mesh.hpp"
#include "mfem/fem/fespace.hpp"
#include "mfem/fem/gridfunc.hpp"
#ifdef MFEM_USE_MPI
#include <mfem/mesh/pmesh.hpp>
#include <mfem/fem/pfespace.hpp>
#include <mfem/fem/pgridfunc.hpp>
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"

namespace mfem_mgis {

  template <std::size_t space_dimension, bool parallel>
  static std::optional<size_type> getDegreeOfFreedomForClosestNode(
      FiniteElementSpace<parallel>& fes,
      const Mesh<parallel>& mesh,
      std::array<real, space_dimension> pt,
      const size_type c) {
    const auto dim = mesh.Dimension();
    if (dim != space_dimension) {
      raise(
          "getDegreeOfFreedomForClosestNode: "
          "unmatched space dimension");
    }
    mfem::GridFunction nodes(&fes);
    mesh.GetNodes(nodes);
    const auto bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
    const auto size = nodes.Size() / dim;
    // Traversal of all dofs to detect which one is the corner
    std::optional<size_type> dof_id = size_type{};
    auto min_distance = std::numeric_limits<real>::max();
#ifdef MFEM_MGIS_DEBUG
    std::array<real, space_dimension> pt_min;
    auto print_node_position = [&pt_min] {
      std::cout << "node position:";
      for (size_type j = 0; j < space_dimension; ++j) {
        std::cout << " " << pt_min[j];
      }
      std::cout << '\n';
    };
#endif /* MFEM_MGIS_DEBUG */
    for (size_type i = 0; i < size; ++i) {
      auto d2 = real{};
      for (size_type j = 0; j < space_dimension; ++j) {
        if (bynodes) {
          const auto pc = nodes[j * size + i];
          d2 += (pc - pt[j]) * (pc - pt[j]);
        } else {
          const auto pc = nodes[i * space_dimension + j];
          d2 += (pc - pt[j]) * (pc - pt[j]);
        }
      }
      const auto d = std::sqrt(d2);
      if (d < min_distance) {
        if (bynodes) {
          if constexpr (parallel) {
            dof_id = fes.GetLocalTDofNumber(c * size + i);
          } else {
            dof_id = c * size + i;
          }
#ifdef MFEM_MGIS_DEBUG
          for (size_type j = 0; j < space_dimension; ++j) {
            pt_min[j] = nodes[j * size + i];
          }
#endif /* MFEM_MGIS_DEBUG */
        } else {
          if constexpr (parallel) {
            dof_id = fes.GetLocalTDofNumber(i * space_dimension + c);
          } else {
            dof_id = i * space_dimension + c;
          }
#ifdef MFEM_MGIS_DEBUG
          for (size_type j = 0; j < space_dimension; ++j) {
            pt_min[j] = nodes[i * space_dimension + j];
          }
#endif /* MFEM_MGIS_DEBUG */
        }
        min_distance = d;
      }
    }
#ifdef MFEM_USE_MPI
    if constexpr (parallel) {
      int nbranks, myrank;
      MPI_Comm_size(MPI_COMM_WORLD, &nbranks);
      MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
      std::vector<double> recv_buf(nbranks,-1);
      double d_min = min_distance;
      // gathering all minimum values overs all processes
      MPI_Allgather(&d_min, 1, MPI_DOUBLE, recv_buf.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
#ifdef MFEM_MGIS_DEBUG
      if (myrank == 0) {
        std::cout << "minimal distance per proc:";
        for (const auto& d : recv_buf) {
          std::cout << " " << d;
        }
        std::cout << '\n';
      }
#endif /* MFEM_MGIS_DEBUG */
      // locate the minimum among the mimum values
      auto result = std::min_element(recv_buf.begin(),recv_buf.end());
      // locate on which process we have the minimum
      const auto target_pid = result - recv_buf.begin();
      // if the min value belongs to my process, I register the associated unknowns
      if (target_pid != myrank) {
        dof_id.reset();
      } else {
#ifdef MFEM_MGIS_DEBUG
        std::cout << "target_pid: " << target_pid << '\n';
print_node_position();
#endif /* MFEM_MGIS_DEBUG */
      }
    } else {
#ifdef MFEM_MGIS_DEBUG
      print_node_position();
#endif /* MFEM_MGIS_DEBUG */
    }
#else /* MFEM_USE_MPI */
#ifdef MFEM_MGIS_DEBUG
    print_node_position();
#endif /* MFEM_MGIS_DEBUG */
#endif /* MFEM_USE_MPI */
    return dof_id;
  }  // end of getDegreeOfFreedomForClosestNode

  template <std::size_t space_dimension>
  static std::optional<size_type> getDegreeOfFreedomForClosestNode(
      FiniteElementDiscretization& fed,
      std::array<real, space_dimension> pt,
      const size_type c) {
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      auto& fes = fed.getFiniteElementSpace<true>();
      const auto& mesh = fed.getMesh<true>();
      return getDegreeOfFreedomForClosestNode<space_dimension, true>(fes, mesh,
                                                                     pt, c);
#else
      raise(
          "getDegreeOfFreedomForClosestNode: "
          "unsupported parallel computations");
#endif
    }
    auto& fes = fed.getFiniteElementSpace<false>();
    const auto& mesh = fed.getMesh<false>();
    return getDegreeOfFreedomForClosestNode<space_dimension, false>(fes, mesh,
                                                                    pt, c);
  }  // end of getDegreeOfFreedomForClosestNode

    ImposedDirichletBoundaryConditionAtClosestNode::
        ImposedDirichletBoundaryConditionAtClosestNode(
            const std::shared_ptr<FiniteElementDiscretization> fed,
            const std::array<real, 2u> pt, const size_type c)
        : ufct([](const real) { return 0; }),
          dof(getDegreeOfFreedomForClosestNode<2u>(*fed, pt, c)) {
    }  // end of  ImposedDirichletBoundaryConditionAtClosestNode

    ImposedDirichletBoundaryConditionAtClosestNode::
        ImposedDirichletBoundaryConditionAtClosestNode(
            const std::shared_ptr<FiniteElementDiscretization> fed,
            const std::array<real, 2u> pt, const size_type c,
            std::function<real(const real)> uvalues)
        : ufct(uvalues),
          dof(getDegreeOfFreedomForClosestNode<2u>(*fed, pt, c)) {
    }  // end of  ImposedDirichletBoundaryConditionAtClosestNode

    ImposedDirichletBoundaryConditionAtClosestNode::
        ImposedDirichletBoundaryConditionAtClosestNode(
            const std::shared_ptr<FiniteElementDiscretization> fed,
            const std::array<real, 3u> pt, const size_type c)
        : ufct([](const real) { return 0; }),
          dof(getDegreeOfFreedomForClosestNode<3u>(*fed, pt, c)) {
    }  // end of  ImposedDirichletBoundaryConditionAtClosestNode

    ImposedDirichletBoundaryConditionAtClosestNode::
        ImposedDirichletBoundaryConditionAtClosestNode(
            const std::shared_ptr<FiniteElementDiscretization> fed,
            const std::array<real, 3u> pt, const size_type c,
            std::function<real(const real)> uvalues)
        : ufct(uvalues),
          dof(getDegreeOfFreedomForClosestNode<3u>(*fed, pt, c)) {
    }  // end of  ImposedDirichletBoundaryConditionAtClosestNode

    std::vector<size_type>
    ImposedDirichletBoundaryConditionAtClosestNode::getHandledDegreesOfFreedom()
        const {
      if (!this->dof.has_value()) {
        return {};
      }
      return {*(this->dof)};
  }  // end of getHandledDegreesOfFreedom

  void ImposedDirichletBoundaryConditionAtClosestNode::updateImposedValues(
      mfem::Vector& u, const real t) const {
    if (!this->dof.has_value()) {
      return;
    }
    const auto uv = this->ufct(t);
    u[*(this->dof)] = uv;
  }  // end of updateImposedValues

  void
  ImposedDirichletBoundaryConditionAtClosestNode::setImposedValuesIncrements(
      mfem::Vector& du, const real ti, const real te) const {
    if (!this->dof.has_value()) {
      return;
    }
    const auto duv = this->ufct(te) - this->ufct(ti);
    du[*(this->dof)] = duv;
  }  // end of setImposedValuesIncrements

  ImposedDirichletBoundaryConditionAtClosestNode::
      ~ImposedDirichletBoundaryConditionAtClosestNode() = default;

  }  // end of namespace mfem_mgis