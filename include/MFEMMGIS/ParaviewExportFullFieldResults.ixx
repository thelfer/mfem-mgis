/*!
 * \file   include/MFEMMGIS/ParaviewExportFullFieldResults.ixx
 * \brief
 * \author Maxence Wangermez
 * \date   10/10/2023
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTFULLFIELDRESULTS_IXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTFULLFIELDRESULTS_IXX

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"

namespace mfem_mgis {

  template <bool parallel>
  ParaviewExportFullFieldResults<parallel>::ParaviewExportFullFieldResults(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const Parameters& params)
      : exporter(get<std::string>(params, "OutputFileName"),
                 p.getFiniteElementSpace().GetMesh()),
        result(&p.getFiniteElementSpace()),
        cycle(0) {
    const auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
    
    // MW : Using assignment operator to copy one vector to other 
    mfem::Vector uTot;
    uTot = u1;

    // const FiniteElementSpace<parallel>& fes = p.getFiniteElementSpace();
    // const auto* const mesh = p.getFiniteElementSpace().GetMesh();
    // const auto dim = mesh->Dimension();
    // mfem::GridFunction nodes(&p.getFiniteElementSpace());
    // bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
    // const auto size = nodes.Size() / dim;
    // for (int i = 0; i < size; ++i) {
    //   const auto coord = getCoordinates(nodes, bynodes, dim, i, size);
    // }

    this->result.MakeTRef(&p.getFiniteElementSpace(), uTot, 0);
    if (contains(params, "OutputFieldName")) {
      this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
                                   &(this->result));
    } else {
      this->exporter.RegisterField("u", &(this->result));
    }
  }  // end of ParaviewExportFullFieldResults

  template <bool parallel>
  void ParaviewExportFullFieldResults<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const real t,
      const real dt) {
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t + dt);

    // MW : calculer le champs total
    std::cout << "execute from: " << t << " to: " << dt << "\n";
    // boucle sur les dofs
    for (auto& r : this->result) {
      // std::cout << "r =  " << r << "\n";
      this->computeFullFieldResult(r, t, dt);
    }
    // SetFromTrueVector needed here in MFEM for at least two rationales:
    //    - it applies prolongation matrix (Non-Conforming mesh, BCs, AMR ...)
    //      to set the values of some unkwown dofs deduced from known dofs
    //    - exchange data between processes in order to retrieve information
    //      needed to perform the previous prolongation step
    this->result.SetFromTrueVector();
    this->exporter.Save();
    ++(this->cycle);
  }  // end of execute

  template <bool parallel>
  void ParaviewExportFullFieldResults<parallel>::computeFullFieldResult(
      real& r,
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const real& t,
      const real& dt) {
    std::cout << "computeFullFieldResult from: " << t << " to: " << dt << "\n";

    const auto mids = p.getMaterialsIdentifiers();
    const auto& m = p.getMaterial(mids[0]);
    const auto& b = m.b;
    if (b->btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
    } else if (b->btype != Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
    } else {
      // mas la merde
    }

    // --- gradient Macro F ou E
    // auto& mGrad = p.getMacroscopicGradients(t, dt);

    // --- coordonnees
    // coord = p.getcoord

    // --- solution macro
    // if size(mGrad)
    // auto& uM = (F-1)*coord;
    // else
    // auto& uM = (E)*coord;

    // --- calculer le nouveau champs total
    // r += uM

  }  // end of computeFullFieldResult

  // mfem_mgis::real* getCoordinates(const mfem::GridFunction& nodes,
  //                       const bool reorder_space,
  //                       const size_t dim,
  //                       const int index,
  //                       const int size) {
  //   real coord[dim];  // coordinates of a node
  //   for (int j = 0; j < dim; ++j) {
  //     if (reorder_space)
  //       coord[j] = (nodes)[j * size + index];
  //     else
  //       coord[j] = (nodes)[index * dim + j];
  //   }
  //   return (coord);
  // }

  template <bool parallel>
  ParaviewExportFullFieldResults<parallel>::~ParaviewExportFullFieldResults() = default;

  }  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTFULLFIELDRESULTS_IXX */
