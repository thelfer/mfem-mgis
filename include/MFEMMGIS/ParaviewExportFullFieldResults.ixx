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
    std::cout << "Start ParaviewExportFullFieldResults" << "\n";
    checkParameters(params, {"OutputFileName", "Materials", "Results"});
    // if Materials exists, use it, otherwise, take all materials
    this->materials_identifiers = getMaterialsIdentifiers(p, params);

    const auto& u1 = p.getUnknownsAtEndOfTheTimeStep();
    this->uTot = u1;

    this->result.MakeTRef(&p.getFiniteElementSpace(), this->uTot, 0);
    if (contains(params, "OutputFieldName")) {
      this->exporter.RegisterField(get<std::string>(params, "OutputFieldName"),
                                   &(this->result));
    } else {
      this->exporter.RegisterField("u", &(this->result));
    }
    std::cout << "End ParaviewExportFullFieldResults" << "\n";
  }  // end of ParaviewExportFullFieldResults

  template <bool parallel>
  void ParaviewExportFullFieldResults<parallel>::execute(
      NonLinearEvolutionProblemImplementation<parallel>& p,
      const real t,
      const real dt) {
    std::cout << "Start ParaviewExportFullFieldResults::execute"
              << "\n";
    this->exporter.SetCycle(this->cycle);
    this->exporter.SetTime(t + dt);

    const auto& mids = this->materials_identifiers;
    // tester si les hypothèses de déformation sont homogènes entre les cristaux
    // ?
    const auto& m = p.getMaterial(mids[0]);
    const auto& b = m.b;
    auto mGrad = m.getMacroscopicGradients();

    if (b.btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
      std::cout << "SMALL STRAIN BEHAVIOUR" << '\n';
      for (auto cpt : mGrad) {
        std::cout << "E = " << cpt << '\n';
      }
    } else if (b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
      for (int i = 0; i < 3; i++) {
        mGrad[i] -= 1;
      }
      std::cout << "FINITE STRAIN BEHAVIOUR" << '\n';
      for (auto cpt : mGrad) {
        std::cout << "F = " << cpt << '\n';
      }
    } else {
      raise("invalid behaviour type");
    }

    this->uTot = p.getUnknownsAtEndOfTheTimeStep();

    const FiniteElementSpace<parallel>& fes = p.getFiniteElementSpace();
    const auto* const mesh = p.getFiniteElementSpace().GetMesh();
    const auto dim = mesh->Dimension();
    mfem::GridFunction nodes(&p.getFiniteElementSpace());
    mesh->GetNodes(nodes);
    bool bynodes = fes.GetOrdering() == mfem::Ordering::byNODES;
    const auto size = nodes.Size() / dim;

    for (int i = 0; i < size; ++i) {
      double coord[dim];  // coordinates of a node
      for (int j = 0; j < dim; ++j) {
        if (bynodes)
          coord[j] = (nodes)[j * size + i];
        else
          coord[j] = (nodes)[i * dim + j];
      }

      for (int j = 0; j < dim; ++j) {
        int id_unk;
        if (bynodes) {
          if constexpr (parallel)
            id_unk = fes.GetLocalTDofNumber(j * size + i);
          else
            id_unk = (j * size + i);
        } else {
          if constexpr (parallel)
            id_unk = fes.GetLocalTDofNumber(i * dim + j);
          else
            id_unk = (i * dim + j);
        }
        if (id_unk >= 0) {
          if (b.btype == Behaviour::STANDARDSTRAINBASEDBEHAVIOUR) {
            if (j == 0) {
              this->uTot[id_unk] += mGrad[0] * coord[0] +
                                    1 / sqrt(2) * mGrad[3] * coord[1] +
                                    1 / sqrt(2) * mGrad[4] * coord[2];
            } else if (j == 1) {
              this->uTot[id_unk] += 1 / sqrt(2) * mGrad[3] * coord[0] +
                                    mGrad[1] * coord[1] +
                                    1 / sqrt(2) * mGrad[5] * coord[2];
            } else if (j == 2) {
              this->uTot[id_unk] += 1 / sqrt(2) * mGrad[4] * coord[0] +
                                    1 / sqrt(2) * mGrad[5] * coord[1] +
                                    mGrad[2] * coord[2];
            }
          } else if (b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) {
            if (j == 0) {
              this->uTot[id_unk] += mGrad[0] * coord[0] + mGrad[3] * coord[1] +
                                    mGrad[5] * coord[2];
            } else if (j == 1) {
              this->uTot[id_unk] += mGrad[4] * coord[0] + mGrad[1] * coord[1] +
                                    mGrad[7] * coord[2];
            } else if (j == 2) {
              this->uTot[id_unk] += mGrad[6] * coord[0] + mGrad[8] * coord[1] +
                                    mGrad[2] * coord[2];
            }
          }
        }
      }
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
  ParaviewExportFullFieldResults<parallel>::~ParaviewExportFullFieldResults() = default;

  }  // end of namespace mfem_mgis

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTFULLFIELDRESULTS_IXX */
