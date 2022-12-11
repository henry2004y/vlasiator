/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/*!\file wave.cpp
 * \brief Implementation of the class SysBoundaryCondition::Wave.
 * This is the interface for customized boundary implementation.
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>

#include "../fieldsolver/fs_common.h"
#include "../object_wrapper.h"
#include "../vlasovmover.h"
#include "wave.h"

#ifndef NDEBUG
#define DEBUG_WAVE
#endif
#ifdef DEBUG_SYSBOUNDARY
#define DEBUG_WAVE
#endif

namespace SBC {

void Wave::addParameters() {
   typedef Readparameters RP;
   RP::add("wave.type", "Prescribed wave type given velocity perturbations (\"Alfven\", \"Fast\", \"None\")", "None");
   RP::add("wave.precedence", "Precedence value of the wave boundary condition (integer), the higher the stronger.", 4);
   RP::add("wave.reapplyUponRestart",
           "If 0 (default), keep going with the state existing in the restart file. If 1, calls again "
           "applyInitialState. Can be used to change boundary condition behaviour during a run.",
           0);
   RP::add("wave.Bx0", "Wave background B field x", 0.0);
   RP::add("wave.By0", "Wave background B field y", 0.0);
   RP::add("wave.Bz0", "Wave background B field z", 0.0);
   RP::add("wave.Bx1", "Perturbed B field magnitude x", 0.0);
   RP::add("wave.By1", "Perturbed B field magnitude y", 0.0);
   RP::add("wave.Bz1", "Perturbed B field magnitude z", 0.0);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      RP::addComposing(pop + "_wave.face",
                       "List of faces on which outflow boundary conditions are to be applied ([xyz][+-]).");
      RP::add(pop + "_wave.kx", "Incoming wave vector x-component.", 1.0);
      RP::add(pop + "_wave.ky", "Incoming wave vector y-component.", 0.0);
      RP::add(pop + "_wave.kz", "Incoming wave vector z-component.", 0.0);
      RP::add(pop + "_wave.frequency", "Incoming wave frequency.", 0.0);

      RP::add(pop + "_wave.vx0", "Background bulk velocity in x", 0.0);
      RP::add(pop + "_wave.vy0", "Background bulk velocity in y", 0.0);
      RP::add(pop + "_wave.vz0", "Background bulk velocuty in z", 0.0);
      RP::add(pop + "_wave.n0", "Background number density", 1.0e6);
      RP::add(pop + "_wave.T0", "Background temperature", 1.0e6);

      RP::add(pop + "_wave.vx1", "Perturbed velocity magnitude in x", 0.0);
      RP::add(pop + "_wave.vy1", "Perturbed velocity magnitude in y", 0.0);
      RP::add(pop + "_wave.vz1", "Perturbed velocuty magnitude in z", 0.0);
      RP::add(pop + "_wave.n1", "Perturbed number density magnitude", 0.0);
      RP::add(pop + "_wave.T1", "Perturbed temperature magnitude", 0.0);
   }
}

void Wave::initBoundary(creal t, Project& project) {
   this->getParameters();

   dynamic = true;

   for (auto& it : faceList) {
      if (it == "x+") {
         facesToProcess[0] = true;
      } else if (it == "x-") {
         facesToProcess[1] = true;
      } else if (it == "y+") {
         facesToProcess[2] = true;
      } else if (it == "y-") {
         facesToProcess[3] = true;
      } else if (it == "z+") {
         facesToProcess[4] = true;
      } else if (it == "z-") {
         facesToProcess[5] = true;
      }
   }
}

void Wave::getParameters() {
   typedef Readparameters RP;
   RP::get("wave.type", wavetype);
   if (wavetype != "None" && wavetype != "Alfven" && wavetype != "Fast") {
      abort_mpi("Unknown boundary wave type " + wavetype);
   }

   RP::get("wave.precedence", precedence);
   uint reapply;
   RP::get("wave.reapplyUponRestart", reapply);

   this->applyUponRestart = (reapply == 1) ? true : false;

   RP::get("wave.Bx0", B0[0]);
   RP::get("wave.By0", B0[1]);
   RP::get("wave.Bz0", B0[2]);

   RP::get("wave.Bx1", B1[0]);
   RP::get("wave.By1", B1[1]);
   RP::get("wave.Bz1", B1[2]);

   // Per-species parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;

      std::vector<std::string> thisSpeciesFaceList;
      RP::get(pop + "_wave.face", thisSpeciesFaceList);

      for (auto& face : thisSpeciesFaceList) {
         if (face == "x+") {
            facesToProcess[0] = true;
         }
         if (face == "x-") {
            facesToProcess[1] = true;
         }
         if (face == "y+") {
            facesToProcess[2] = true;
         }
         if (face == "y-") {
            facesToProcess[3] = true;
         }
         if (face == "z+") {
            facesToProcess[4] = true;
         }
         if (face == "z-") {
            facesToProcess[5] = true;
         }
      }

      RP::get(pop + "_wave.kx", waveVector[0]);
      RP::get(pop + "_wave.ky", waveVector[1]);
      RP::get(pop + "_wave.kz", waveVector[2]);
      RP::get(pop + "_wave.frequency", freq);
      RP::get(pop + "_wave.n0", n0);
      RP::get(pop + "_wave.T0", T0);
      RP::get(pop + "_wave.vx0", V0[0]);
      RP::get(pop + "_wave.vy0", V0[1]);
      RP::get(pop + "_wave.vz0", V0[2]);
      RP::get(pop + "_wave.n1", n1);
      RP::get(pop + "_wave.T1", T1);
      RP::get(pop + "_wave.vx1", V1[0]);
      RP::get(pop + "_wave.vy1", V1[1]);
      RP::get(pop + "_wave.vz1", V1[2]);
   }
}

void Wave::assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                          FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   bool doAssign;
   bool isCellOnFace[6];

   // Assign boundary flags to local DCCRG cells
   const vector<CellID>& cells = getLocalCells();
   for (auto& id : cells) {
      if (mpiGrid[id]->sysBoundaryFlag == sysboundarytype::DO_NOT_COMPUTE) {
         continue;
      }
      creal* const cellParams = &(mpiGrid[id]->parameters[0]);
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5 * dx;
      creal y = cellParams[CellParams::YCRD] + 0.5 * dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5 * dz;

      determineFace(isCellOnFace, x, y, z, dx, dy, dz);

      doAssign = false;
      for (int iface = 0; iface < 6; iface++) {
         doAssign = doAssign || (facesToProcess[iface] && isCellOnFace[iface]);
      }
      if (doAssign) {
         mpiGrid[id]->sysBoundaryFlag = this->getIndex();
      }
   }

   // Assign boundary flags to local fsgrid cells
   const std::array<int, 3> gridDims(technicalGrid.getLocalSize());
   for (int k = 0; k < gridDims[2]; k++) {
      for (int j = 0; j < gridDims[1]; j++) {
         for (int i = 0; i < gridDims[0]; i++) {
            const auto coords = technicalGrid.getPhysicalCoords(i, j, k);

            // Shift to the center of the fsgrid cell
            auto cellCenterCoords = coords;
            cellCenterCoords[0] += 0.5 * technicalGrid.DX;
            cellCenterCoords[1] += 0.5 * technicalGrid.DY;
            cellCenterCoords[2] += 0.5 * technicalGrid.DZ;
            const auto refLvl = mpiGrid.get_refinement_level(mpiGrid.get_existing_cell(cellCenterCoords));

            if (refLvl == -1) {
               abort_mpi("Error, could not get refinement level of remote DCCRG cell!", 1);
            }

            creal dx = P::dx_ini * pow(2, -refLvl);
            creal dy = P::dy_ini * pow(2, -refLvl);
            creal dz = P::dz_ini * pow(2, -refLvl);

            doAssign = false;

            determineFace(isCellOnFace, cellCenterCoords[0], cellCenterCoords[1], cellCenterCoords[2], dx, dy, dz);
            for (int iface = 0; iface < 6; iface++) {
               doAssign = doAssign || (facesToProcess[iface] && isCellOnFace[iface]);
            }
            if (doAssign) {
               technicalGrid.get(i, j, k)->sysBoundaryFlag = this->getIndex();
            }
         }
      }
   }
}

void Wave::applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                             FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                             Project& project) {
   const std::vector<CellID>& cells = getLocalCells();
   namespace CP = CellParams;

   for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
      for (auto& id : cells) {
         SpatialCell* cell = mpiGrid[id];
         if (cell->sysBoundaryFlag != this->getIndex()) {
            continue;
         }

         for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
            std::vector<vmesh::GlobalID> blocksToInit = findBlocksToInit(popID, cell, n0, T0, V0[0], V0[1], V0[2]);
            Realf* f = cell->get_data(popID);

            for (vmesh::GlobalID i = 0; i < blocksToInit.size(); ++i) {
               const vmesh::GlobalID blockGID = blocksToInit[i];
               const vmesh::LocalID blockLID = cell->get_velocity_block_local_id(blockGID, popID);
               const Real* block_parameters = cell->get_block_parameters(blockLID, popID);
               creal vxBlock = block_parameters[BlockParams::VXCRD];
               creal vyBlock = block_parameters[BlockParams::VYCRD];
               creal vzBlock = block_parameters[BlockParams::VZCRD];
               creal dvx = block_parameters[BlockParams::DVX];
               creal dvy = block_parameters[BlockParams::DVY];
               creal dvz = block_parameters[BlockParams::DVZ];

               // Calculate volume-averaged distribution function for each cell in the block.
               for (uint kc = 0; kc < WID; ++kc) {
                  creal vzC = vzBlock + kc * dvz;
                  for (uint jc = 0; jc < WID; ++jc) {
                     creal vyC = vyBlock + jc * dvy;
                     for (uint ic = 0; ic < WID; ++ic) {
                        creal vxC = vxBlock + ic * dvx;
                        // hyzhou: most values are very close to 0!
                        f[blockLID * WID3 + cellIndex(ic, jc, kc)] = maxwellian(
                            popID, n0, T0, vxC + 0.5 * dvx - V0[0], vyC + 0.5 * dvy - V0[1], vzC + 0.5 * dvz - V0[2]);
                        ;
                     }
                  }
               } // loop over cells in velocity block
            }    // loop over velocity blocks

            // Store blocks according to the sparsity threshold
            cell->adjustSingleCellVelocityBlocks(popID);
         } // loop over particle species

         // For all species
         calculateCellMoments(cell, true, true);

         for (uint popID = 0; popID < getObjectWrapper().particleSpecies.size(); ++popID) {
            cell->get_population(popID).RHO_R = cell->get_population(popID).RHO;
            cell->get_population(popID).RHO_V = cell->get_population(popID).RHO;
            for (int i = 0; i < 3; i++) {
               cell->get_population(popID).V_R[i] = cell->get_population(popID).V[i];
               cell->get_population(popID).V_V[i] = cell->get_population(popID).V[i];
               cell->get_population(popID).P_R[i] = cell->get_population(popID).P[i];
               cell->get_population(popID).P_V[i] = cell->get_population(popID).P[i];
            }
         }

         cell->parameters[CP::RHOM_R] = cell->parameters[CP::RHOM];
         cell->parameters[CP::VX_R] = cell->parameters[CP::VX];
         cell->parameters[CP::VY_R] = cell->parameters[CP::VY];
         cell->parameters[CP::VZ_R] = cell->parameters[CP::VZ];
         cell->parameters[CP::RHOQ_R] = cell->parameters[CP::RHOQ];
         cell->parameters[CP::P_11_R] = cell->parameters[CP::P_11];
         cell->parameters[CP::P_22_R] = cell->parameters[CP::P_22];
         cell->parameters[CP::P_33_R] = cell->parameters[CP::P_33];
         cell->parameters[CP::RHOM_V] = cell->parameters[CP::RHOM];
         cell->parameters[CP::VX_V] = cell->parameters[CP::VX];
         cell->parameters[CP::VY_V] = cell->parameters[CP::VY];
         cell->parameters[CP::VZ_V] = cell->parameters[CP::VZ];
         cell->parameters[CP::RHOQ_V] = cell->parameters[CP::RHOQ];
         cell->parameters[CP::P_11_V] = cell->parameters[CP::P_11];
         cell->parameters[CP::P_22_V] = cell->parameters[CP::P_22];
         cell->parameters[CP::P_33_V] = cell->parameters[CP::P_33];
      }
   }
}

void Wave::updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                       FsGrid<array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, creal t) {}

Real Wave::maxwellian(const uint popID, creal rho, creal T, creal vx, creal vy, creal vz) {
   const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
   return rho * MASS / (2.0 * M_PI * physicalconstants::K_B * T) *
          sqrt(MASS / (2.0 * M_PI * physicalconstants::K_B * T)) *
          exp(-MASS * (vx * vx + vy * vy + vz * vz) / (2.0 * physicalconstants::K_B * T));
}

void Wave::vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
                                   const CellID& cellID, const uint popID, creal t, const bool doCalcMomentsV) {
   SpatialCell* cell = mpiGrid[cellID];

   creal dx = cell->parameters[CellParams::DX];
   creal dy = cell->parameters[CellParams::DY];
   creal dz = cell->parameters[CellParams::DZ];
   creal x = cell->parameters[CellParams::XCRD] + 0.5 * dx;
   creal y = cell->parameters[CellParams::YCRD] + 0.5 * dy;
   creal z = cell->parameters[CellParams::ZCRD] + 0.5 * dz;

   bool isCellOnFace[6];
   determineFace(isCellOnFace, x, y, z, dx, dy, dz, true);

   if (wavetype == "None") {
      for (uint iface = 0; iface < 6; iface++) {
         if (facesToProcess[iface] && isCellOnFace[iface]) {
            Real rho, T, Vx, Vy, Vz, tmp;
            tmp = std::sin(2 * PI * freq * t + waveVector[0] * x + waveVector[1] * y + waveVector[2] * z);
            rho = n0 + n1 * tmp;
            T = T0 + T1 * tmp;
            Vx = V0[0] + V1[0] * tmp;
            Vy = V0[1] + V1[1] * tmp;
            Vz = V0[2] + V1[2] * tmp;

            updateVDF(cell, popID, rho, T, Vx, Vy, Vz);

            break; // Effectively sets the precedence of faces through the order of faces.
         }
      }
   } else if (wavetype == "Alfven") {
      for (uint iface = 0; iface < 6; iface++) {
         if (facesToProcess[iface] && isCellOnFace[iface]) {
            Real rho, T, Vx, Vy, Vz, tmp;
            tmp = std::sin(2 * PI * freq * t + waveVector[0] * x + waveVector[1] * y + waveVector[2] * z);
            rho = n0;
            T = T0;
            Vx = V0[0] + V1[0] * tmp;
            Vy = V0[1] + V1[1] * tmp;
            Vz = V0[2] + V1[2] * tmp;

            updateVDF(cell, popID, rho, T, Vx, Vy, Vz);

            break; // Effectively sets the precedence of faces through the order of faces.
         }
      }
   } else if (wavetype == "Fast") {
      for (uint iface = 0; iface < 6; iface++) {
         if (facesToProcess[iface] && isCellOnFace[iface]) {
            const Real MASS = getObjectWrapper().particleSpecies[popID].mass;
            Real rho, T, Vx, Vy, Vz, tmp, kdotv, omega;
            omega = 2 * PI * freq;
            tmp = std::sin(2 * PI * freq * t + waveVector[0] * x + waveVector[1] * y + waveVector[2] * z);
            kdotv = std::inner_product(std::begin(waveVector), std::end(waveVector), std::begin(V1), 0.0);
            rho = n0 / omega * kdotv;
            T = T0 * (P::electronPTindex - 1) * kdotv / (MASS * (omega + kdotv));
            Vx = V0[0] + V1[0] * tmp;
            Vy = V0[1] + V1[1] * tmp;
            Vz = V0[2] + V1[2] * tmp;

            updateVDF(cell, popID, rho, T, Vx, Vy, Vz);

            break; // Effectively sets the precedence of faces through the order of faces.
         }
      }
   }
}

inline void Wave::updateVDF(spatial_cell::SpatialCell* cell, const uint popID, creal rho, creal T, creal Vx, creal Vy,
                            creal Vz) {
   std::vector<vmesh::GlobalID> blocksToInit = findBlocksToInit(popID, cell, rho, T, Vx, Vy, Vz);
   Realf* f = cell->get_data(popID);

   for (vmesh::GlobalID i = 0; i < blocksToInit.size(); ++i) {
      const vmesh::GlobalID blockGID = blocksToInit[i];
      const vmesh::LocalID blockLID = cell->get_velocity_block_local_id(blockGID, popID);
      const Real* block_parameters = cell->get_block_parameters(blockLID, popID);
      creal vxBlock = block_parameters[BlockParams::VXCRD];
      creal vyBlock = block_parameters[BlockParams::VYCRD];
      creal vzBlock = block_parameters[BlockParams::VZCRD];
      creal dvx = block_parameters[BlockParams::DVX];
      creal dvy = block_parameters[BlockParams::DVY];
      creal dvz = block_parameters[BlockParams::DVZ];

      // Calculate volume-averaged distribution function for each cell in the block.
      for (uint kc = 0; kc < WID; ++kc) {
         creal vzC = vzBlock + kc * dvz;
         for (uint jc = 0; jc < WID; ++jc) {
            creal vyC = vyBlock + jc * dvy;
            for (uint ic = 0; ic < WID; ++ic) {
               creal vxC = vxBlock + ic * dvx;
               // hyzhou: most values are very close to 0!
               f[blockLID * WID3 + cellIndex(ic, jc, kc)] =
                   maxwellian(popID, rho, T, vxC + 0.5 * dvx - Vx, vyC + 0.5 * dvy - Vy, vzC + 0.5 * dvz - Vz);
            }
         }
      } // loop over cells in velocity block
   }    // loop over velocity blocks

   // For all species
   calculateCellMoments(cell, true);

   cell->get_population(popID).RHO_R = cell->get_population(popID).RHO;
   cell->get_population(popID).RHO_V = cell->get_population(popID).RHO;
   for (int idim = 0; idim < 3; idim++) {
      cell->get_population(popID).V_R[idim] = cell->get_population(popID).V[idim];
      cell->get_population(popID).V_V[idim] = cell->get_population(popID).V[idim];
      cell->get_population(popID).P_R[idim] = cell->get_population(popID).P[idim];
      cell->get_population(popID).P_V[idim] = cell->get_population(popID).P[idim];
   }

   // Store blocks according to the sparsity threshold
   cell->adjustSingleCellVelocityBlocks(popID);
}

// Find velocity blocks to initialize.
std::vector<vmesh::GlobalID> Wave::findBlocksToInit(const uint popID, spatial_cell::SpatialCell* cell, creal rho,
                                                    creal T, creal Vx0, creal Vy0, creal Vz0) {
   std::vector<vmesh::GlobalID> blocksToInit;
   bool search = true;
   uint counter = 0;
   const uint8_t refLevel = 0;

   const vmesh::LocalID* vblocks_ini = cell->get_velocity_grid_length(popID, refLevel);
   const Real dvx = cell->get_velocity_grid_cell_size(popID, refLevel)[0];
   const Real dvy = cell->get_velocity_grid_cell_size(popID, refLevel)[1];
   const Real dvz = cell->get_velocity_grid_cell_size(popID, refLevel)[2];

   while (search) {
#warning TODO: add SpatialCell::getVelocityBlockMinValue() in place of sparseMinValue?
      if (0.1 * getObjectWrapper().particleSpecies[popID].sparseMinValue >
              maxwellian(popID, rho, T, counter * cell->get_velocity_grid_block_size(popID, refLevel)[0] + 0.5 * dvx,
                         0.5 * dvy, 0.5 * dvz) ||
          counter > vblocks_ini[0]) {
         search = false;
      }
      counter++;
   }
   counter += 2;

   Real vRadiusSquared = (Real)counter * (Real)counter * cell->get_velocity_grid_block_size(popID, refLevel)[0] *
                         cell->get_velocity_grid_block_size(popID, refLevel)[0];

   for (uint kv = 0; kv < vblocks_ini[2]; ++kv)
      for (uint jv = 0; jv < vblocks_ini[1]; ++jv)
         for (uint iv = 0; iv < vblocks_ini[0]; ++iv) {
            vmesh::GlobalID blockIndices[3];
            blockIndices[0] = iv;
            blockIndices[1] = jv;
            blockIndices[2] = kv;
            const vmesh::GlobalID blockGID = cell->get_velocity_block(popID, blockIndices, refLevel);

            Real V_crds[3];
            cell->get_velocity_block_coordinates(popID, blockGID, V_crds);
            Real dV[3];
            cell->get_velocity_block_size(popID, blockGID, dV);
            V_crds[0] += 0.5 * dV[0];
            V_crds[1] += 0.5 * dV[1];
            V_crds[2] += 0.5 * dV[2];
            Real R2 = ((V_crds[0] - Vx0) * (V_crds[0] - Vx0) + (V_crds[1] - Vy0) * (V_crds[1] - Vy0) +
                       (V_crds[2] - Vz0) * (V_crds[2] - Vz0));

            if (R2 < vRadiusSquared) {
               cell->add_velocity_block(blockGID, popID);
               blocksToInit.push_back(blockGID);
            }
         }

   return blocksToInit;
}

// Set the non-dipole magnetic field boundary condition.
Real Wave::fieldSolverBoundaryCondMagneticField(
    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& bGrid,
    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, cint i, cint j, cint k, cuint component) {
   typedef Parameters P;
   Real dB = 0.0;
   creal dx = P::dx_ini;
   creal dy = P::dy_ini;
   creal dz = P::dz_ini;
   const array<int, 3> globalIndices = technicalGrid.getGlobalIndices(i, j, k);
   creal x = (convert<Real>(globalIndices[0]) + 0.5) * technicalGrid.DX + P::xmin;
   creal y = (convert<Real>(globalIndices[1]) + 0.5) * technicalGrid.DY + P::ymin;
   creal z = (convert<Real>(globalIndices[2]) + 0.5) * technicalGrid.DZ + P::zmin;

   auto coords = bGrid.getPhysicalCoords(i, j, k);
   coords[0] += 0.5 * bGrid.DX;
   coords[1] += 0.5 * bGrid.DY;
   coords[2] += 0.5 * bGrid.DZ;

   bool isCellOnFace[6];
   determineFace(isCellOnFace, x, y, z, dx, dy, dz, true);

   if (wavetype == "None") {
      for (uint iface = 0; iface < 6; iface++) {
         if (isCellOnFace[iface]) {
            Real tmp, t;
            t = P::t + 0.5 * P::dt;
            tmp = std::sin(2 * PI * freq * t + waveVector[0] * x + waveVector[1] * y + waveVector[2] * z);
            dB = B0[component] + B1[component] * tmp;
            break; // This effectively sets the precedence of faces through the order of faces.
         }
      }
   } else if (wavetype == "Alfven" || wavetype == "Fast") {
      for (uint iface = 0; iface < 6; iface++) {
         if (isCellOnFace[iface]) {
            Real E[3], omega;
            // E = B0 x v
            E[0] = B0[1] * V1[2] - B0[2] * V1[1];
            E[1] = -(B0[0] * V1[2] - B0[2] * V1[0]);
            E[2] = B0[0] * V1[1] - B0[1] * V1[0];
            // B1 = k x E / omega
            omega = 2 * PI * freq;
            B1[0] = (waveVector[1] * E[2] - waveVector[2] * E[1]) / omega;
            B1[1] = -(waveVector[0] * E[2] - waveVector[2] * E[0]) / omega;
            B1[2] = (waveVector[0] * E[1] - waveVector[1] * E[0]) / omega;
            dB = B0[component] + B1[component];
            break; // This effectively sets the precedence of faces through the order of faces.
         }
      }
   }

   return dB;
}

void Wave::fieldSolverBoundaryCondMagneticFieldProjection(
    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& bGrid,
    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid, cint i, cint j, cint k) {}

void Wave::fieldSolverBoundaryCondElectricField(
    FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH>& EGrid, cint i, cint j, cint k,
    cuint component) {
   EGrid.get(i, j, k)->at(fsgrids::efield::EX + component) = 0.0;
}

void Wave::fieldSolverBoundaryCondHallElectricField(
    FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH>& EHallGrid, cint i, cint j, cint k,
    cuint component) {
   std::array<Real, fsgrids::ehall::N_EHALL>* cp = EHallGrid.get(i, j, k);
   switch (component) {
   case 0:
      cp->at(fsgrids::ehall::EXHALL_000_100) = 0.0;
      cp->at(fsgrids::ehall::EXHALL_010_110) = 0.0;
      cp->at(fsgrids::ehall::EXHALL_001_101) = 0.0;
      cp->at(fsgrids::ehall::EXHALL_011_111) = 0.0;
      break;
   case 1:
      cp->at(fsgrids::ehall::EYHALL_000_010) = 0.0;
      cp->at(fsgrids::ehall::EYHALL_100_110) = 0.0;
      cp->at(fsgrids::ehall::EYHALL_001_011) = 0.0;
      cp->at(fsgrids::ehall::EYHALL_101_111) = 0.0;
      break;
   case 2:
      cp->at(fsgrids::ehall::EZHALL_000_001) = 0.0;
      cp->at(fsgrids::ehall::EZHALL_100_101) = 0.0;
      cp->at(fsgrids::ehall::EZHALL_010_011) = 0.0;
      cp->at(fsgrids::ehall::EZHALL_110_111) = 0.0;
      break;
   default:
      abort_mpi("Invalid component", 1);
   }
}

void Wave::fieldSolverBoundaryCondGradPeElectricField(
    FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, 2>& EGradPeGrid, cint i, cint j, cint k, cuint component) {
   EGradPeGrid.get(i, j, k)->at(fsgrids::egradpe::EXGRADPE + component) = 0.0;
}

void Wave::fieldSolverBoundaryCondDerivatives(
    FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH>& dPerBGrid,
    FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH>& dMomentsGrid, cint i, cint j, cint k,
    cuint RKCase, cuint component) {
   // this->setCellDerivativesToZero(dPerBGrid, dMomentsGrid, i, j, k, component);
   // hyzhou: do I need to set the derivatives? What are they used for?
}

void Wave::fieldSolverBoundaryCondBVOLDerivatives(
    FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH>& volGrid, cint i, cint j, cint k,
    cuint component) {
   // this->setCellBVOLDerivativesToZero(volGrid, i, j, k, component);
   // hyzhou: do I need to set the derivatives? What are they used for?
}

void Wave::getFaces(std::array<bool, 6>& faces) { faces = facesToProcess; }

string Wave::getName() const { return "Wave"; }
uint Wave::getIndex() const { return sysboundarytype::WAVE; }

} // namespace SBC
