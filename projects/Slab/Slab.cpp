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

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

%
Interplanetary shock project by Markus Battarbee (markus.battarbee@gmail.com)
Based on SilvaShock project by Urs Ganse
Previous development version name was UtuShock
*/

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../../backgroundfield/backgroundfield.h"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "Slab.h"

using namespace spatial_cell;

namespace projects {
Slab::Slab() : TriAxisSearch() {}
Slab::~Slab() {}

bool Slab::initialize() { return Project::initialize(); }

void Slab::addParameters() {
   typedef Readparameters RP;

   RP::add("Slab.bx0L", "Left domain B field value (T)", 1.0e-9);
   RP::add("Slab.by0L", "Left domain B field value (T)", 2.0e-9);
   RP::add("Slab.bz0L", "Left domain B field value (T)", 3.0e-9);
   RP::add("Slab.bx0R", "Right domain B field value (T)", 1.0e-9);
   RP::add("Slab.by0R", "Right domain B field value (T)", 2.0e-9);
   RP::add("Slab.bz0R", "Right domain B field value (T)", 3.0e-9);
   RP::add("Slab.halfWidth", "Halfwidth of the current sheet (m)", 0.0);

   RP::add("Slab.AMR_L1width", "L1 AMR region width (m)", 0);
   RP::add("Slab.AMR_L2width", "L2 AMR region width (m)", 0);
   RP::add("Slab.AMR_L3width", "L3 AMR region width (m)", 0);
   RP::add("Slab.AMR_L4width", "L4 AMR region width (m)", 0);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      RP::add(pop + "_Slab.vx0L", "Left bulk velocity in x", 0.0);
      RP::add(pop + "_Slab.vy0L", "Left bulk velocity in y", 0.0);
      RP::add(pop + "_Slab.vz0L", "Left bulk velocuty in z", 0.0);
      RP::add(pop + "_Slab.nL", "Left number density (m^-3)", 1.0e7);
      RP::add(pop + "_Slab.TL", "Left temperature (K)", 2.0e6);

      RP::add(pop + "_Slab.vx0R", "Right bulk velocity in x", 0.0);
      RP::add(pop + "_Slab.vy0R", "Right bulk velocity in y", 0.0);
      RP::add(pop + "_Slab.vz0R", "Right bulk velocuty in z", 0.0);
      RP::add(pop + "_Slab.nR", "Right number density (m^-3)", 1.0e7);
      RP::add(pop + "_Slab.TR", "Right temperature (K)", 2.0e6);

      RP::add(pop + "_Slab.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add(pop + "_Slab.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
   }
}

void Slab::getParameters() {
   Project::getParameters();

   typedef Readparameters RP;
   RP::get("Slab.bx0L", this->b0L[0]);
   RP::get("Slab.by0L", this->b0L[1]);
   RP::get("Slab.bz0L", this->b0L[2]);
   RP::get("Slab.bx0R", this->b0R[0]);
   RP::get("Slab.by0R", this->b0R[1]);
   RP::get("Slab.bz0R", this->b0R[2]);
   RP::get("Slab.halfWidth", this->halfWidth);

   RP::get("Slab.AMR_L1width", this->AMR_L1width);
   RP::get("Slab.AMR_L2width", this->AMR_L2width);
   RP::get("Slab.AMR_L3width", this->AMR_L3width);
   RP::get("Slab.AMR_L4width", this->AMR_L4width);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      SlabSpeciesParameters sP;

      RP::get(pop + "_Slab.vx0L", sP.v0L[0]);
      RP::get(pop + "_Slab.vy0L", sP.v0L[1]);
      RP::get(pop + "_Slab.vz0L", sP.v0L[2]);
      RP::get(pop + "_Slab.nL", sP.nL);
      RP::get(pop + "_Slab.TL", sP.TL);

      RP::get(pop + "_Slab.vx0R", sP.v0R[0]);
      RP::get(pop + "_Slab.vy0R", sP.v0R[1]);
      RP::get(pop + "_Slab.vz0R", sP.v0R[2]);
      RP::get(pop + "_Slab.nR", sP.nR);
      RP::get(pop + "_Slab.TR", sP.TR);

      RP::get(pop + "_Slab.nSpaceSamples", sP.nSpaceSamples);
      RP::get(pop + "_Slab.nVelocitySamples", sP.nVelocitySamples);

      speciesParams.push_back(sP);
   }
}

std::vector<std::array<Real, 3>> Slab::getV0(creal x, creal y, creal z, const uint popID) const {
   Real mass = getObjectWrapper().particleSpecies[popID].mass;
   Real mu0 = physicalconstants::MU_0;
   const SlabSpeciesParameters& sP = this->speciesParams[popID];

   Real vx, vy, vz;
   if (x < 0) {
      vx = sP.v0L[0];
      vy = sP.v0L[1];
      vz = sP.v0L[2];
   } else {
      vx = sP.v0R[0];
      vy = sP.v0R[1];
      vz = sP.v0R[2];
   }

   std::array<Real, 3> V0{{vx, vy, vz}};
   std::vector<std::array<Real, 3>> retval;
   retval.push_back(V0);

   return retval;
}

Real Slab::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                           creal& dvz, const uint popID) const {
   Real mass = getObjectWrapper().particleSpecies[popID].mass;
   Real KB = physicalconstants::K_B;
   const SlabSpeciesParameters& sP = this->speciesParams[popID];

   Real f, n, T;
   creal(&v0)[3] = x < 0 ? sP.v0L : sP.v0R;
   if (x < 0) {
      n = sP.nL;
      T = sP.TL;
   } else {
      n = sP.nR;
      T = sP.TR;
   }
   if (halfWidth != 0) {
      n = 0.5 * ((sP.nL + sP.nR) + (sP.nL - sP.nR) * tanh(x / halfWidth));
   }

   f = n * std::pow(mass / (2.0 * M_PI * KB * T), 1.5) *
       exp(-mass * ((vx - v0[0]) * (vx - v0[0]) + (vy - v0[1]) * (vy - v0[1]) + (vz - v0[2]) * (vz - v0[2])) /
           (2.0 * KB * T));

   return f;
}

Real Slab::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy,
                                 creal& vz, creal& dvx, creal& dvy, creal& dvz, const uint popID) const {
   const SlabSpeciesParameters& sP = this->speciesParams[popID];
   Real f = 0.0;
   if ((sP.nSpaceSamples > 1) && (sP.nVelocitySamples > 1)) {
      creal d_x = dx / (sP.nSpaceSamples - 1);
      creal d_y = dy / (sP.nSpaceSamples - 1);
      creal d_z = dz / (sP.nSpaceSamples - 1);
      creal d_vx = dvx / (sP.nVelocitySamples - 1);
      creal d_vy = dvy / (sP.nVelocitySamples - 1);
      creal d_vz = dvz / (sP.nVelocitySamples - 1);

      Real avg = 0.0;

      for (uint i = 0; i < sP.nSpaceSamples; ++i)
         for (uint j = 0; j < sP.nSpaceSamples; ++j)
            for (uint k = 0; k < sP.nSpaceSamples; ++k)
               for (uint vi = 0; vi < sP.nVelocitySamples; ++vi)
                  for (uint vj = 0; vj < sP.nVelocitySamples; ++vj)
                     for (uint vk = 0; vk < sP.nVelocitySamples; ++vk) {
                        avg += getDistribValue(x + i * d_x, y + j * d_y, z + k * d_z, vx + vi * d_vx, vy + vj * d_vy,
                                               vz + vk * d_vz, dvx, dvy, dvz, popID);
                     }

      f = avg / (sP.nSpaceSamples * sP.nSpaceSamples * sP.nSpaceSamples) /
          (sP.nVelocitySamples * sP.nVelocitySamples * sP.nVelocitySamples);
   } else {
      f = getDistribValue(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx,
                          dvy, dvz, popID);
   }

   return f;
}

void Slab::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

void Slab::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                            FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                            FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   setBackgroundFieldToZero(BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
      for (int x = 0; x < localSize[0]; ++x) {
         for (int y = 0; y < localSize[1]; ++y) {
            for (int z = 0; z < localSize[2]; ++z) {
               const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(x, y, z);
               std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(x, y, z);

               if (xyz[0] < 0.0) {
                  cell->at(fsgrids::bfield::PERBX) = this->b0L[0];
                  cell->at(fsgrids::bfield::PERBY) = this->b0L[1];
                  cell->at(fsgrids::bfield::PERBZ) = this->b0L[2];
               } else {
                  cell->at(fsgrids::bfield::PERBX) = this->b0R[0];
                  cell->at(fsgrids::bfield::PERBY) = this->b0R[1];
                  cell->at(fsgrids::bfield::PERBZ) = this->b0R[2];
               }
            }
         }
      }
   }
}

bool Slab::refineSpatialCells(dccrg::Dccrg<spatial_cell::SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid) const {
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   std::vector<CellID> refinedCells;

   if (myRank == MASTER_RANK)
      std::cout << "Maximum refinement level is " << mpiGrid.mapping.get_maximum_refinement_level() << std::endl;

   // Calculate regions for refinement
   if (P::amrMaxSpatialRefLevel > 0) {
      // L1 refinement.
      for (uint i = 0; i < P::xcells_ini; ++i) {
         for (uint j = 0; j < P::ycells_ini; ++j) {
            for (uint k = 0; k < P::zcells_ini; ++k) {
               std::array<double, 3> xyz;
               xyz[0] = P::xmin + (i + 0.5) * P::dx_ini;
               xyz[1] = P::ymin + (j + 0.5) * P::dy_ini;
               xyz[2] = P::zmin + (k + 0.5) * P::dz_ini;

               if (abs(xyz[0]) < AMR_L1width) {
                  CellID myCell = mpiGrid.get_existing_cell(xyz);
                  mpiGrid.refine_completely(myCell);
               }
            }
         }
      }
      refinedCells = mpiGrid.stop_refining(true);
      if (myRank == MASTER_RANK)
         std::cout << "Finished first level of refinement" << std::endl;
      mpiGrid.balance_load();
   }

   if (P::amrMaxSpatialRefLevel > 1) {
      // L2 refinement.
      for (uint i = 0; i < 2 * P::xcells_ini; ++i) {
         for (uint j = 0; j < 2 * P::ycells_ini; ++j) {
            for (uint k = 0; k < 2 * P::zcells_ini; ++k) {

               std::array<double, 3> xyz;
               xyz[0] = P::xmin + (i + 0.5) * 0.5 * P::dx_ini;
               xyz[1] = P::ymin + (j + 0.5) * 0.5 * P::dy_ini;
               xyz[2] = P::zmin + (k + 0.5) * 0.5 * P::dz_ini;

               if (abs(xyz[0]) < AMR_L2width) {
                  CellID myCell = mpiGrid.get_existing_cell(xyz);
                  mpiGrid.refine_completely(myCell);
               }
            }
         }
      }
      refinedCells = mpiGrid.stop_refining(true);
      if (myRank == MASTER_RANK)
         std::cout << "Finished second level of refinement" << std::endl;
      mpiGrid.balance_load();
   }

   if (P::amrMaxSpatialRefLevel > 2) {
      // L3 refinement.
      for (uint i = 0; i < 4 * P::xcells_ini; ++i) {
         for (uint j = 0; j < 4 * P::ycells_ini; ++j) {
            for (uint k = 0; k < 4 * P::zcells_ini; ++k) {

               std::array<double, 3> xyz;
               xyz[0] = P::xmin + (i + 0.5) * 0.25 * P::dx_ini;
               xyz[1] = P::ymin + (j + 0.5) * 0.25 * P::dy_ini;
               xyz[2] = P::zmin + (k + 0.5) * 0.25 * P::dz_ini;

               if (abs(xyz[0]) < AMR_L3width) {
                  CellID myCell = mpiGrid.get_existing_cell(xyz);
                  mpiGrid.refine_completely(myCell);
               }
            }
         }
      }
      refinedCells = mpiGrid.stop_refining(true);
      if (myRank == MASTER_RANK)
         std::cout << "Finished third level of refinement" << std::endl;
      mpiGrid.balance_load();
   }

   if (P::amrMaxSpatialRefLevel > 3) {
      // L4 refinement.
      for (uint i = 0; i < 8 * P::xcells_ini; ++i) {
         for (uint j = 0; j < 8 * P::ycells_ini; ++j) {
            for (uint k = 0; k < 8 * P::zcells_ini; ++k) {

               std::array<double, 3> xyz;
               xyz[0] = P::xmin + (i + 0.5) * 0.125 * P::dx_ini;
               xyz[1] = P::ymin + (j + 0.5) * 0.125 * P::dy_ini;
               xyz[2] = P::zmin + (k + 0.5) * 0.125 * P::dz_ini;

               if (abs(xyz[0]) < AMR_L4width) {
                  CellID myCell = mpiGrid.get_existing_cell(xyz);
                  mpiGrid.refine_completely(myCell);
               }
            }
         }
      }
      refinedCells = mpiGrid.stop_refining(true);
      if (myRank == MASTER_RANK)
         std::cout << "Finished fourth level of refinement" << std::endl;
      mpiGrid.balance_load();
   }

   return true;
}

} // namespace projects
