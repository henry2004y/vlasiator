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

   RP::add("Slab.bx0L", "Left domain B field value", 0.0);
   RP::add("Slab.by0L", "Left domain B field value", 0.0);
   RP::add("Slab.bz0L", "Left domain B field value", 1.0);
   RP::add("Slab.bx0R", "Right domain B field value", 0.0);
   RP::add("Slab.by0R", "Right domain B field value", 0.0);
   RP::add("Slab.bz0R", "Right domain B field value", 1.0);
   RP::add("Slab.halfWidth", "Halfwidth of the current sheet (m)", 0.0);
   RP::add("Slab.Tratio", "Te/Ti for setting B field from pressure balance", 0.4);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      RP::add(pop + "_Slab.vx0L", "Left bulk velocity in x", 0.0);
      RP::add(pop + "_Slab.vy0L", "Left bulk velocity in y", 0.0);
      RP::add(pop + "_Slab.vz0L", "Left bulk velocuty in z", 0.0);
      RP::add(pop + "_Slab.nL", "Left number density", 10.0);
      RP::add(pop + "_Slab.TL", "Left temperature (K)", 2.0e6);

      RP::add(pop + "_Slab.vx0R", "Right bulk velocity in x", 0.0);
      RP::add(pop + "_Slab.vy0R", "Right bulk velocity in y", 0.0);
      RP::add(pop + "_Slab.vz0R", "Right bulk velocuty in z", 0.0);
      RP::add(pop + "_Slab.nR", "Right number density", 1.0);
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
   RP::get("Slab.Tratio", this->Tratio); // redundant in field solver!

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
   const Real xmid = 0.5 * (P::xmax + P::xmin);
   Real vx, vy, vz;
   if (x < xmid) {
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
   creal mass = getObjectWrapper().particleSpecies[popID].mass;
   creal KB = physicalconstants::K_B;
   const SlabSpeciesParameters& sP = this->speciesParams[popID];

   Real f, n, T;
   creal xmid = 0.5 * (P::xmax + P::xmin);
   creal(&v0)[3] = x < xmid ? sP.v0L : sP.v0R;
   if (x < xmid) {
      n = sP.nL;
      T = sP.TL;
   } else {
      n = sP.nR;
      T = sP.TR;
   }
   if (halfWidth != 0) {
      n = 0.5 * ((sP.nL + sP.nR) + (sP.nR - sP.nL) * tanh((x - xmid) / halfWidth));
   }

   creal coef1 = 2.0 * KB * T;
   creal coef2 = mass / (M_PI * coef1);
   f = n * coef2 * sqrt(coef2) * exp(-mass * (sqr(vx - v0[0]) + sqr(vy - v0[1]) + sqr(vz - v0[2])) / coef1);

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
   // It would be better to put all B field into the perturbed term since it's not uniform at t = 0!
   setBackgroundFieldToZero(BgBGrid);

   const Real xmid = 0.5 * (P::xmax + P::xmin);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

      if (this->halfWidth != 0) {
#pragma omp parallel for collapse(3)
         for (int i = 0; i < localSize[0]; ++i) {
            for (int j = 0; j < localSize[1]; ++j) {
               for (int k = 0; k < localSize[2]; ++k) {
                  const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);

                  // Assume only one ion species
                  const SlabSpeciesParameters& sP = this->speciesParams[0];
                  Real KB = physicalconstants::K_B;
                  Real MU_0 = physicalconstants::MU_0;
                  Real n, Pi, Pe, PB, Bz;

                  n = 0.5 * ((sP.nL + sP.nR) + (sP.nR - sP.nL) * tanh((x[0] - xmid) / halfWidth));

                  Pi = sP.nL * KB * sP.TL;
                  Pe = sP.nL * KB * sP.TL * this->Tratio;
                  PB = Pi + Pe + sqr(this->b0L[2]) / (2 * MU_0) - n * KB * sP.TL * (1 + this->Tratio);

                  Bz = std::sqrt(PB * 2 * MU_0);
                  // hyzhou: currently limited to Bz only!
                  cell->at(fsgrids::bfield::PERBX) = 0.0;
                  cell->at(fsgrids::bfield::PERBY) = 0.0;
                  cell->at(fsgrids::bfield::PERBZ) = Bz;
               }
            }
         }
      } else {
#pragma omp parallel for collapse(3)
         for (int i = 0; i < localSize[0]; ++i) {
            for (int j = 0; j < localSize[1]; ++j) {
               for (int k = 0; k < localSize[2]; ++k) {
                  const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);

                  if (x[0] < xmid) {
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
}

} // namespace projects
