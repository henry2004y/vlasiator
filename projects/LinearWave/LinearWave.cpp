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

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "LinearWave.h"

using namespace spatial_cell;

namespace projects {
LinearWave::LinearWave() : TriAxisSearch() {}
LinearWave::~LinearWave() {}

bool LinearWave::initialize(void) {
   bool success = Project::initialize();

   creal m = physicalconstants::MASS_PROTON;
   creal e = physicalconstants::CHARGE;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;
   creal mu0 = physicalconstants::MU_0;

   creal p0 = 1.0 / gamma;

   // Three reference values: BRef, nRef, mRef=m
   BRef = 1e-9;
   nRef = 1e6;
   // All other reference values can be derived.
   lRef = std::cbrt(nRef);
   uRef = BRef / sqrt(mu0 * m * nRef);
   tRef = lRef / uRef;
   PRef = 0.5 * m * nRef * sqr(uRef);
   // In the hybrid model, electrons are massless so there is no contribution to pressure.
   // Maybe I shouldn't include the factor 2 here?
   TRef = PRef / (kB * nRef);
   // TRef = PRef / (2 * kB * nRef);
   // TRef = sqr(BRef) / (2 * mu0 * kB * nRef);
   // lRef = P::xmax - P::xmin;
   // tRef = P::t_max;

   Real Vth = 0.0;
   // Defined in dimensionless units
   if (wavetype == "fast") {
      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));
      rho1 = 4.472135954999580e-1;
      rhoux1 = -8.944271909999160e-1;
      rhouy1 = 4.216370213557840e-1;
      rhouz1 = 1.490711984999860e-1;
      E1 = 2.012457825664615e+0;
      Bx1 = 0.0;
      By1 = 8.432740427115680e-1;
      Bz1 = 2.981423969999720e-1;
   } else if (wavetype == "alfven") {
      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));
      rho1 = 0.0;
      rhoux1 = 0.0;
      rhouy1 = -3.333333333333333e-1;
      rhouz1 = 9.428090415820634e-1;
      E1 = 0.0;
      Bx1 = 0.0;
      By1 = -3.333333333333333e-1;
      Bz1 = 9.428090415820634e-1;
   } else if (wavetype == "slow") {
      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));
      rho1 = 8.944271909999159e-1;
      rhoux1 = -4.472135954999579e-1;
      rhouy1 = -8.432740427115680e-1;
      rhouz1 = -2.981423969999720e-1;
      E1 = 6.708136850795449e-01;
      Bx1 = 0.0;
      By1 = -4.216370213557841e-1;
      Bz1 = -1.490711984999860e-1;
   } else if (wavetype == "alfvenx") {
      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 0.0;
      Bz0 = 0.0;
      E0 = p0 / (gamma - 1) + 0.5 * sqr(Bx0);
      rho1 = 0.0;
      rhoux1 = 0.0;
      rhouy1 = -3.333333333333333e-1;
      rhouz1 = 0.0;
      E1 = 0.0;
      Bx1 = 0.0;
      By1 = -3.333333333333333e-1;
      Bz1 = 0.0;
   } else if (wavetype == "nowave"){
      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));
      rho1 = 0.0;
      rhoux1 = 0.0;
      rhouy1 = 0.0;
      rhouz1 = 0.0;
      E1 = 0.0;
      Bx1 = 0.0;
      By1 = 0.0;
      Bz1 = 0.0;
   } else if (wavetype == "nobwave"){
      rho0 = 1.0;
      Bx0 = 0.0;
      By0 = 0.0;
      Bz0 = 0.0;
      E0 = p0 / (gamma - 1);
      rho1 = 0.0;
      rhoux1 = 0.0;
      rhouy1 = 0.0;
      rhouz1 = 0.0;
      E1 = 0.0;
      Bx1 = 0.0;
      By1 = 0.0;
      Bz1 = 0.0;
   }

   Vth = sqrt(gamma * p0 / rho0) * uRef;

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         std::cout << "Unit conversion factors from dimensionless to SI:" << std::endl;
         std::cout << "lRef: " << lRef << " [m]" << std::endl;
         std::cout << "tRef: " << tRef << " [s]" << std::endl;
         std::cout << "vRef: " << uRef << " [m/s]" << std::endl;
         std::cout << "nRef: " << nRef << " [#]" << std::endl;
         std::cout << "TRef: " << TRef << " [K]" << std::endl;
         std::cout << "PRef: " << PRef << " [Pa]" << std::endl;
         std::cout << "BRef: " << BRef << " [T]" << std::endl;
         // Useful for setting the Vspace ~ 10*Vth width
         std::cout << "Vth: " << Vth << " [m/s]" << std::endl;
         const vmesh::MeshParameters& mp = getObjectWrapper().velocityMeshes[0];
         Real dvx = (mp.meshLimits[1] - mp.meshLimits[0]) / (mp.gridLength[0] * mp.blockLength[0]);
         std::cout << "Vth/dvx = " << Vth / dvx << std::endl;
         std::cout << "Vth/Vmax = " << Vth / mp.meshLimits[1] << std::endl;
         if (dvx > Vth) {
            std::cout << "vmesh maybe too coarse!" << dvx << std::endl;
         }
         if (Vth > mp.meshLimits[1]) {
            std::cout << "vmesh maybe too small!" << dvx << std::endl;
         }
      }
   }

   return success;
}

void LinearWave::addParameters() {
   typedef Readparameters RP;
   // Dimensionless inputs
   RP::add("LinearWave.Bx0", "Background x-magnetic field", 1.0);
   RP::add("LinearWave.By0", "Background y-magnetic field", 1.4142135623730951);
   RP::add("LinearWave.Bz0", "Background z-magnetic field", 0.5);
   RP::add("LinearWave.Bx1", "Perturbation x-magnetic field", 0.0);
   RP::add("LinearWave.By1", "Perturbation y-magnetic field", 0.0);
   RP::add("LinearWave.Bz1", "Perturbation z-magnetic field", 0.0);
   RP::add("LinearWave.Amplitude", "Amplitude of perturbations", 1e-6);
   RP::add("LinearWave.verbose", "Turn on/off detailed information", 0);
   RP::add("LinearWave.wavetype", "Type of MHD waves with preset parameters (\"fast\", \"alfven\", \"slow\", \"none\")",
           "none");
   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      RP::add(pop + "_LinearWave.rho0", "Background mass density", 1.0);
      RP::add(pop + "_LinearWave.E0", "Background energy density", 0.6);
      // Perturbation right-eigenvectors
      RP::add(pop + "_LinearWave.rho1", "Perturbation mass density", 0.0);
      RP::add(pop + "_LinearWave.rhoux1", "Perturbation x-velocity", 0.0);
      RP::add(pop + "_LinearWave.rhouy1", "Perturbation y-velocity", 0.0);
      RP::add(pop + "_LinearWave.rhouz1", "Perturbation z-velocity", 0.0);
      RP::add(pop + "_LinearWave.E1", "Perturbation energy density", 0.0);
   }
}

void LinearWave::getParameters() {
   Project::getParameters();

   typedef Readparameters RP;
   RP::get("LinearWave.Bx0", Bx0);
   RP::get("LinearWave.By0", By0);
   RP::get("LinearWave.Bz0", Bz0);
   RP::get("LinearWave.Bx1", Bx1);
   RP::get("LinearWave.By1", By1);
   RP::get("LinearWave.Bz1", Bz1);
   RP::get("LinearWave.Amplitude", A);
   RP::get("LinearWave.verbose", verbose);
   RP::get("LinearWave.wavetype", wavetype);

   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;

      if (pop != "proton") {
         int myRank;
         MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
         if (myRank == MASTER_RANK) {
            std::cerr << "Linear wave test only supports proton!" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
         }
      }

      RP::get(pop + "_LinearWave.rho0", rho0);
      RP::get(pop + "_LinearWave.E0", E0);
      RP::get(pop + "_LinearWave.rho1", rho1);
      RP::get(pop + "_LinearWave.rhoux1", rhoux1);
      RP::get(pop + "_LinearWave.rhouy1", rhouy1);
      RP::get(pop + "_LinearWave.rhouz1", rhouz1);
      RP::get(pop + "_LinearWave.E1", E1);
   }
}

Real LinearWave::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                               creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;

   creal perturb = A * sin(2 * M_PI * x / lRef);
   creal drho = perturb * rho1;
   creal rho = rho0 + drho;
   Real n = rho;
   // Assume background velocity 0
   Real ux = perturb * rhoux1 / rho;
   Real uy = perturb * rhouy1 / rho;
   Real uz = perturb * rhouz1 / rho;

   creal p = (gamma - 1) * ((E0 + E1 * perturb) - 0.5 * rho * (sqr(ux) + sqr(uy) + sqr(uz)) -
                            0.5 * (sqr(Bx0 + Bx1 * perturb) + sqr(By0 + By1 * perturb) + sqr(Bz0 + Bz1 * perturb)));
   Real T = p / rho;

   // Convert to SI units
   n *= nRef;
   ux *= uRef;
   uy *= uRef;
   uz *= uRef;
   T *= TRef;

   creal coef = m / (2 * M_PI * kB * T);
   // Maxwellian f(v) = f(v;n,u,T)
   creal f = n * sqrt(coef) * coef * exp(-coef * M_PI * (sqr(vx - ux) + sqr(vy - uy) + sqr(vz - uz)));

   return f;
}

Real LinearWave::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                       creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                       const uint popID) const {
   Real f = 0.0;
   f = getMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx, dvy,
                     dvz, popID);

   return f;
}

std::vector<std::array<Real, 3>> LinearWave::getV0(creal x, creal y, creal z, const uint popID) const {
   std::vector<std::array<Real, 3>> V0;
   std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
   V0.push_back(v);
   return V0;
}

void LinearWave::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

void LinearWave::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                  FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                  FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   // Background field
   ConstantField bgField;
   bgField.initialize(Bx0 * BRef, By0 * BRef, Bz0 * BRef);
   setBackgroundField(bgField, BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
      for (int i = 0; i < localSize[0]; ++i) {
         for (int j = 0; j < localSize[1]; ++j) {
            for (int k = 0; k < localSize[2]; ++k) {
               const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
               std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);
               creal perturb = A * sin(2 * M_PI * x[0] / lRef);
               cell->at(fsgrids::bfield::PERBX) = perturb * Bx1 * BRef;
               cell->at(fsgrids::bfield::PERBY) = perturb * By1 * BRef;
               cell->at(fsgrids::bfield::PERBZ) = perturb * Bz1 * BRef;
            }
         }
      }
   }
}

} // namespace projects
