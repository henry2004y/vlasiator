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

   // Three reference values: BRef, nRef, mRef=m
   // All other reference values can be derived.
   lRef = sqrt(m / (mu0 * nRef)) / e;
   uRef = BRef / sqrt(mu0 * m * nRef);
   tRef = lRef / uRef;
   pRef = 0.5 * m * nRef * sqr(uRef);
   TRef = pRef / (kB * nRef);

   // Defined in dimensionless units
   if (wavetype == "fast") {
      isBiMaxwellian = false;
      isPrimitive = false;

      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      p0 = 1.0 / gamma;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));

      rho1 = A * 4.472135954999580e-1;
      rhoux1 = A * -8.944271909999160e-1;
      rhouy1 = A * 4.216370213557840e-1;
      rhouz1 = A * 1.490711984999860e-1;
      E1 = A * 2.012457825664615e+0;
      Bx1 = A * 0.0;
      By1 = A * 8.432740427115680e-1;
      Bz1 = A * 2.981423969999720e-1;
   } else if (wavetype == "fastr") {
      isBiMaxwellian = false;
      isPrimitive = false;

      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      p0 = 1.0 / gamma;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));

      rho1 = A * 4.472135954999580e-1;
      rhoux1 = A * 8.944271909999160e-1;
      rhouy1 = A * -4.216370213557840e-1;
      rhouz1 = A * -1.490711984999860e-1;
      E1 = A * 2.012457825664615e+0;
      Bx1 = A * 0.0;
      By1 = A * 8.432740427115680e-1;
      Bz1 = A * 2.981423969999720e-1;
   } else if (wavetype == "alfven") {
      isBiMaxwellian = false;
      isPrimitive = false;

      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      p0 = 1.0 / gamma;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));

      rho1 = A * 0.0;
      rhoux1 = A * 0.0;
      rhouy1 = A * -3.333333333333333e-1;
      rhouz1 = A * 9.428090415820634e-1;
      E1 = A * 0.0;
      Bx1 = A * 0.0;
      By1 = A * -3.333333333333333e-1;
      Bz1 = A * 9.428090415820634e-1;
   } else if (wavetype == "slow") {
      isBiMaxwellian = false;
      isPrimitive = false;

      rho0 = 1.0;
      Bx0 = 1.0;
      By0 = 1.4142135623730951;
      Bz0 = 0.5;
      p0 = 1.0 / gamma;
      E0 = p0 / (gamma - 1) + 0.5 * (sqr(Bx0) + sqr(By0) + sqr(Bz0));
      rho1 = 8.944271909999159e-1;
      rhoux1 = -4.472135954999579e-1;
      rhouy1 = -8.432740427115680e-1;
      rhouz1 = -2.981423969999720e-1;
      E1 = 6.708136850795449e-01;
      Bx1 = 0.0;
      By1 = -4.216370213557841e-1;
      Bz1 = -1.490711984999860e-1;
   } else if (wavetype == "alfven_ipic3d") {
      isBiMaxwellian = false;
      isPrimitive = true;

      rho0 = 1.0;
      Bx0 = 0.2;
      By0 = 0.0;
      Bz0 = 0.0;
      p0 = 5.12e-4;
      rho1 = 0.0;
      ux1 = 0.0;
      uy1 = A * Bx0;
      uz1 = A * Bx0;
      p1 = 0.0;
      Bx1 = 0.0;
      By1 = A * Bx0;
      Bz1 = A * Bx0;
   } else if (wavetype == "fast_ipic3d") {
      isBiMaxwellian = false;
      isPrimitive = true;

      rho0 = 1.0;
      Bx0 = 0.0;
      By0 = 0.04;
      Bz0 = 0.0;
      p0 = 4.5e-4;
      rho1 = A * rho0;
      ux1 = A * sqrt((sqr(By0) + 2 * p0) / rho0);
      uy1 = 0.0;
      uz1 = 0.0;
      ppar1 = A * p0;
      pperp1 = A * 0.5 * (3 * gamma - 1) * p0;
      Bx1 = 0.0;
      By1 = A * By0;
      Bz1 = 0.0;
   }

   // Convert to SI units
   rho0 *= nRef * m;
   p0 *= pRef;
   E0 *= pRef;
   Bx0 *= BRef;
   By0 *= BRef;
   Bz0 *= BRef;

   rho1 *= nRef * m;
   ux1 *= uRef;
   uy1 *= uRef;
   uz1 *= uRef;
   rhoux1 *= nRef * m * uRef;
   rhouy1 *= nRef * m * uRef;
   rhouz1 *= nRef * m * uRef;
   p1 *= pRef;
   ppar1 *= pRef;
   pperp1 *= pRef;
   E1 *= pRef;

   Bx1 *= BRef;
   By1 *= BRef;
   Bz1 *= BRef;

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         std::cout << "wavetype = " << wavetype << std::endl;
         std::cout << "amplitude = " << A << std::endl;
         std::cout << "-------------------------------------------------\n";
         std::cout << "Unit conversion factors from dimensionless to SI:" << std::endl;
         std::cout << "lRef: " << lRef << " [m]" << std::endl;
         std::cout << "tRef: " << tRef << " [s]" << std::endl;
         std::cout << "vRef: " << uRef << " [m/s]" << std::endl;
         std::cout << "nRef: " << nRef << " [#]" << std::endl;
         std::cout << "TRef: " << TRef << " [K]" << std::endl;
         std::cout << "pRef: " << pRef << " [Pa]" << std::endl;
         std::cout << "BRef: " << BRef << " [T]" << std::endl;
         std::cout << "-------------------------------------------------\n";
         creal Lx = P::xmax - P::xmin;
         std::cout << "wavelength / di = " << Lx / lRef << std::endl;
         // Useful for setting the Vspace ~ 10*Vth width
         Real Vth = sqrt(gamma * p0 / rho0);
         std::cout << "Vth: " << Vth << " [m/s]" << std::endl;
         const vmesh::MeshParameters& mp = getObjectWrapper().velocityMeshes[0];
         Real dvx = (mp.meshLimits[1] - mp.meshLimits[0]) / (mp.gridLength[0] * mp.blockLength[0]);
         std::cout << "Vmax:" << mp.meshLimits[1] << std::endl;
         std::cout << "Vmin:" << mp.meshLimits[0] << std::endl;
         std::cout << "nvcells:" << mp.gridLength[0] * mp.blockLength[0] << std::endl;
         std::cout << "dvx = " << dvx << std::endl;
         std::cout << "Vth/dvx = " << Vth / dvx << std::endl;
         std::cout << "Vmax/Vth = " << mp.meshLimits[1] / Vth << std::endl;
         if (10 * dvx > Vth) {
            std::cout << "vmesh maybe too coarse!" << dvx << std::endl;
         }
         if (5 * Vth > mp.meshLimits[1]) {
            std::cout << "vmesh extent maybe too small!" << dvx << std::endl;
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
   RP::add("LinearWave.isPrimitive", "Use primitive variable inputs", 1);
   RP::add("LinearWave.isBiMaxwellian", "Use BiMaxwellian distribution", 0);
   RP::add("LinearWave.verbose", "Turn on/off detailed information", 0);
   RP::add("LinearWave.wavetype", "Type of MHD waves with preset parameters (\"fast\", \"alfven\", \"slow\", \"none\")",
           "none");
   // SI units inputs
   RP::add("LinearWave.BRef", "Reference magnetic field strength", 1e-8);
   RP::add("LinearWave.nRef", "Reference number density", 1e6);
   // Per-population parameters
   for (uint i = 0; i < getObjectWrapper().particleSpecies.size(); i++) {
      const std::string& pop = getObjectWrapper().particleSpecies[i].name;
      RP::add(pop + "_LinearWave.rho0", "Background mass density", 1.0);
      RP::add(pop + "_LinearWave.p0", "Background pressure", 0.0);
      RP::add(pop + "_LinearWave.E0", "Background energy density", 0.0);
      // Perturbation right-eigenvectors
      RP::add(pop + "_LinearWave.rho1", "Perturbation mass density", 0.0);
      RP::add(pop + "_LinearWave.ux1", "Perturbation x-velocity", 0.0);
      RP::add(pop + "_LinearWave.uy1", "Perturbation y-velocity", 0.0);
      RP::add(pop + "_LinearWave.uz1", "Perturbation z-velocity", 0.0);
      RP::add(pop + "_LinearWave.rhoux1", "Perturbation x-momentum", 0.0);
      RP::add(pop + "_LinearWave.rhouy1", "Perturbation y-momentum", 0.0);
      RP::add(pop + "_LinearWave.rhouz1", "Perturbation z-momentum", 0.0);
      RP::add(pop + "_LinearWave.ppar1", "Perturbation parallel pressure", 0.0);
      RP::add(pop + "_LinearWave.pperp1", "Perturbation perpendicular pressure", 0.0);
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
   RP::get("LinearWave.isPrimitive", isPrimitive);
   RP::get("LinearWave.isBiMaxwellian", isBiMaxwellian);
   RP::get("LinearWave.wavetype", wavetype);
   RP::get("LinearWave.BRef", BRef);
   RP::get("LinearWave.nRef", nRef);

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
      RP::get(pop + "_LinearWave.ux1", ux1);
      RP::get(pop + "_LinearWave.uy1", uy1);
      RP::get(pop + "_LinearWave.uz1", uz1);
      RP::get(pop + "_LinearWave.rhoux1", rhoux1);
      RP::get(pop + "_LinearWave.rhouy1", rhouy1);
      RP::get(pop + "_LinearWave.rhouz1", rhouz1);
      RP::get(pop + "_LinearWave.p1", p1);
      RP::get(pop + "_LinearWave.ppar1", ppar1);
      RP::get(pop + "_LinearWave.pperp1", pperp1);
      RP::get(pop + "_LinearWave.E1", E1);
   }
}

Real LinearWave::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                               creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;
   creal Lx = P::xmax - P::xmin;
   creal kwave = 2 * M_PI / Lx;

   Real n, ux, uy, uz, T;

   if (!isPrimitive) {
      creal gamma = 5.0 / 3.0;
      creal mu0 = physicalconstants::MU_0;
      creal perturb = sin(kwave * x);
      creal rho = rho0 + perturb * rho1;

      n = rho / m;
      // Assume background velocity 0
      ux = perturb * rhoux1 / rho;
      uy = perturb * rhouy1 / rho;
      uz = perturb * rhouz1 / rho;

      creal p =
          (gamma - 1) * ((E0 + E1 * perturb) - 0.5 * rho * (sqr(ux) + sqr(uy) + sqr(uz)) -
                         0.5 / mu0 * (sqr(Bx0 + Bx1 * perturb) + sqr(By0 + By1 * perturb) + sqr(Bz0 + Bz1 * perturb)));

      T = p / (n * kB);
   } else {
      if (wavetype == "alfven_ipic3d") {
         n = rho0 / m;
         ux = 0.0;
         uy = -uy1 * cos(kwave * x);
         uz = uz1 * sin(kwave * x);
         T = p0 / (n * kB);
      } else {
         creal perturb = sin(kwave * x);
         creal rho = rho0 + perturb * rho1;

         n = rho / m;
         // Assume background velocity 0
         ux = perturb * ux1;
         uy = perturb * uy1;
         uz = perturb * uz1;

         creal p = p0 + perturb * p1;
         T = p / (n * kB);
      }
   }

   creal coef = m / (2 * M_PI * kB * T);
   // Maxwellian f(v) = f(v;n,u,T)
   creal f = n * sqrt(coef) * coef * exp(-coef * M_PI * (sqr(vx - ux) + sqr(vy - uy) + sqr(vz - uz)));

   return f;
}

Real LinearWave::getBiMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy,
                                 creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;
   creal Lx = P::xmax - P::xmin;
   creal kwave = 2 * M_PI / Lx;
   creal perturb = sin(kwave * x);
   creal rho = rho0 + perturb * rho1;

   Real n, ux, uy, uz, Tpar, Tperp;

   n = rho / m;
   // Assume background velocity 0
   ux = ux1 * perturb;
   uy = uy1 * perturb;
   uz = uz1 * perturb;

   Tpar = (p0 + perturb * ppar1) / (n * kB);
   Tperp = (p0 + perturb * pperp1) / (n * kB);

   creal coef = m / (2 * M_PI * kB);
   // BiMaxwellian f(v) = f(v;n,u,Tpar, Tperp) with B0 in the y-direction
   creal f = n * sqrt(coef) * coef / (Tperp * sqrt(Tpar)) *
             exp(-coef * M_PI / Tperp * (sqr(vx - ux) + sqr(vz - uz)) - coef * M_PI / Tpar * sqr(vy - uy));

   return f;
}

Real LinearWave::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                       creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                       const uint popID) const {
   Real f = 0.0;
   if (!isBiMaxwellian) {
      f = getMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx,
                        dvy, dvz, popID);
   } else {
      f = getBiMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx,
                          dvy, dvz, popID);
   }

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
   bgField.initialize(Bx0, By0, Bz0);
   setBackgroundField(bgField, BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();
      creal Lx = P::xmax - P::xmin;
      creal kwave = 2 * M_PI / Lx;

      if (wavetype == "alfven_ipic3d") {
#pragma omp parallel for collapse(3)
         for (int i = 0; i < localSize[0]; ++i) {
            for (int j = 0; j < localSize[1]; ++j) {
               for (int k = 0; k < localSize[2]; ++k) {
                  const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
                  std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);
                  cell->at(fsgrids::bfield::PERBX) = 0.0;
                  cell->at(fsgrids::bfield::PERBY) = By1 * cos(kwave * x[0]);
                  cell->at(fsgrids::bfield::PERBZ) = -Bz1 * sin(kwave * x[0]);
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
                  creal perturb = sin(kwave * x[0]);
                  cell->at(fsgrids::bfield::PERBX) = perturb * Bx1;
                  cell->at(fsgrids::bfield::PERBY) = perturb * By1;
                  cell->at(fsgrids::bfield::PERBZ) = perturb * Bz1;
               }
            }
         }
      }
   }
}

} // namespace projects
