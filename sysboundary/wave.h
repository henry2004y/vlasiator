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

#ifndef WAVE_H
#define WAVE_H

#include <vector>
#include "../definitions.h"
#include "../readparameters.h"
#include "../spatial_cell.hpp"
#include "sysboundarycondition.h"

namespace SBC {

/*!\brief Class for boundary conditions with wave settings.
 */
class Wave : public SysBoundaryCondition {
public:
   Wave();
   ~Wave() override;

   static void addParameters();
   void getParameters() override;

   void initBoundary(creal t, Project &project) override;
   void assignBoundary(dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid) override;
   void applyInitialState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                          FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid,
                          Project &project) override;
   void updateState(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                    FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &perBGrid, creal t) override;
   Real
   fieldSolverBoundaryCondMagneticField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
                                        FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j,
                                        cint k, creal dt, cuint component) override;
   void fieldSolverBoundaryCondMagneticFieldProjection(
       FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH> &bGrid,
       FsGrid<fsgrids::technical, FS_STENCIL_WIDTH> &technicalGrid, cint i, cint j, cint k) override;
   void
   fieldSolverBoundaryCondElectricField(FsGrid<std::array<Real, fsgrids::efield::N_EFIELD>, FS_STENCIL_WIDTH> &EGrid,
                                        cint i, cint j, cint k, cuint component) override;
   void fieldSolverBoundaryCondHallElectricField(
       FsGrid<std::array<Real, fsgrids::ehall::N_EHALL>, FS_STENCIL_WIDTH> &EHallGrid, cint i, cint j, cint k,
       cuint component) override;
   void fieldSolverBoundaryCondGradPeElectricField(
       FsGrid<std::array<Real, fsgrids::egradpe::N_EGRADPE>, FS_STENCIL_WIDTH> &EGradPeGrid, cint i, cint j, cint k,
       cuint component) override;
   void fieldSolverBoundaryCondDerivatives(
       FsGrid<std::array<Real, fsgrids::dperb::N_DPERB>, FS_STENCIL_WIDTH> &dPerBGrid,
       FsGrid<std::array<Real, fsgrids::dmoments::N_DMOMENTS>, FS_STENCIL_WIDTH> &dMomentsGrid, cint i, cint j, cint k,
       cuint RKCase, cuint component) override;
   void fieldSolverBoundaryCondBVOLDerivatives(
       FsGrid<std::array<Real, fsgrids::volfields::N_VOL>, FS_STENCIL_WIDTH> &volGrid, cint i, cint j, cint k,
       cuint component) override;
   void vlasovBoundaryCondition(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry> &mpiGrid,
                                const CellID &cellID, const uint popID, creal t, const bool doCalcMomentsV) override;

   void getFaces(std::array<bool, 6>& faces) override;

   std::string getName() const override;
   uint getIndex() const override;

protected:
   /*! Faces processed by this boundary condition.*/
   std::array<bool, 6> facesToProcess = {false};
   /*! List of faces on which wave boundary conditions are to be applied ([xyz][+-]). */
   std::vector<std::string> faceList;

   static constexpr Real PI = 3.1415926;
   /*! Wave vector.*/
   std::array<Real, 3> waveVector;
   /*! Wave frequency.*/
   Real freq, n0, T0, n1, T1;
   std::array<Real, 3> B0;
   std::array<Real, 3> V0;
   std::array<Real, 3> B1;
   std::array<Real, 3> V1;

   // changing Bx, By, Bz?

   void setCell(SpatialCell* cell, const uint popID, creal t, const bool doCalcMomentsV);
   void setB(const dccrg::Dccrg<SpatialCell, dccrg::Cartesian_Geometry>& mpiGrid,
             FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid, creal t);
   Real maxwellian(const uint popID, creal rho, creal T, creal Vx, creal Vy, creal Vz);
   std::vector<vmesh::GlobalID> findBlocksToInit(const uint popID, SpatialCell* cell, creal rho, creal T, creal Vx,
                                                 creal Vy, creal Vz);
};
} // namespace SBC

#endif
