/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

// natural convection of air in a square cavity in 2D


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

// #define SMAGORINSKY

#ifdef SMAGORINSKY
typedef D2Q9<FORCE,TAU_EFF> NSDESCRIPTOR;
typedef D2Q5<VELOCITY,TAU_EFF> TDESCRIPTOR;
#else
typedef D2Q9<FORCE> NSDESCRIPTOR;
typedef D2Q5<VELOCITY> TDESCRIPTOR;
#endif

// Parameters for the simulation setup
T Ra = 1e6;  // Rayleigh-Zahl
const T Pr = 0.5; // Prandtl-Zahl
const T epsilon = 5.e-3;  // precision of the convergence (residuum)

#ifdef SMAGORINSKY
const int statisticsIntervall = 10; // take the turbulent statistics every 10 time steps after convergence
const int statisticsEnsembles = 200; // take 20 ensembles for the turbulent statistics
#endif

const T Thot = 280.15;
const T Tin = 300;



const int N = 400;       // resolution of the model
const T Re = 20.;       // Reynolds number
const T maxPhysT = 3000;  // max. simulation time in s, SI unit
const T L = 0.1 / N;      // latticeL
const T lengthX = 2.2;
const T lx = 2.2;
const T lengthY = 2.2;
const T centerCylinderX = 1.1;
const T centerCylinderY = 1.1;
const T radiusCylinder = 0.1;


/// Compute the nusselt number at the left wall
T computeNusselt(SuperGeometry2D<T>& superGeometry,
                 SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice2D<T, TDESCRIPTOR>& ADlattice)
{
  int voxel = 0, material = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0;
  T q = 0;

  for (int iC = 0; iC < NSlattice.getLoadBalancer().size(); iC++) {
    int ny = NSlattice.getBlockLattice(iC).getNy();
    int iX = 0;
    for (int iY = 0; iY < ny; ++iY) {
      material = superGeometry.getBlockGeometry(iC).getMaterial(iX,iY);

      T_x = ADlattice.getBlockLattice(iC).get(iX,iY).computeRho();
      T_xplus1 = ADlattice.getBlockLattice(iC).get(iX+1,iY).computeRho();
      T_xplus2 = ADlattice.getBlockLattice(iC).get(iX+2,iY).computeRho();

      if ( material == 2 ) {
        q += (3.0*T_x - 4.0*T_xplus1 + 1.0*T_xplus2)/2.0*N;
        voxel++;
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(q, MPI_SUM);
  singleton::mpi().reduceAndBcast(voxel, MPI_SUM);
#endif

  return q / (T)voxel;
}

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry2D<T>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{

 OstreamManager clout(std::cout, "prepareGeometry");
 clout << "Prepare Geometry ..." << std::endl;

 Vector<T, 2> extend(lengthX, lengthY);
 Vector<T, 2> center(centerCylinderX, centerCylinderY);
 Vector<T, 2> origin;
 IndicatorCircle2D<T> circle(center, radiusCylinder);

 superGeometry.rename(0, 2);

 superGeometry.rename(2, 1, 1, 1);

 // Set material number for inflow
 extend[0] = lengthX/N;
 origin[0] = -lengthX / N;
 IndicatorCuboid2D<T> inflow(extend, origin);
 superGeometry.rename(2, 3, 1, inflow);
 // Set material number for outflow
 origin[0] = lengthX - L;
 IndicatorCuboid2D<T> outflow(extend, origin);
 superGeometry.rename(2, 4, 1, outflow);
 // Set material number for cylinder
 superGeometry.rename(1, 5, circle);

 // Removes all not needed boundary voxels outside the surface
 superGeometry.clean();
 superGeometry.checkForErrors();

 superGeometry.print();

 clout << "Prepare Geometry ... OK" << std::endl;

}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                     ForcedBGKdynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  T omega  =  converter.getLatticeRelaxationFrequency();
  T Tomega  =  converter.getLatticeThermalRelaxationFrequency();

  ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>());

  ADlattice.defineDynamics(superGeometry, 5, &instances::getNoDynamics<T, TDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 5, &instances::getNoDynamics<T, NSDESCRIPTOR>());

  ADlattice.defineDynamics(superGeometry.getMaterialIndicator({1,2,3,4}), &advectionDiffusionBulkDynamics);
 // ADlattice.defineDynamics(superGeometry, 5, &instances::getBounceBack<T, TDESCRIPTOR>());

  NSlattice.defineDynamics(superGeometry.getMaterialIndicator({1,2,3,4}), &bulkDynamics);
 // NSlattice.defineDynamics(superGeometry, 5, &instances::getBounceBack<T, NSDESCRIPTOR>());

  Vector<T, 2> center(centerCylinderX, centerCylinderY);
  IndicatorCircle2D<T> circle(center, radiusCylinder);

  setBouzidiZeroVelocityBoundary<T, NSDESCRIPTOR>(NSlattice, superGeometry, 5, circle);
  setAdvectionDiffusionTemperatureBoundary<T, TDESCRIPTOR>(ADlattice, Tomega, superGeometry.getMaterialIndicator({3,5}));
  setLocalPressureBoundary<T, NSDESCRIPTOR>(NSlattice, omega, superGeometry, 4);
  setLocalVelocityBoundary<T, NSDESCRIPTOR>(NSlattice, omega, superGeometry, 3);

 // setLocalVelocityBoundary<T,NSDESCRIPTOR>(NSlattice, omega, superGeometry.getMaterialIndicator({5}));

  /// define initial conditions

  AnalyticalConst2D<T,T> T_in(converter.getLatticeTemperature(Tin));
  AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));


  AnalyticalConst2D<T, T> rho(10.);
  AnalyticalConst2D<T, T> u0(0.05, 0.05);

  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry.getMaterialIndicator({1,3,4}), rho, u0);
  NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1,3,4}), rho, u0);

  AnalyticalConst2D<T, T> rho1(1.0);
  AnalyticalConst2D<T, T> u1(0.0, 0.0);

  /// for each material set Rho, U and the Equilibrium
  NSlattice.defineRhoU(superGeometry.getMaterialIndicator({5}), rho1, u1);
  NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({5}), rho1, u1);


  ADlattice.defineRho(superGeometry, 3, T_in);
  ADlattice.iniEquilibrium(superGeometry, 1, T_in, u0);
  ADlattice.defineRho(superGeometry, 4, T_in);
  ADlattice.iniEquilibrium(superGeometry, 4, T_in, u0);
  ADlattice.defineRho(superGeometry, 5, T_hot);
  ADlattice.iniEquilibrium(superGeometry, 5, T_hot, u1);

#ifdef SMAGORINSKY
  AnalyticalConst2D<T,T> tauNS(1./omega);
  AnalyticalConst2D<T,T> tauAD(1./Tomega);

  NSlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3}), tauNS );
  ADlattice.defineField<descriptors::TAU_EFF>( superGeometry.getMaterialIndicator({1, 2, 3}), tauAD );
#endif

  /// Make the lattice ready for simulation
  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                        SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                        SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                        int iT, SuperGeometry2D<T>& superGeometry)
{

  // nothing to do here

}

void getResults( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                 SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                 SuperLattice2D<T, TDESCRIPTOR>& ADlattice, int iT,
                 SuperGeometry2D<T>& superGeometry,
                 Timer<T>& timer,
                 bool converged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter2D<T> vtkWriter("thermalNaturalConvection2D");
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
  SuperLatticePhysTemperature2D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );

  AnalyticalFfromSuperF2D<T> interpolation(velocity, true);

  const int statIter = 10.;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT % statIter == 0 || converged) {

    timer.update(iT);
    timer.printStep();

    /// NSLattice statistics console output
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADLattice statistics console output
    ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));

    vtkWriter.write(iT);

    BlockReduction2D2D<T> planeReduction(temperature, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tin-.1, Thot+.1, iT, "temperature");

    SuperEuklidNorm2D<T, NSDESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction2(normVel, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );

  }

  if ( converged ) {

    T nusselt = computeNusselt(superGeometry, NSlattice, ADlattice);

    /// Initialize vectors for data output
    T xVelocity[2] = { T() };
    T outputVelX[2] = { T() };
    T yVelocity[2] = { T() };
    T outputVelY[2] = { T() };
    const int outputSize = 512;
    Vector<T, outputSize> velX;
    Vector<T, outputSize> posX;
    Vector<T, outputSize> velY;
    Vector<T, outputSize> posY;

    /// loop for the resolution of the cavity at x = lx/2 in yDirection and vice versa
    for (int n = 0; n < outputSize; ++n) {
      T yPosition[2] = { lx / 2, lx * n / (T) outputSize };
      T xPosition[2] = { lx * n / (T) outputSize, lx / 2 };

      /// Interpolate xVelocity at x = lx/2 for each yPosition
      interpolation(xVelocity, yPosition);
      interpolation(yVelocity, xPosition);
      /// Store the interpolated values to compare them among each other in order to detect the maximum
      velX[n] = xVelocity[0];
      posY[n] = yPosition[1];
      velY[n] = yVelocity[1];
      posX[n] = xPosition[0];

      /// Initialize output with the corresponding velocities and positions at the origin
      if (n == 0) {
        outputVelX[0] = velX[0];
        outputVelX[1] = posY[0];
        outputVelY[0] = velY[0];
        outputVelY[1] = posX[0];
      }
      /// look for the maximum velocity in xDirection and the corresponding position in yDirection
      if (n > 0 && velX[n] > outputVelX[0]) {
        outputVelX[0] = velX[n];
        outputVelX[1] = posY[n];
      }
      /// look for the maximum velocity in yDirection and the corresponding position in xDirection
      if (n > 0 && velY[n] > outputVelY[0]) {
        outputVelY[0] = velY[n];
        outputVelY[1] = posX[n];
      }
    }

    // compare to De Vahl Davis' benchmark solutions
    clout << "Comparison against De Vahl Davis (1983):" << endl;
  
    }
  }


int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

#ifndef SMAGORINSKY
  T tau = 0.6;
#endif
  T charU = 0.001;
  //Pr = 0.05263;

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> converter(
    (T) lx / N,
#ifdef SMAGORINSKY
    (T) 2.*0.056/charU*lx/N,
#else
    (T) (tau - 0.5) / descriptors::invCs2<T,NSDESCRIPTOR>() * pow((lx/N),2) / 15.126e-6,
#endif
    (T) lx,
    (T) charU,
    (T) 15.126e-6,
    (T) 1.0,
    (T) 25.684e-3,
    (T) Pr * 25.684e-3 / 15.126e-6 / 1.0,
    (T) 0.00341,
    (T) Tin,
    (T) Thot
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx + 2*converter.getPhysLength(1);
  extend[1] = lx + converter.getPhysLength(1);
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidGeometry
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize());

  cuboidGeometry.setPeriodicity(false, true);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice2D<T, TDESCRIPTOR> ADlattice(superGeometry);
  SuperLattice2D<T, NSDESCRIPTOR> NSlattice(superGeometry);



#ifdef SMAGORINSKY
  ExternalTauEffLESForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>(), 0.1);

  ExternalTauEffLESBGKdynamics<T, TDESCRIPTOR> TbulkDynamics(
    converter.getLatticeThermalRelaxationFrequency(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>(), 0.1);
#else
  ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>());

  AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> TbulkDynamics(
    converter.getLatticeThermalRelaxationFrequency(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>());
#endif
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<T> dir{0.0, 1.0};

  T boussinesqForcePrefactor = 0; // 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                              // converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();


  NavierStokesAdvectionDiffusionCouplingGenerator2D<T,NSDESCRIPTOR>
  coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(lx),
           boussinesqForcePrefactor, converter.getLatticeTemperature(Tin), 1., dir);

  NSlattice.addLatticeCoupling(superGeometry, 1, coupling, ADlattice);

  //prepareLattice and setBoundaryCondition
  prepareLattice(converter,
                 NSlattice, ADlattice,
                 NSbulkDynamics, TbulkDynamics,
                 superGeometry );



  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(6,epsilon);
  bool converged = false;
  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

    if (converge.hasConverged() && !converged) {
      converged = true;
      clout << "Simulation converged." << endl;
      clout << "Time " << iT << "." << std::endl;

      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    }

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    NSlattice.executeCoupling();
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();

    /// === 7th Step: Computation and Output of the Results ===
    if ( !converged ) {
      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    }
    if (!converged && iT % 1000 == 0) {
      converge.takeValue(computeNusselt(superGeometry, NSlattice, ADlattice),true);
    }
  }



  timer.stop();
  timer.printSummary();

}
