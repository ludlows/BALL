// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: molecularDynamics.sip,v 1.6 2005/01/24 16:59:23 amoll Exp $
//

class MolecularDynamics
{
%TypeHeaderCode
  #include <BALL/MOLMEC/MDSIMULATION/molecularDynamics.h>
  using namespace BALL;
%End                                                                                                                                                            
	public:
  //struct Option
  //{
  //  static const char *MAXIMAL_NUMBER_OF_ITERATIONS;
  //  static const char *MAXIMAL_SIMULATION_TIME; 
  //  static const char *NUMBER_OF_ITERATION;
  //  static const char *ENERGY_OUTPUT_FREQUENCY;
  //  static const char *SNAPSHOT_FREQUENCY;
  //  static const char *TIME_STEP;
  //  static const char *REFERENCE_TEMPERATURE;
  //  static const char *BATH_RELAXATION_TIME;
 //  static const char *CURRENT_TIME;
  // };
  // struct Default
  // {
  //  static const int MAXIMAL_NUMBER_OF_ITERATIONS;
  //  static const double  MAXIMAL_SIMULATION_TIME; 
  //  static const int NUMBER_OF_ITERATION;
  //  static const int ENERGY_OUTPUT_FREQUENCY;
  //  static const int SNAPSHOT_FREQUENCY;
  //  static const double TIME_STEP;
  //  static const double REFERENCE_TEMPERATURE;
  //  static const double BATH_RELAXATION_TIME;
  // static const double CURRENT_TIME;
  // };
  // BALL_CREATE (MolecularDynamics)
  MolecularDynamics ();
  MolecularDynamics (ForceField& /forcefield/);
  MolecularDynamics (const MolecularDynamics& /rhs/);
  ~MolecularDynamics ();
  // const MolecularDynamics& operator = (const MolecularDynamics& /rhs/);
	void set(const MolecularDynamics& /rhs/);
  bool isValid () const;
  virtual bool setup (ForceField & /myforcefield/, SnapShotManager* /snapshot_man/);
  virtual bool setup (ForceField & /myforcefield/, SnapShotManager* /snapshot_man/, const Options & /myoptions/);
  virtual bool specificSetup();
  void setNumberOfIteration (Size /number/);
  void setMaximalNumberOfIterations (Size /number/);
  void setMaximalSimulationTime(double /time/); 
  virtual void setTimeStep(double /step/);
  void setReferenceTemperature (double /temperature/);
  void setCurrentTime (double /time/);
  void setEnergyOutputFrequency (Size /number/);
  void setSnapShotFrequency (Size /number/);
  Size getEnergyOutputFrequency () const;
  Size getNumberOfIterations () const;
  Size getMaximalNumberOfIterations () const;
  double getMaximalSimulationTime () const;
  double getTimeStep () const; 
  Size getSnapShotFrequency () const;
  double getTemperature () const;
  double getTime () const;
  double getTotalEnergy () const;
  double getPotentialEnergy () const;
  double getKineticEnergy () const;
  ForceField *getForceField () const;
  virtual bool simulate (bool /restart/ = false);
  virtual bool simulateIterations (Size /number/, bool /restart/ = false);
  virtual bool simulateTime (double /simulation_time/, bool /restart/ = false);
  Options options;
};  
  
