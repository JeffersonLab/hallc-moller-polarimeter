#include "MolPolEMFieldSetup.hh"
#include "MolPolEMFieldMessenger.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "MolPolQuad.hh"
#include "MolPolSolenoid.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4ios.hh"

#include <iostream>
using namespace std;


MolPolEMFieldSetup::MolPolEMFieldSetup()
  : fFieldManager(0),
    fChordFinder(0),
    fStepper(0),
    fIntgrDriver(0),
    fFieldMessenger(0),
    fStepperType(0),
    fMinStep(0)
{

  fMagSourceMode = 0;
  fQ1A = 0;
  fQ2A = 0;
  fQ3A = 0;
  fQ6A = 0;

  fQ1T = 0;
  fQ2T = 0;
  fQ3T = 0;
  fQ6T = 0;

  InitialseAll();

}

void MolPolEMFieldSetup::InitialseAll()
{

  fFieldMessenger = new MolPolEMFieldMessenger(this);

  fEMfield = new MolPolEMField();
  fEquation = new G4EqMagElectricField(fEMfield);
  fMinStep  = 0.01*mm ; // minimal step of 1 miron, default is 0.01 mm: Doesn't seem to make much difference here
  fStepperType = 4 ;    // ClassicalRK4 -- the default stepper

  fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  fChordFinder = 0;   //will be set in UpdateField()
  UpdateField();

  G4RotationMatrix* NOROT = new G4RotationMatrix;

  G4double ORIGINQ1 =  85.0 * cm;
  G4double ORIGINQ2 = 297.8 * cm;
  G4double ORIGINQ3 = 431.1 * cm;
  G4double ORIGINQ6 = 0.0    * cm;

  G4double BORERADIUS = 5.08 * cm;//what is this?

  fMagSourceMode = 1;

  G4double KAPPA1 = 0.;
  G4double KAPPA2 = 0.;
  G4double KAPPA3 = 0.;
  G4double SOLENOID = 0.;

  if( fMagSourceMode == 0 ){
      KAPPA1 = CalA2T(fQ1A, 1) / BORERADIUS;
      KAPPA2 = CalA2T(fQ2A, 2) / BORERADIUS;
      KAPPA3 = CalA2T(fQ3A, 3) / BORERADIUS;
  } else if( fMagSourceMode == 1){
      KAPPA1 = fQ1T * tesla / BORERADIUS;
      KAPPA2 = fQ2T * tesla / BORERADIUS;
      KAPPA3 = fQ3T * tesla / BORERADIUS;
  }

  SOLENOID = fQ6T * tesla;

  G4cout << "Received values from macro: " << G4endl;
  G4cout << "fQ1T: " << fQ1T << G4endl;
  G4cout << "fQ2T: " << fQ2T << G4endl;
  G4cout << "fQ3T: " << fQ3T << G4endl;
  G4cout << "fQ6T: " << fQ6T << G4endl;

  G4cout << __PRETTY_FUNCTION__ <<"\t at line: "<<__LINE__<<G4endl;
  G4cout << "\tfMagSourceMode: "<<fMagSourceMode<<G4endl
	 << "\tKAPPA1: "<<KAPPA1/(tesla / m)<< " tesla/m"<<G4endl
	 << "\tKAPPA2: "<<KAPPA2/(tesla / m)<< " tesla/m"<<G4endl
	 << "\tKAPPA3: "<<KAPPA3/(tesla / m)<< " tesla/m"<<G4endl
   << "\tSOLENOID: "<<SOLENOID/tesla<< " tesla"<<G4endl;

  //MolPolQuad(G4double pGradient, G4ThreeVector pOrigin, G4RotationMatrix* pMatrix, G4double pRadius)
  fMagFieldFZB1 = new MolPolQuad(KAPPA1, G4ThreeVector(0.0, 0.0, ORIGINQ1), NOROT, BORERADIUS);
  fEquationFZB1 = new G4Mag_UsualEqRhs(fMagFieldFZB1);
  fStepperFZB1  = new G4ClassicalRK4(fEquationFZB1);
  fLocalFieldManagerFZB1 = new G4FieldManager();
  fChordFinderFZB1 = 0;
  UpdateFieldFZB1();

  fMagFieldFZB2 = new MolPolQuad(KAPPA2, G4ThreeVector(0.0, 0.0, ORIGINQ2), NOROT, BORERADIUS);
  fEquationFZB2 = new G4Mag_UsualEqRhs(fMagFieldFZB2);
  fStepperFZB2  = new G4ClassicalRK4(fEquationFZB2);
  fLocalFieldManagerFZB2 = new G4FieldManager();
  fChordFinderFZB2 = 0;
  UpdateFieldFZB2();

  fMagFieldFZB3 = new MolPolQuad(KAPPA3, G4ThreeVector(0.0, 0.0, ORIGINQ3), NOROT, BORERADIUS);
  fEquationFZB3 = new G4Mag_UsualEqRhs(fMagFieldFZB3);
  fStepperFZB3  = new G4ClassicalRK4(fEquationFZB3);
  fLocalFieldManagerFZB3 = new G4FieldManager();
  fChordFinderFZB3 = 0;
  UpdateFieldFZB3();

  fMagFieldFZB6 = new MolPolSolenoid(SOLENOID, 0, G4ThreeVector(0.0, 0.0, ORIGINQ6));
  fEquationFZB6 = new G4Mag_UsualEqRhs(fMagFieldFZB6);
  fStepperFZB6  = new G4ClassicalRK4(fEquationFZB6);
  fLocalFieldManagerFZB6 = new G4FieldManager();
  fChordFinderFZB6 = 0;
  UpdateFieldFZB6();


}

/////////////////////////////////////////////////////////////////////////////////
//

MolPolEMFieldSetup::~MolPolEMFieldSetup()
{
  if(fChordFinder)    delete fChordFinder;
  if(fStepper)        delete fStepper;
  if(fEquation)       delete fEquation;
  if(fEMfield)        delete fEMfield;
  if(fFieldMessenger) delete fFieldMessenger;
}
/////////////////////////////////////////////////////////////////////////////////
//

void MolPolEMFieldSetup::UpdateConfiguration(){

  G4RotationMatrix* NOROT = new G4RotationMatrix;

  G4double ORIGINQ1 =  85.0 * cm;
  G4double ORIGINQ2 = 297.8 * cm;
  G4double ORIGINQ3 = 431.1 * cm;
  G4double ORIGINQ6 = 0.0 * cm;

  G4double BORERADIUS = 5.08 * cm;//what is this for our magnets?

  G4double KAPPA1 = 0.;
  G4double KAPPA2 = 0.;
  G4double KAPPA3 = 0.;
  G4double SOLENOID = 0.;

  if( fMagSourceMode == 0 ){
      KAPPA1 = CalA2T(fQ1A, 1) / BORERADIUS;
      KAPPA2 = CalA2T(fQ2A, 2) / BORERADIUS;
      KAPPA3 = CalA2T(fQ3A, 3) / BORERADIUS;
  } else if( fMagSourceMode == 1){
      KAPPA1 = fQ1T * tesla / BORERADIUS;
      KAPPA2 = fQ2T * tesla / BORERADIUS;
      KAPPA3 = fQ3T * tesla / BORERADIUS;
  }

  SOLENOID = fQ6T * tesla;

  G4cout << "Received values from macro: " << G4endl;
  G4cout << "fQ1T: " << fQ1T << G4endl;
  G4cout << "fQ2T: " << fQ2T << G4endl;
  G4cout << "fQ3T: " << fQ3T << G4endl;
  G4cout << "fQ6T: " << fQ6T << G4endl;

  G4cout << __PRETTY_FUNCTION__ <<"\t at line: "<<__LINE__<<G4endl;
  G4cout << "\tfMagSourceMode: "<<fMagSourceMode<<G4endl
	 << "\tKAPPA1: "<<KAPPA1/(tesla / m)<< " tesla/m"<<G4endl
	 << "\tKAPPA2: "<<KAPPA2/(tesla / m)<< " tesla/m"<<G4endl
	 << "\tKAPPA3: "<<KAPPA3/(tesla / m)<< " tesla/m"<<G4endl
   << "\tSOLENOID: "<<SOLENOID/tesla<< " tesla"<<G4endl;

  fMagFieldFZB1->UpdateQuad(KAPPA1, G4ThreeVector(0.0, 0.0, ORIGINQ1), NOROT, BORERADIUS);
  fMagFieldFZB2->UpdateQuad(KAPPA2, G4ThreeVector(0.0, 0.0, ORIGINQ2), NOROT, BORERADIUS);
  fMagFieldFZB3->UpdateQuad(KAPPA3, G4ThreeVector(0.0, 0.0, ORIGINQ3), NOROT, BORERADIUS);
  fMagFieldFZB6->UpdateSolenoid(SOLENOID, 0, G4ThreeVector(0.0, 0.0, ORIGINQ6));
}

/////////////////////////////////////////////////////////////////////////////
//
// Register this field to 'global' Field Manager and
// Create Stepper and Chord Finder with predefined type, minstep (resp.)
//

void MolPolEMFieldSetup::UpdateField()
{
  fStepper = new G4ClassicalRK4( fEquation, 8 );

  fFieldManager->SetDetectorField(fEMfield);

  if(fChordFinder) delete fChordFinder;
  fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  fFieldManager->SetChordFinder( fChordFinder );
}


/////////////////////////////////////////////////////////////////////////////
void MolPolEMFieldSetup::UpdateFieldFZB1()
{

  fLocalFieldManagerFZB1->SetDetectorField(fMagFieldFZB1);

  if(fChordFinderFZB1) delete fChordFinderFZB1;
  fIntgrDriverFZB1 = new G4MagInt_Driver(fMinStep,fStepperFZB1,fStepperFZB1->GetNumberOfVariables());
  fChordFinderFZB1 = new G4ChordFinder((G4MagneticField*) fMagFieldFZB1, fMinStep, fStepperFZB1);
  fLocalFieldManagerFZB1->SetChordFinder( fChordFinderFZB1 );

}

/////////////////////////////////////////////////////////////////////////////
void MolPolEMFieldSetup::UpdateFieldFZB2()
{

  fLocalFieldManagerFZB2->SetDetectorField(fMagFieldFZB2);

  if(fChordFinderFZB2) delete fChordFinderFZB2;
  fIntgrDriverFZB2 = new G4MagInt_Driver(fMinStep,fStepperFZB2,fStepperFZB2->GetNumberOfVariables());
  fChordFinderFZB2 = new G4ChordFinder((G4MagneticField*) fMagFieldFZB2, fMinStep, fStepperFZB2);
  fLocalFieldManagerFZB2->SetChordFinder( fChordFinderFZB2 );

}

/////////////////////////////////////////////////////////////////////////////
void MolPolEMFieldSetup::UpdateFieldFZB3()
{

  fLocalFieldManagerFZB3->SetDetectorField(fMagFieldFZB3);

  if(fChordFinderFZB3) delete fChordFinderFZB3;
  fIntgrDriverFZB3 = new G4MagInt_Driver(fMinStep,fStepperFZB3,fStepperFZB3->GetNumberOfVariables());
  fChordFinderFZB3 = new G4ChordFinder((G4MagneticField*) fMagFieldFZB3, fMinStep, fStepperFZB3);
  fLocalFieldManagerFZB3->SetChordFinder( fChordFinderFZB3 );

}

/////////////////////////////////////////////////////////////////////////////

void MolPolEMFieldSetup::UpdateFieldFZB6()
{

  fLocalFieldManagerFZB6->SetDetectorField(fMagFieldFZB6);

  if(fChordFinderFZB6) delete fChordFinderFZB6;
  fIntgrDriverFZB6 = new G4MagInt_Driver(fMinStep,fStepperFZB6,fStepperFZB6->GetNumberOfVariables());
  fChordFinderFZB6 = new G4ChordFinder((G4MagneticField*) fMagFieldFZB6, fMinStep, fStepperFZB6);
  fLocalFieldManagerFZB6->SetChordFinder( fChordFinderFZB6 );

}


/////////////////////////////////////////////////////////////////////////////
//
// Set stepper according to the stepper type
//


void MolPolEMFieldSetup::SetStepper()
{
  G4int nvar = 8;

  if(fStepper) delete fStepper;

  switch ( fStepperType )
    {
    case 0:
      fStepper = new G4ExplicitEuler( fEquation, nvar );
      break;
    case 1:
      fStepper = new G4ImplicitEuler( fEquation, nvar );
      break;
    case 2:
      fStepper = new G4SimpleRunge( fEquation, nvar );
      break;
    case 3:
      fStepper = new G4SimpleHeum( fEquation, nvar );
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation, nvar );
      break;
    case 5:
      fStepper = new G4CashKarpRKF45( fEquation, nvar );
      break;
    default: fStepper = 0;
    }
}


///////////////////////////////////////////////////////////////////////////////
// Current to Field Calculation (paprameters from quads_subs.f)
// Return field at the pole tip in Tesla (Quadrupole)

G4double MolPolEMFieldSetup::CalA2T(G4double current, G4int magnet)
{

 
  G4double fld = 0;
  G4double a = 0.;
  G4double b = 0.;
  G4double c = 0.;
  G4double d = 0.;
  G4double f = 0.;

  if(magnet == 1)
    {
      //Los Alamos quad
      a = -0.17085*pow(10.0,-6.);
      b = 0.166073*pow(10.0,-2.);
      c = 0.12525*pow(10.0,-1.); 

      fld = a*current*current + b*current + c;
    }
  else if(magnet == 2 || magnet ==3)
    {
      //Argonne quads
      f = -0.86527*pow(10.0,-13.);
      a = 0.23769*pow(10.0,-10.);
      b = 0.79144*pow(10.0,-7.);
      c = 0.10542*pow(10.0,-2.);
      d = 0.70794*pow(10.0,-2.);

      fld =  f*current*current*current*current + a*current*current*current + b*current*current + c*current + d;
    }
 
  else
    {
      //wrong magnet setup
      fld = 0.0;
    }

  return fld; //units in tesla

}
