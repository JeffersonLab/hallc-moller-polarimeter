#include "MolPolDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "MolPolDetector.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Para.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpticalSurface.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "MolPolEMFieldSetup.hh"
#include "G4NistManager.hh"

G4int nel;

void MolPolDetectorConstruction::DetModeSet(G4int detMode = 1) {
}

void MolPolDetectorConstruction::StandModeSet(G4int standMode = 0) {
}


MolPolDetectorConstruction::MolPolDetectorConstruction():
  fCheckOverlaps(false)
{
  mEMFieldSetup = 0;
}

MolPolDetectorConstruction::~MolPolDetectorConstruction(){
  if(mEMFieldSetup) delete mEMFieldSetup;
}

G4VPhysicalVolume* MolPolDetectorConstruction::Construct() {

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //// Define Visual Attributes
  G4double alphaVacuum = 0.15;
  G4double alphaMatStd = 0.50;
  G4double alphaTarget = 0.85;
  G4VisAttributes* IronVisAtt = new G4VisAttributes( G4Colour( 10./255., 10./255.,10./255.,alphaTarget) );
  //G4VisAttributes* LeadVisAtt = new G4VisAttributes( G4Colour(149./255.,149./255.,100./255.,alphaMatStd) );
  G4VisAttributes* SteelVisAtt= new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaMatStd) );
  G4VisAttributes* AlumVisAtt = new G4VisAttributes( G4Colour(  0./255.,237./255.,  0./255.,alphaMatStd) );
  G4VisAttributes* VacVisAtt  = new G4VisAttributes( G4Colour(255./255.,255./255.,255./255.,alphaVacuum) );
  G4VisAttributes* CuVisAtt   = new G4VisAttributes( G4Colour(178./255.,102./255., 26./255.,alphaMatStd) );
  G4VisAttributes* ScintVisAtt= new G4VisAttributes( G4Colour(  0./255.,100./255.,100./255.,alphaMatStd) );
  //G4VisAttributes* DipVisAtt  = new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaVacuum) );


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //// Materials Setup
  G4double a, z, density, pressure, temperature;
  G4int nelements, natoms;

  G4Element* N  = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  //G4Element* H  = new G4Element("Hydrogen", "H", z=1 , a=1.01 *g/mole);
  //G4Element* C  = new G4Element("Carbon",   "C", z=6, a=12.01*g/mole);
  G4Element* As  = new G4Element("Arsen","As", z=33, a=74.92*g/mole);

  G4Element* Al = new G4Element("Aluminum", "Al", z=13, a=26.98*g/mole);
  G4Element* Fe = new G4Element("Iron"   , "Fe", z=26, a=55.845*g/mole);
  G4Element* Si = new G4Element("Silicon", "Si", z=14, a=28.09 *g/mole);
  //G4Element* Pb = new G4Element("Lead", "Pb", z=82., a=207.19*g/mole);

  density = 0.787 * g/cm3;
  a = 55.85 * g /mole;
  G4Material* iron = new G4Material("iron", z=26, a, density);

  density = 7.65 *g/cm3;
  G4Material* siliconsteel = new G4Material("SiliconSteel", density, nelements=2);
  siliconsteel->AddElement(Fe, natoms=11);
  siliconsteel->AddElement(Si, natoms=1);
  density = 2.70 *g/cm3;
  G4Material* aluminum = new G4Material("aluminum", density, nelements=1);
  aluminum->AddElement(Al, natoms=1);

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  //G4Material* lead = new G4Material("lead", z=82, a, density);

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  G4Material* Cu = new G4Material("Copper" , z=29., a, density);

  density = 1.032*g/cm3;
  a = 12.01*g/mole;
  G4Material* scint = new G4Material("scint", z=6., a, density);

  density = 1.e-6/760.0 * 1.29*mg/cm3; //0.001 of air density
  pressure = 1.e-6/760.0 *atmosphere;
  temperature = 293.15 *kelvin;  //room temperature
  a = 28.97 *g/mole;

  G4Material* Vacuum = new G4Material("Vacuum",z=1,a,density,kStateGas,temperature,pressure);
  G4Material* Air    = new G4Material("Air",    density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 79.0*perCent);
  Air->AddElement(O, 21.0*perCent);

  G4NistManager* matman = G4NistManager::Instance();

  G4Material* LeadOxide   = matman->FindOrBuildMaterial("G4_LEAD_OXIDE"); 
  G4Material* SilicOxide  = matman->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); 
  G4Material* PotasOxide  = matman->FindOrBuildMaterial("G4_POTASSIUM_OXIDE"); 
  G4Material* SodMonOxide = matman->FindOrBuildMaterial("G4_SODIUM_MONOXIDE");  
  G4Material* As2O3 = new G4Material("As2O3", density= 3.738*g/cm3, nel=2);
  As2O3->AddElement( As, 2 );
  As2O3->AddElement( O,  3 );

  //Lead glass:
  G4Material* LgTF1 = new G4Material("LgTF1", density= 3.86*g/cm3 , nel=5 ); 
  LgTF1->AddMaterial( LeadOxide   , 0.5120 );
  LgTF1->AddMaterial( SilicOxide  , 0.4130 );
  LgTF1->AddMaterial( PotasOxide  , 0.0422 );
  LgTF1->AddMaterial( SodMonOxide , 0.0278 );
  LgTF1->AddMaterial( As2O3       , 0.0050 );

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Build world
  G4double world_x = 14*m;  G4double world_y = 14*m;  G4double world_z = 14*m;
  G4Box* world_box = new G4Box("World",world_x,world_y,world_z);
  G4LogicalVolume* world_log = new G4LogicalVolume(world_box,Vacuum,"World",0,0,0);
  world_log->SetVisAttributes(G4VisAttributes::Invisible);
  G4VPhysicalVolume* world_phys = new G4PVPlacement(0,G4ThreeVector(),world_log,"World",0,false,fCheckOverlaps);

  /////////////////////////////////////////////////////////////
  //Defining Variables
  G4double theta = 3.25 * deg;

  G4RotationMatrix* RotPos = new G4RotationMatrix();
  RotPos->rotateY(theta);

  G4RotationMatrix* RotNeg = new G4RotationMatrix();
  RotNeg->rotateY(-theta);

  G4double ZeroR = 0.0 * cm;

  G4double pQ2HL = 61.595 * cm;

  G4double DistTrgtQ3 = 431.099 * cm; //Distance between the target and Q3. Useful for writing distances in terms of angles



  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Field Setup
  mEMFieldSetup = new MolPolEMFieldSetup();  //setup the field,
  G4FieldManager* Q1FieldManager = mEMFieldSetup->GetFieldManagerFZB1();
  G4FieldManager* Q2FieldManager = mEMFieldSetup->GetFieldManagerFZB2();
  G4FieldManager* Q3FieldManager = mEMFieldSetup->GetFieldManagerFZB3();
  G4FieldManager* Q6FieldManager = mEMFieldSetup->GetFieldManagerFZB6();
  G4bool allLocal = true;

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Target
  G4double pMTATRin   = 0.0 * cm; G4double pMTATRout  = 1.5 * cm;   G4double pMTATHLZ = 0.0062 * cm;
  G4double pMTATPos_X = 0.0 * cm; G4double pMTATPos_Y = 0.0 * cm; G4double pMTATPos_Z = 0.0 * cm;
  G4VSolid* MTATSolid = new G4Tubs( "MTATTube", pMTATRin, pMTATRout, pMTATHLZ, 0.0, 360.0 * deg );
  G4LogicalVolume* MTATLogical = new G4LogicalVolume(MTATSolid, iron, "Target", 0, 0, 0);
  MTATLogical->SetVisAttributes(IronVisAtt);
  new G4PVPlacement(0, G4ThreeVector( pMTATPos_X,pMTATPos_Y,pMTATPos_Z), MTATLogical, "Target", world_log, 0, 0, fCheckOverlaps);


  ///////////////////////////////////////////////////////////////
  // Target Chamber
  G4double ChamberRI = 41.592 * cm;
  G4double ChamberRO = 42.545 * cm;
  G4double ChamberH = 31.120 * cm;
  
  G4RotationMatrix * ChamRot = new G4RotationMatrix();
  ChamRot -> rotateX(90.0*deg);
  
  G4double ChamberPos_X = 0.0 * cm;
  G4double ChamberPos_Y = 0.0 * cm;
  G4double ChamberPos_Z = 0.0 * cm;

  G4VSolid*Cham1 = new G4Tubs( "ChamTub", ZeroR, ChamberRO,ChamberH, 0.0, 360.0*deg);
  G4VSolid*Cham2 = new G4Tubs( "InnerTub", ZeroR, ChamberRI, ChamberH, 0.0, 360.0*deg);
  G4SubtractionSolid*Cham_Tub = new G4SubtractionSolid( "Cham_Tub", Cham1, Cham2, 0, G4ThreeVector(ChamberPos_X,ChamberPos_Y,ChamberPos_Z) );
  G4LogicalVolume*ChamLogical = new G4LogicalVolume(Cham_Tub, aluminum, "ChamLogical" , 0, 0, 0 );
  G4LogicalVolume*Cham2Logical = new G4LogicalVolume( Cham2, Vacuum, "Cham2Logical", 0, 0, 0);
  ChamLogical -> SetVisAttributes(AlumVisAtt);
  Cham2Logical -> SetVisAttributes(VacVisAtt);
  new G4PVPlacement(ChamRot, G4ThreeVector(ChamberPos_X, ChamberPos_Y, ChamberPos_Z),ChamLogical, "TrgtCham", world_log, 0, 0, fCheckOverlaps);
  // new G4PVPlacement(ChamRot, G4ThreeVector(ChamberPos_X, ChamberPos_Y, ChamberPos_Z), Cham2Logical, "VacTrgtCham",world_log, 0, 0, fCheckOverlaps);

  G4double ChamExitRO = 10.477 * cm;
  G4double ChamExitRI = 10.160 * cm;
  G4double ChamExitH = 5.0 * cm;
  
  G4double ChamExit_X = 0.0 * cm;
  G4double ChamExit_Y = 0.0 * cm; 
  G4double ChamExit_Z = (ChamberRO/cm + ChamExitH/cm) * cm;

  G4VSolid*ChamExit = new G4Tubs( "ChamTubExit", ChamExitRI, ChamExitRO,ChamExitH, 0.0, 360.0*deg);
  //G4VSolid*ChamExit2 = new G4Tubs( "ChamTubExit2", ZeroR, ChamExitRI, ChamExitH, 0.0, 360.0*deg);
  //G4SubtractionSolid*Cham_Exit = new G4SubtractionSolid("Cham_Exit",ChamExit, ChamExit2,0,G4ThreeVector(ChamExit_X,ChamExit_Y,ChamExit_Z) );
  G4LogicalVolume*ChamExitLogical = new G4LogicalVolume(ChamExit, aluminum, "ChamExitLogical", 0,0,0);
  //G4LogicalVolume*ChamExit2Logical = new G4LogicalVolume(ChamExit2, Vacuum, "ChamExit2Logical",0,0,0);
  ChamExitLogical -> SetVisAttributes(AlumVisAtt);
  //ChamExit2Logical -> SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0, G4ThreeVector(ChamExit_X, ChamExit_Y, ChamExit_Z),ChamExitLogical,"TrgtChamExit",world_log, 0, 0, fCheckOverlaps);
  // new G4PVPlacement(0, G4ThreeVector(ChamExit_X, ChamExit_Y, ChamExit_Z), ChamExit2Logical, "TrgtChamExitInner", world_log,0,0,fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // HELMHOLTZ COIL Magnetic Volume
  // TODO: don't know pQ6Rout, pQ6HL (look this up)
  G4double pQ6Rin  =  0    * cm;  G4double pQ6Rout =  25.4 * cm;  G4double pQ6HL   = 38.1  * cm;  G4double pQ6Pos_z=   0.0 * cm;
  G4VSolid* Q6MagSolid = new G4Tubs( "Q6MagTubs", pQ6Rin, pQ6Rout, pQ6HL, 0.0, 360.0 * deg);
  G4LogicalVolume* Q6MagLogical = new G4LogicalVolume(Q6MagSolid, Vacuum, "Q6Mag", 0,0,0);
  Q6MagLogical->SetFieldManager(Q6FieldManager, allLocal);
  Q6MagLogical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0, G4ThreeVector(0, 0, pQ6Pos_z), Q6MagLogical, "Q6Mag", world_log, 0, 0, fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //Solenoid 'Physical' Volume
  //TODO: Check this physical colume, this is the same as before
  G4double pHLMZRin = 5.10 * cm;   G4double pHLMZRout = 15.0 * cm;   G4double pHLMZHLZ = 5.0 * cm;
  G4double pHLMZ1Pos_X = 0.0 * cm;   G4double pHLMZ1Pos_Y = 0.0 * cm;   G4double pHLMZ1Pos_Z = -12.6 * cm;
  G4double pHLMZ2Pos_X = 0.0 * cm;   G4double pHLMZ2Pos_Y = 0.0 * cm;   G4double pHLMZ2Pos_Z = 12.6 * cm;

  G4VSolid* HLMZ1Solid = new G4Tubs( "HLMZ1Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );
  G4VSolid* HLMZ2Solid = new G4Tubs( "HLMZ2Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );

  G4LogicalVolume* HLMZ1Logical = new G4LogicalVolume(HLMZ1Solid, Cu, "Helmholtz1", 0, 0, 0);
  G4LogicalVolume* HLMZ2Logical = new G4LogicalVolume(HLMZ2Solid, Cu, "Helmholtz2", 0, 0, 0);
  HLMZ1Logical->SetVisAttributes(CuVisAtt);
  HLMZ2Logical->SetVisAttributes(CuVisAtt);

  new G4PVPlacement(0, G4ThreeVector( pHLMZ1Pos_X, pHLMZ1Pos_Y, pHLMZ1Pos_Z), HLMZ1Logical, "Helmholtz1", world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector( pHLMZ2Pos_X,  pHLMZ2Pos_Y, pHLMZ2Pos_Z), HLMZ2Logical, "Helmholtz2", world_log, 0, 0, fCheckOverlaps);
 

  ///////////////////////////////////////////////////////////////
  //Collimator Vacuum Can (after Los Alamos Quad)

  G4double CollVacRO = 17.78*cm;
  G4double CollVacRI = 16.51*cm;
  G4double CollVacH = 34.29 * cm;
  
  G4double CollVacPosX = 0.0 * cm;
  G4double CollVacPosY = 0.0 * cm;
  G4double CollVacPosZ = 196.85 * cm; //not sure what this actually is 
  
  G4VSolid*CollVac = new G4Tubs("Coll_Vac",CollVacRI,CollVacRO, CollVacH, 0.0, 360.0* deg);
  G4LogicalVolume*CollVacLogical = new G4LogicalVolume(CollVac, aluminum, "CollVacLogical",0,0,0);

  CollVacLogical -> SetVisAttributes(AlumVisAtt);

  new G4PVPlacement(0,G4ThreeVector(CollVacPosX,CollVacPosY,CollVacPosZ), CollVacLogical,"CollVacChamber",world_log, 0, 0, fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // MAGNETS (Los Alamos and Argonne Quads)

  //Los Alamos Quad
  G4double pQ1Rin  =  4.7625 * cm;  G4double pQ1Rout = 20.00 * cm;  G4double pQ1HL   = 18.29 * cm;  G4double pQ1Pos_Z= 85.0 * cm;//check z position

  //Argonne Quads
  G4double pQ2Rin  =  12.7 * cm;  G4double pQ2Rout = 86.36 * cm;  G4double pQ2Pos_Z= 297.8 * cm;
  G4double pQ3Pos_Z= 431.1 * cm;

  G4VSolid* Q1Solid = new G4Tubs( "Q1Tubs", pQ1Rin, pQ1Rout, pQ1HL, 0.0, 360.0 * deg );
  G4VSolid* Q2Solid = new G4Tubs( "Q2Tubs", pQ2Rin, pQ2Rout, pQ2HL, 0.0, 360.0 * deg );
  G4VSolid* Q3Solid = new G4Tubs( "Q3Tubs", pQ2Rin, pQ2Rout, pQ2HL, 0.0, 360.0 * deg );
  G4VSolid* Q1MagSolid = new G4Tubs( "Q1MagTubs", 0., pQ1Rin, pQ1HL, 0.0, 360.0 * deg );
  G4VSolid* Q2MagSolid = new G4Tubs( "Q2MagTubs", 0., pQ2Rin, pQ2HL, 0.0, 360.0 * deg );
  G4VSolid* Q3MagSolid = new G4Tubs( "Q3MagTubs", 0., pQ2Rin, pQ2HL, 0.0, 360.0 * deg );

  G4LogicalVolume* Q1Logical = new G4LogicalVolume(Q1Solid,siliconsteel,"Q1Logical",0,0,0);
  G4LogicalVolume* Q2Logical = new G4LogicalVolume(Q2Solid,siliconsteel,"Q2Logical",0,0,0);
  G4LogicalVolume* Q3Logical = new G4LogicalVolume(Q3Solid,siliconsteel,"Q3Logical",0,0,0);

  G4LogicalVolume* Q1MagLogical = new G4LogicalVolume(Q1MagSolid,Vacuum,"Q1MagLogical",0,0,0);
  G4LogicalVolume* Q2MagLogical = new G4LogicalVolume(Q2MagSolid,Vacuum,"Q2MagLogical",0,0,0);
  G4LogicalVolume* Q3MagLogical = new G4LogicalVolume(Q3MagSolid,Vacuum,"Q3MagLogical",0,0,0);

  Q1Logical->SetVisAttributes(SteelVisAtt);
  Q2Logical->SetVisAttributes(SteelVisAtt);
  Q3Logical->SetVisAttributes(SteelVisAtt);

  Q1MagLogical->SetVisAttributes(VacVisAtt);
  Q1MagLogical->SetFieldManager(Q1FieldManager,allLocal);

  Q2MagLogical->SetVisAttributes(VacVisAtt);
  Q2MagLogical->SetFieldManager(Q2FieldManager,allLocal);

  Q3MagLogical->SetVisAttributes(VacVisAtt);
  Q3MagLogical->SetFieldManager(Q3FieldManager,allLocal);

  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z),Q1Logical,"Q1Phys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z),Q2Logical,"Q2Phys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z),Q3Logical,"Q3Phys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z),Q1MagLogical,"Q1MagPhys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z),Q2MagLogical,"Q2MagPhys",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z),Q3MagLogical,"Q3MagPhys",world_log,0,0,fCheckOverlaps);

  ///////////////////////////////////////////////////////////////
  //Y-Vacuum Can

  //trapezoid dimensions
  G4double dx2 = 41.0 * cm;
  G4double dx1 = 22.5 * cm;
  G4double dy2 = 17.5 * cm;
  G4double dy1 = 17.5 * cm;
  G4double dz = 114.59 * cm;

  //inner trapezoid/vacuum dimensions
  G4double dx2_2 = 40.0 * cm;
  G4double dx1_2 = 21.5 * cm;
  G4double dy2_2 = 16.5 * cm;
  G4double dy1_2 = 16.5 * cm;

  G4double Trap_X = 0.0 * cm;
  G4double Trap_Y = 0.0 * cm;
  G4double Trap_Z = 618.9167 *cm;//633.2 * cm;//check

  G4VSolid*OuterTrap = new G4Trd("YCanOut", dx1,dx2,dy1,dy2,dz);
  G4VSolid*InnerTrap = new G4Trd("YCanIn", dx1_2,dx2_2,dy1_2,dy2_2,dz+1.0*cm);
  G4SubtractionSolid* YVacCan = new G4SubtractionSolid("YVacCan", OuterTrap, InnerTrap, 0, G4ThreeVector(Trap_X, Trap_Y, Trap_Z) );
  G4LogicalVolume* YVacCanLogical = new G4LogicalVolume( YVacCan, aluminum, "YVacCanLogical", 0, 0, 0);
  G4LogicalVolume* YVacCanInLogical  = new G4LogicalVolume( InnerTrap, Vacuum, "YVacCanInLogical", 0, 0, 0);
  YVacCanLogical->SetVisAttributes(AlumVisAtt);
  YVacCanInLogical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0,G4ThreeVector(Trap_X,Trap_Y,Trap_Z), YVacCanLogical,"YVacCanOut",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(Trap_X,Trap_Y,Trap_Z), YVacCanInLogical,"YVacCanIn",world_log, 0, 0, fCheckOverlaps);


  //Entrance flange  
  G4double YVOpeningRO = 17.8 * cm;
  G4double YVOpeningRI = 17.15 * cm; 
  G4double YVOpeningH = 4.25 * cm;

  G4double YVOpening_X = 0.0 * cm;
  G4double YVOpening_Y = 0.0 * cm;
  G4double YVOpening_Z = Trap_Z - dz -YVOpeningH;

  G4VSolid* YVOpeningOut = new G4Tubs( "YVOpenOut", YVOpeningRI, YVOpeningRO, YVOpeningH, 0.0, 360.0 * deg );
  //G4VSolid* YVOpeningIn = new G4Tubs( "YVOpenIn", ZeroR, YVOpeningRI, YVOpeningH, 0.0, 360.0 * deg );
  //G4SubtractionSolid* YVacOpen = new G4SubtractionSolid("YVacOpen", YVOpeningOut, YVOpeningIn, 0, G4ThreeVector(YVOpening_X, YVOpening_Y,YVOpening_Z) );
  G4LogicalVolume* YVacOpenLogical = new G4LogicalVolume( YVOpeningOut, aluminum, "YVacOpenLogical", 0, 0, 0);
  //G4LogicalVolume* YVacOpenInLogical = new G4LogicalVolume( YVOpeningIn, Vacuum, "YVacOpenInLogical", 0, 0, 0);
  YVacOpenLogical->SetVisAttributes(AlumVisAtt);
  // YVacOpenInLogical->SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0,G4ThreeVector(YVOpening_X,YVOpening_Y,YVOpening_Z), YVacOpenLogical,"YVacCanOpening",world_log, 0, 0, fCheckOverlaps);
  // new G4PVPlacement(0,G4ThreeVector(YVOpening_X,YVOpening_Y,YVOpening_Z), YVacOpenInLogical,"YVacCanOpeningInner",world_log, 0, 0, fCheckOverlaps);

  //Exit flanges: LRC - left rear cylinder; RRC - right rear cylinder; MRC - middle rear cylinder; RC - rear cylinder (general)
  G4double RCOR = 12.7 * cm;
  G4double RCIR = 12.5476 * cm;
  G4double RCH = 10.0 * cm;
 
  G4double RC_PosZ = Trap_Z + dz + RCH;
  G4double RC_X = 26.31 * cm;
  G4double RC_Y = 0.0 * cm;
  
    

  G4VSolid* LRC = new G4Tubs( "LRC", RCIR, RCOR, RCH, 0.0, 360.0 * deg);
  //G4VSolid* LRCInner = new G4Tubs( "LRCIn", ZeroR, RCIR, RCH, 0.0, 360.0 * deg); 
  //G4SubtractionSolid* LRC = new G4SubtractionSolid("LRC", LRCOuter, LRCInner, 0, G4ThreeVector(RC_X, RC_Y,RC_PosZ) );
  G4LogicalVolume* LRCLogical = new G4LogicalVolume( LRC, aluminum, "LRCLogical", 0, 0, 0);
  // G4LogicalVolume* LRCInLogical = new G4LogicalVolume( LRCInner, Vacuum, "LRCInLogical", 0, 0, 0);
  LRCLogical -> SetVisAttributes(AlumVisAtt);
  // LRCInLogical -> SetVisAttributes(VacVisAtt);
  new G4PVPlacement(RotPos,G4ThreeVector(-RC_X,RC_Y,RC_PosZ), LRCLogical,"LRC",world_log, 0, 0, fCheckOverlaps);
  // new G4PVPlacement(RotPos,G4ThreeVector(-RC_X,RC_Y,RC_PosZ), LRCInLogical,"LRCInner",world_log, 0, 0, fCheckOverlaps);


  G4VSolid* RRC = new G4Tubs( "RRC",RCIR, RCOR, RCH, 0.0, 360.0 * deg);
  //G4VSolid* RRCInner = new G4Tubs( "RRCIn", ZeroR, RCIR, RCH, 0.0, 360.0 * deg); 
  //G4SubtractionSolid* RRC = new G4SubtractionSolid("RRC", RRCOuter, RRCInner, 0, G4ThreeVector(RC_X, RC_Y,RC_PosZ) );
  G4LogicalVolume* RRCLogical = new G4LogicalVolume( RRC, aluminum, "RRCLogical", 0, 0, 0);
  //G4LogicalVolume* RRCInLogical = new G4LogicalVolume( RRCInner, Vacuum, "RRCInLogical", 0, 0, 0);
  RRCLogical -> SetVisAttributes(AlumVisAtt);
  //RRCInLogical -> SetVisAttributes(VacVisAtt);
  new G4PVPlacement(RotNeg,G4ThreeVector(RC_X,RC_Y,RC_PosZ), RRCLogical,"RRC",world_log, 0, 0, fCheckOverlaps);
  //new G4PVPlacement(RotNeg,G4ThreeVector(RC_X,RC_Y,RC_PosZ), RRCInLogical,"RRCInner",world_log, 0, 0, fCheckOverlaps);

  G4double MRCOR = 5.08 * cm;
  G4double MRCIR = 4.921 * cm;
  G4double MRCH = 6.0 * cm;
 
  G4double MRC_X = 0 * cm;
  G4double MRC_Y = 0.0 * cm;
  G4double MRC_Z = Trap_Z + dz + MRCH;

  G4VSolid* MRC = new G4Tubs( "MRCOut", MRCIR, MRCOR, MRCH, 0.0, 360.0 * deg);
  //G4VSolid* MRCInner = new G4Tubs( "MRCIn", ZeroR, MRCIR, MRCH, 0.0, 360.0 * deg); 
  // G4SubtractionSolid* MRC = new G4SubtractionSolid("MRC", MRCOuter, MRCInner, 0, G4ThreeVector(MRC_X, MRC_Y,MRC_Z) );
  G4LogicalVolume* MRCLogical = new G4LogicalVolume( MRC, aluminum, "MRCLogical", 0, 0, 0);
  //G4LogicalVolume* MRCInLogical = new G4LogicalVolume( MRCInner, Vacuum, "MRCInLogical", 0, 0, 0);
  MRCLogical -> SetVisAttributes(AlumVisAtt);
  //MRCInLogical -> SetVisAttributes(VacVisAtt);
  new G4PVPlacement(0,G4ThreeVector(MRC_X,MRC_Y,MRC_Z), MRCLogical,"MRC",world_log, 0, 0, fCheckOverlaps);
  // new G4PVPlacement(0,G4ThreeVector(MRC_X,MRC_Y,MRC_Z), MRCInLogical,"MRCInner",world_log, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Planes for Virtual Detectors
  G4VSolid*        VBSolid   = new G4Tubs("VBSolid",0,pQ1Rout, 0.00001 * mm ,0.0,360.0*deg);
  G4LogicalVolume* Q1ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q1ENLogical",0,0,0);
  G4LogicalVolume* Q1EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q1EXLogical",0,0,0);
  G4LogicalVolume* Q2ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q2ENLogical",0,0,0);
  G4LogicalVolume* Q2EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q2EXLogical",0,0,0);
  G4LogicalVolume* Q3ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q3ENLogical",0,0,0);
  G4LogicalVolume* Q3EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q3EXLogical",0,0,0);
 
  Q1ENLogical->SetVisAttributes(VacVisAtt);
  Q1EXLogical->SetVisAttributes(VacVisAtt);
  Q2ENLogical->SetVisAttributes(VacVisAtt);
  Q2EXLogical->SetVisAttributes(VacVisAtt);
  Q3ENLogical->SetVisAttributes(VacVisAtt);
  Q3EXLogical->SetVisAttributes(VacVisAtt);

  MolPolDetector* Q1ENSD = new MolPolDetector("q1en", 1);
  MolPolDetector* Q1EXSD = new MolPolDetector("q1ex", 2);
  MolPolDetector* Q2ENSD = new MolPolDetector("q2en", 3);
  MolPolDetector* Q2EXSD = new MolPolDetector("q2ex", 4);
  MolPolDetector* Q3ENSD = new MolPolDetector("q3en", 5);
  MolPolDetector* Q3EXSD = new MolPolDetector("q3ex", 6);

  MolPolDetector* DETSD  = new MolPolDetector("det",  9);
  MolPolDetector* DETSD2 = new MolPolDetector("det2", 10);
  MolPolDetector* APPSD1 = new MolPolDetector("app1", 11);
  MolPolDetector* APPSD2 = new MolPolDetector("app2", 12);
  MolPolDetector* DETVP  = new MolPolDetector("vp",   13);
  MolPolDetector* DPIN   = new MolPolDetector("dpin", 14);
  MolPolDetector* DPOUT  = new MolPolDetector("dpout",15);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  SDman->AddNewDetector(Q1ENSD);
  SDman->AddNewDetector(Q1EXSD);
  SDman->AddNewDetector(Q2ENSD);
  SDman->AddNewDetector(Q2EXSD);
  SDman->AddNewDetector(Q3ENSD);
  SDman->AddNewDetector(Q3EXSD);
  SDman->AddNewDetector(DETSD );
  SDman->AddNewDetector(DETSD2);
  SDman->AddNewDetector(APPSD1);
  SDman->AddNewDetector(APPSD2);
  SDman->AddNewDetector(DETVP );
  SDman->AddNewDetector(DPIN  );
  SDman->AddNewDetector(DPOUT );

  Q1ENLogical->SetSensitiveDetector(Q1ENSD);
  Q1EXLogical->SetSensitiveDetector(Q1EXSD);
  Q2ENLogical->SetSensitiveDetector(Q2ENSD);
  Q2EXLogical->SetSensitiveDetector(Q2EXSD);
  Q3ENLogical->SetSensitiveDetector(Q3ENSD);
  Q3EXLogical->SetSensitiveDetector(Q3EXSD);


  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z - pQ1HL ), Q1ENLogical,"VP.Q1.Entr",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z + pQ1HL ), Q1EXLogical,"VP.Q1.Exit",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z - pQ2HL ), Q2ENLogical,"VP.Q2.Entr",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z + pQ2HL ), Q2EXLogical,"VP.Q2.Exit",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z - pQ2HL ), Q3ENLogical,"VP.Q3.Entr",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z + pQ2HL ), Q3EXLogical,"VP.Q3.Exit",world_log,0,0,fCheckOverlaps);

  ///////////////////////////////////////////////////////////////
  //Collimaters

  //Rotation Matrix
  G4RotationMatrix* LeftCollRot = new G4RotationMatrix();
  LeftCollRot->rotateY(86.75*deg);
  LeftCollRot->rotateZ(90.0*deg);

  G4RotationMatrix* RightCollRot = new G4RotationMatrix();
  RightCollRot->rotateY(93.25*deg);
  RightCollRot->rotateZ(90.0*deg);  

  //Collimater Box Dimmensions
  G4double Box_Y = 5.08 * cm; //Coll width 
  G4double Box_X = 7.0 * cm;
  G4double CollThick = 10.0 * cm;

  //LCT = Left Collimater trapezoid
  G4double LCT_dx2 = 2.80 * cm;
  G4double LCT_dx1 = 4.00 * cm;
  G4double LCT_dy = 5.08 * cm;
  G4double LCT_dz = 7.2 * cm;

  //RCT = Right Collimater trapezoid
  G4double RCT_dx1 = 2.0 * cm;
  G4double RCT_dx2 = 3.0*cm;
  G4double RCT_dy = 5.08 *cm;
  G4double RCT_dz = 7.0*cm;

  //Placement
  G4double Coll_PosZ = 1124.08 *cm; 
  G4double Coll_PosX = 49.00 * cm;
  G4double Coll_PosY = 0.0 * cm;

  G4VSolid*CollBox = new G4Box("CollBox",Box_X, Box_Y,CollThick);
  G4VSolid*LCT = new G4Trd("LCT", LCT_dx1, LCT_dx2, LCT_dy, LCT_dy,LCT_dz);
  G4VSolid*RCT = new G4Trd("RCT", RCT_dx1,RCT_dx2,RCT_dy,RCT_dy,RCT_dz);

  G4SubtractionSolid* LeftCollimater = new G4SubtractionSolid("Left_Collimater", CollBox, LCT, 0, G4ThreeVector(Coll_PosX, Coll_PosY, Coll_PosZ) );
  G4SubtractionSolid* RightCollimater = new G4SubtractionSolid("Right_Collimater", CollBox, RCT, 0, G4ThreeVector(-Coll_PosX, Coll_PosY, Coll_PosZ) );
  
  G4LogicalVolume*LeftCollLogic = new G4LogicalVolume(LeftCollimater, LgTF1, "LeftColl_Logical",0,0,0);
  G4LogicalVolume*LCT_Logic = new G4LogicalVolume(LCT, Vacuum, "LCT_Logical",0,0,0);

  G4LogicalVolume*RightCollLogic = new G4LogicalVolume(RightCollimater, LgTF1, "RightColl_Logical",0,0,0);
  G4LogicalVolume*RCT_Logic = new G4LogicalVolume(RCT, Vacuum, "RCT_Logical",0,0,0);

  LeftCollLogic->SetVisAttributes(AlumVisAtt);
  RightCollLogic->SetVisAttributes(AlumVisAtt);
  LCT_Logic->SetVisAttributes(IronVisAtt);
  RCT_Logic->SetVisAttributes(IronVisAtt);
  
  new G4PVPlacement(LeftCollRot,G4ThreeVector(Coll_PosX, Coll_PosY, Coll_PosZ), LeftCollLogic,"Left Collimater",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(LeftCollRot,G4ThreeVector(Coll_PosX, Coll_PosY, Coll_PosZ), LCT_Logic,"Inner Left Collimater",world_log, 0, 0, fCheckOverlaps);  

  new G4PVPlacement(RightCollRot,G4ThreeVector(-Coll_PosX, Coll_PosY, Coll_PosZ), RightCollLogic,"Right Collimater",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(RightCollRot,G4ThreeVector(-Coll_PosX, Coll_PosY, Coll_PosZ), RCT_Logic,"Inner Right Collimater",world_log, 0, 0, fCheckOverlaps);  


 
  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    //Lead Glass Detectors
    G4double LG_PosZ  = 1140.66 * cm;
    G4double LG_Y = 7.00 * cm;
    G4double LG_X  = 10.00 * cm;
    G4double LG_Z = 11.5 * cm;
    G4double LG_PosY = 0 * cm;  
 
    G4double LG_PosX   = 49.941 * cm; 
  
    G4VSolid* MDBXSolid3  = new G4Box ( "MDBXBox3 "  , LG_X, LG_Y, LG_Z );
    G4LogicalVolume* DETLogical3 = new G4LogicalVolume(MDBXSolid3, LgTF1, "DETLogical3",0,0,0);
    DETLogical3->SetVisAttributes(CuVisAtt);
    MolPolDetector* MDBXBox3 = new MolPolDetector("PbDetector",20);
    SDman->AddNewDetector(MDBXBox3);
    DETLogical3->SetSensitiveDetector(MDBXBox3);

    new G4PVPlacement(RotNeg,G4ThreeVector(LG_PosX, LG_PosY, LG_PosZ),DETLogical3,"LgDet1",world_log,0,0,fCheckOverlaps);

    //Second detector
   
    new G4PVPlacement(RotPos,G4ThreeVector(-LG_PosX, LG_PosY, LG_PosZ),DETLogical3,"LgDet2",world_log,0,0,fCheckOverlaps);
 
//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    // HODOSCOPE 1

 //X, Y and Z Dimmensions
	
    G4double HODOX_DIM = 0.6 * cm;
    G4double HODOY_DIM = 4.0 * cm;
    G4double HODOZ_DIM = 0.4 * cm;
    G4double HODOSplitY = 2.0 * cm;
	
    //X positions of each hodoscope segment
    //The middle of the hodoscope is located at x= 20.411cm 
    //The first hodoscope to the left (L1) of the center is centered at x=41.391+0.6, to the right of the center (R1) x=41.391-0.6
    //The formula to find the center point of the nth hodoscope to the left: 21.068 + 1.2*(n-1)
    // For the right: 19.868 - 1.2*(n-1)
    
    G4double distBetween = 2.0*HODOX_DIM/cm; 

    //TODO: What is the z position??? For now assumption is 1109.0cm
    G4double HODO_ZPOS = (Coll_PosZ/cm-Box_Y/cm-HODOZ_DIM/cm);//cm
    G4double L1_Z = HODO_ZPOS - HODOX_DIM/cm*sin(theta); //initial offset because there is no centeral box (#of boxes is an even number)
    G4double R1_Z = HODO_ZPOS + HODOX_DIM/cm*sin(theta);
    // G4double sinCheck = sin(Htheta);
    
    G4double L1 = 48.941 + HODOX_DIM/cm;
    G4double R1 = 48.941 - HODOX_DIM/cm;

    G4double L1X   = L1 * cm;                       G4double L6X = (L1+distBetween*(6-1)) * cm;    G4double R4X   = (R1-distBetween*(4-1)) * cm;
    G4double L2X   = (L1 +distBetween*(2-1)) * cm;  G4double L7X   = (L1+distBetween*(7-1)) * cm;  G4double R5X   = (R1-distBetween*(5-1)) * cm;
    G4double L3X   = (L1+distBetween*(3-1)) * cm;   G4double R1X   = R1 * cm;                      G4double R6X   = (R1-distBetween*(6-1)) * cm;
    G4double L4X   = (L1+distBetween*(4-1)) * cm;   G4double R2X   = (R1-distBetween*(2-1)) * cm;  G4double R7X   = (R1-distBetween*(7-1)) * cm;	
    G4double L5X   = (L1+distBetween*(5-1)) * cm;   G4double R3X   = (R1-distBetween*(3-1)) * cm;  
  
	
    //Y positions
    G4double HODOY   = 0.0 * cm;
    G4double HODOSplitNegY = -2.0 * cm;
    G4double HODOSplitPosY = 2.0 * cm;

    //Z positions
    //Z position of the middle of the entire hodoscope: 1119.8 cm
    //nth box has a z-offset of 1.2*sin(3.25)*(n-1)
    // z position of nth left box = 1119.8 - HODOX_DIM*sin(pRotH1)*(n-1)
    //'' '' ''  nth right box = 1119.8 + HODOX_DIM*sin(pRotH1)*(n-1)

  
    
    G4double L1Z = L1_Z * cm;                                  
    G4double L2Z = (L1_Z - (distBetween*sin(theta)*(2-1))) * cm;      
    G4double L3Z = (L1_Z - (distBetween*sin(theta)*(3-1))) * cm;      
    G4double L4Z = (L1_Z - (distBetween*sin(theta)*(4-1))) * cm;      
    G4double L5Z = (L1_Z - (distBetween*sin(theta)*(5-1))) * cm; 
    G4double L6Z = (L1_Z - (distBetween*sin(theta)*(6-1))) * cm;
    G4double L7Z = (L1_Z - (distBetween*sin(theta)*(7-1))) * cm;
  
    G4double R1Z = R1_Z * cm;
    G4double R2Z = (R1_Z + (distBetween*sin(theta)*(2-1))) * cm;
    G4double R3Z = (R1_Z + (distBetween*sin(theta)*(3-1))) * cm;
    G4double R4Z = (R1_Z + (distBetween*sin(theta)*(4-1))) * cm;
    G4double R5Z = (R1_Z + (distBetween*sin(theta)*(5-1))) * cm;
    G4double R6Z = (R1_Z + (distBetween*sin(theta)*(6-1))) * cm;
    G4double R7Z = (R1_Z + (distBetween*sin(theta)*(7-1))) * cm;

    //Creating the 14 hodoscope boxes
    G4VSolid* HODOBOX1  = new G4Box ( "HODOBOX1 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX2  = new G4Box ( "HODOBOX2 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX3  = new G4Box ( "HODOBOX3 "  , HODOX_DIM, HODOSplitY, HODOZ_DIM );
    G4VSolid* HODOBOX3Split = new G4Box ( "HODOBOX3Split" , HODOX_DIM, HODOSplitY,  HODOZ_DIM);
    G4VSolid* HODOBOX4  = new G4Box ( "HODOBOX4 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX5  = new G4Box ( "HODOBOX5 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX6  = new G4Box ( "HODOBOX6 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX7  = new G4Box ( "HODOBOX7 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX8  = new G4Box ( "HODOBOX8 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX9  = new G4Box ( "HODOBOX9 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX10  = new G4Box ( "HODOBOX10 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX11  = new G4Box ( "HODOBOX11 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX12  = new G4Box ( "HODOBOX12 "  , HODOX_DIM, HODOSplitY, HODOZ_DIM );
    G4VSolid* HODOBOX12Split = new G4Box ( "HODOBOX12Split" , HODOX_DIM, HODOSplitY,  HODOZ_DIM);
    G4VSolid* HODOBOX13  = new G4Box ( "HODOBOX13 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
    G4VSolid* HODOBOX14  = new G4Box ( "HODOBOX14 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	
    //Giving each box a volume	
    G4LogicalVolume* HODO1Logical = new G4LogicalVolume(HODOBOX1, scint, "HODO1Logical",0,0,0);
    HODO1Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO1 = new MolPolDetector("HODO1",25);
    SDman -> AddNewDetector(HODO1);
    HODO1Logical -> SetSensitiveDetector(HODO1);

    G4LogicalVolume* HODO2Logical = new G4LogicalVolume(HODOBOX2, scint, "HODO2",0,0,0);
    HODO2Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO2 = new MolPolDetector("HODO2",26);
    SDman -> AddNewDetector(HODO2);
    HODO2Logical -> SetSensitiveDetector(HODO2);

    G4LogicalVolume* HODO3Logical = new G4LogicalVolume(HODOBOX3, scint, "HODO3",0,0,0);
    HODO3Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO3 = new MolPolDetector("HODO3",27);
    SDman -> AddNewDetector(HODO3);
    HODO3Logical -> SetSensitiveDetector(HODO3); 
 
    G4LogicalVolume*HODO3SplitLogical = new G4LogicalVolume(HODOBOX3Split,scint, "HODO3Split",0,0,0);
    HODO3SplitLogical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO3Split = new MolPolDetector("HODO3Split",28);
    SDman -> AddNewDetector(HODO3Split);
    HODO3SplitLogical -> SetSensitiveDetector(HODO3Split);

    G4LogicalVolume* HODO4Logical = new G4LogicalVolume(HODOBOX4, scint, "HODO4",0,0,0);
    HODO4Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO4 = new MolPolDetector("HODO4",29);
    SDman -> AddNewDetector(HODO4);
    HODO4Logical -> SetSensitiveDetector(HODO4);

    G4LogicalVolume* HODO5Logical = new G4LogicalVolume(HODOBOX5, scint, "HODO5",0,0,0);
    HODO5Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO5 = new MolPolDetector("HODO5",30);
    SDman -> AddNewDetector(HODO5);
    HODO5Logical -> SetSensitiveDetector(HODO5);
	
    G4LogicalVolume* HODO6Logical = new G4LogicalVolume(HODOBOX6, scint, "HODO6",0,0,0);
    HODO6Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO6 = new MolPolDetector("HODO6",31);
    SDman -> AddNewDetector(HODO6);
    HODO1Logical -> SetSensitiveDetector(HODO6);
	
    G4LogicalVolume* HODO7Logical = new G4LogicalVolume(HODOBOX7, scint, "HODO7",0,0,0);
    HODO7Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO7 = new MolPolDetector("HODO7",32);
    SDman -> AddNewDetector(HODO7);
    HODO7Logical -> SetSensitiveDetector(HODO7);
	
    G4LogicalVolume* HODO8Logical = new G4LogicalVolume(HODOBOX8, scint, "HODO8",0,0,0);
    HODO8Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO8 = new MolPolDetector("HODO8",33);
    SDman -> AddNewDetector(HODO8);
    HODO8Logical -> SetSensitiveDetector(HODO8);

    G4LogicalVolume* HODO9Logical = new G4LogicalVolume(HODOBOX9, scint, "HODO9",0,0,0);
    HODO9Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO9 = new MolPolDetector("HODO9",34);
    SDman -> AddNewDetector(HODO9);
    HODO9Logical -> SetSensitiveDetector(HODO9);

    G4LogicalVolume* HODO10Logical = new G4LogicalVolume(HODOBOX10, scint, "HODO10",0,0,0);
    HODO10Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO10 = new MolPolDetector("HODO10",35);
    SDman -> AddNewDetector(HODO10);
    HODO10Logical -> SetSensitiveDetector(HODO10);
	
    G4LogicalVolume* HODO11Logical = new G4LogicalVolume(HODOBOX11, scint, "HODO11",0,0,0);
    HODO11Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO11 = new MolPolDetector("HODO11",36);
    SDman -> AddNewDetector(HODO11);
    HODO11Logical -> SetSensitiveDetector(HODO11);

    G4LogicalVolume* HODO12Logical = new G4LogicalVolume(HODOBOX12, scint, "HODO12",0,0,0);
    HODO12Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO12 = new MolPolDetector("HODO12",37);
    SDman -> AddNewDetector(HODO12);
    HODO12Logical -> SetSensitiveDetector(HODO12);

    G4LogicalVolume*HODO12SplitLogical = new G4LogicalVolume(HODOBOX12Split,scint, "HODO12Split",0,0,0);
    HODO12SplitLogical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO12Split = new MolPolDetector("HODO12Split",38);
    SDman -> AddNewDetector(HODO12Split);
    HODO12SplitLogical -> SetSensitiveDetector(HODO12Split);

    G4LogicalVolume* HODO13Logical = new G4LogicalVolume(HODOBOX13, scint, "HODO13",0,0,0);
    HODO13Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO13 = new MolPolDetector("HODO13",39);
    SDman -> AddNewDetector(HODO13);
    HODO13Logical -> SetSensitiveDetector(HODO13);

    G4LogicalVolume* HODO14Logical = new G4LogicalVolume(HODOBOX14, scint, "HODO14",0,0,0);
    HODO14Logical-> SetVisAttributes(ScintVisAtt);
    MolPolDetector*HODO14 = new MolPolDetector("HODO14",40);
    SDman -> AddNewDetector(HODO14);
    HODO14Logical -> SetSensitiveDetector(HODO14);				
		
    	
    new G4PVPlacement(RotPos,G4ThreeVector(L7X, HODOY, L7Z),HODO1Logical,"HODO_1",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(L6X, HODOY, L6Z),HODO2Logical,"HODO_2",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(L5X, HODOSplitPosY, L5Z),HODO3Logical,"HODO_3",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(L5X, HODOSplitNegY, L5Z),HODO3SplitLogical,"HODO_3Lower",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(L4X, HODOY, L4Z),HODO4Logical,"HODO_4",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(L3X, HODOY, L3Z),HODO5Logical,"HODO_5",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(L2X, HODOY, L2Z),HODO6Logical,"HODO_6",world_log,0,0,fCheckOverlaps);	
    new G4PVPlacement(RotPos,G4ThreeVector(L1X, HODOY, L1Z),HODO7Logical,"HODO_7",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R1X, HODOY, R1Z),HODO8Logical,"HODO_8",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R2X, HODOY, R2Z),HODO9Logical,"HODO_9",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R3X, HODOY, R3Z),HODO10Logical,"HODO_10",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R4X, HODOY, R4Z),HODO11Logical,"HODO_11",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R5X, HODOSplitPosY, R5Z),HODO12Logical,"HODO_12",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R5X, HODOSplitNegY, R5Z),HODO12SplitLogical,"HODO_12Lower",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R6X, HODOY, R6Z),HODO13Logical,"HODO_13",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(RotPos,G4ThreeVector(R7X, HODOY, R7Z),HODO14Logical,"HODO_14",world_log,0,0,fCheckOverlaps);
   

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    // HODOSCOPE 2
	
   
 

	//Creating the 14 hodoscope boxes
        G4VSolid* HODOBOX012  = new G4Box ( "HODOBOX012 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX22  = new G4Box ( "HODOBOX22 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX32  = new G4Box ( "HODOBOX32 "  , HODOX_DIM, HODOSplitY, HODOZ_DIM );
	G4VSolid* HODOBOX32Split = new G4Box ( "HODOBOX32Split" , HODOX_DIM, HODOSplitY,  HODOZ_DIM);
	G4VSolid* HODOBOX42  = new G4Box ( "HODOBOX42"  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX52  = new G4Box ( "HODOBOX52 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX62  = new G4Box ( "HODOBOX62 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX72  = new G4Box ( "HODOBOX72 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX82  = new G4Box ( "HODOBOX82 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX92  = new G4Box ( "HODOBOX92 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX102  = new G4Box ( "HODOBOX102"  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX112  = new G4Box ( "HODOBOX112 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX122  = new G4Box ( "HODOBOX122 "  , HODOX_DIM, HODOSplitY, HODOZ_DIM );
	G4VSolid* HODOBOX122Split = new G4Box ( "HODOBOX122Split" , HODOX_DIM, HODOSplitY,  HODOZ_DIM);
	G4VSolid* HODOBOX132  = new G4Box ( "HODOBOX132 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	G4VSolid* HODOBOX142  = new G4Box ( "HODOBOX142 "  , HODOX_DIM, HODOY_DIM, HODOZ_DIM );
	
	
	G4LogicalVolume* HODO012Logical = new G4LogicalVolume(HODOBOX012, scint, "HODO012",0,0,0);
	HODO012Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO012 = new MolPolDetector("HODO012",41);
	SDman -> AddNewDetector(HODO012);
	HODO012Logical -> SetSensitiveDetector(HODO012);        

	G4LogicalVolume* HODO22Logical = new G4LogicalVolume(HODOBOX22, scint, "HODO22",0,0,0);
        HODO22Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO22 = new MolPolDetector("HODO22",42);
	SDman -> AddNewDetector(HODO22);
	HODO22Logical -> SetSensitiveDetector(HODO22); 

	G4LogicalVolume* HODO32Logical = new G4LogicalVolume(HODOBOX32, scint, "HODO32",0,0,0);
        HODO32Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO32 = new MolPolDetector("HODO32",43);
	SDman -> AddNewDetector(HODO32);
	HODO32Logical -> SetSensitiveDetector(HODO32);

	G4LogicalVolume* HODO32SplitLogical = new G4LogicalVolume(HODOBOX32Split,scint, "HODO32Split",0,0,0);
	HODO32SplitLogical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO32Split = new MolPolDetector("HODO32Split",44);
	SDman -> AddNewDetector(HODO32Split);
	HODO32SplitLogical -> SetSensitiveDetector(HODO32Split);

	G4LogicalVolume* HODO42Logical = new G4LogicalVolume(HODOBOX42, scint, "HODO42",0,0,0);
	HODO42Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO42 = new MolPolDetector("HODO42",45);
	SDman -> AddNewDetector(HODO42);
	HODO42Logical -> SetSensitiveDetector(HODO42);        

	G4LogicalVolume* HODO52Logical = new G4LogicalVolume(HODOBOX52, scint, "HODO52",0,0,0);
        HODO52Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO52 = new MolPolDetector("HODO52",46);
	SDman -> AddNewDetector(HODO52);
	HODO52Logical -> SetSensitiveDetector(HODO52);
	
	G4LogicalVolume* HODO62Logical = new G4LogicalVolume(HODOBOX62, scint, "HODO62",0,0,0);
	HODO62Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO62 = new MolPolDetector("HODO62",47);
	SDman -> AddNewDetector(HODO62);
	HODO62Logical -> SetSensitiveDetector(HODO62);
	
	G4LogicalVolume* HODO72Logical = new G4LogicalVolume(HODOBOX72, scint, "HODO72",0,0,0);
	HODO72Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO72 = new MolPolDetector("HODO72",48);
	SDman -> AddNewDetector(HODO72);
	HODO72Logical -> SetSensitiveDetector(HODO72);
	
	G4LogicalVolume* HODO82Logical = new G4LogicalVolume(HODOBOX82, scint, "HODO82",0,0,0);
	HODO82Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO82 = new MolPolDetector("HODO82",49);
	SDman -> AddNewDetector(HODO82);
	HODO82Logical -> SetSensitiveDetector(HODO82);

	G4LogicalVolume* HODO92Logical = new G4LogicalVolume(HODOBOX92, scint, "HODO92",0,0,0);
	HODO92Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO92 = new MolPolDetector("HODO92",50);
	SDman -> AddNewDetector(HODO92);
	HODO92Logical -> SetSensitiveDetector(HODO92);

	G4LogicalVolume* HODO102Logical = new G4LogicalVolume(HODOBOX102, scint, "HODO102",0,0,0);
	HODO102Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO102 = new MolPolDetector("HODO102",51);
	SDman -> AddNewDetector(HODO102);
	HODO102Logical -> SetSensitiveDetector(HODO102);
		
	G4LogicalVolume* HODO112Logical = new G4LogicalVolume(HODOBOX112, scint, "HODO112",0,0,0);
	HODO112Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO112 = new MolPolDetector("HODO112",52);
	SDman -> AddNewDetector(HODO112);
	HODO112Logical -> SetSensitiveDetector(HODO112);

	G4LogicalVolume* HODO122Logical = new G4LogicalVolume(HODOBOX122, scint, "HODO122",0,0,0);
        HODO122Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO122 = new MolPolDetector("HODO122",53);
	SDman -> AddNewDetector(HODO122);
	HODO122Logical -> SetSensitiveDetector(HODO122);

	G4LogicalVolume*HODO122SplitLogical = new G4LogicalVolume(HODOBOX122Split,scint, "HODO122Split",0,0,0);
	HODO122SplitLogical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO122Split = new MolPolDetector("HODO122Split",54);
	SDman -> AddNewDetector(HODO122Split);
	HODO122SplitLogical -> SetSensitiveDetector(HODO122Split);

	G4LogicalVolume* HODO132Logical = new G4LogicalVolume(HODOBOX132, scint, "HODO132",0,0,0);
	HODO132Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO132 = new MolPolDetector("HODO132",55);
	SDman -> AddNewDetector(HODO132);
	HODO132Logical -> SetSensitiveDetector(HODO132);

	G4LogicalVolume* HODO142Logical = new G4LogicalVolume(HODOBOX142, scint, "HODO142",0,0,0);
	HODO142Logical-> SetVisAttributes(ScintVisAtt);
	MolPolDetector*HODO142 = new MolPolDetector("HODO142",56);
	SDman -> AddNewDetector(HODO142);
	HODO142Logical -> SetSensitiveDetector(HODO142);        
    
		
		
	new G4PVPlacement(RotNeg,G4ThreeVector(-L7X, HODOY, L7Z),HODO012Logical,"HODO_012",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-L6X, HODOY, L6Z),HODO22Logical,"HODO2",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-L5X, HODOSplitPosY, L5Z),HODO32Logical,"HODO_32",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-L5X, HODOSplitNegY, L5Z),HODO32SplitLogical,"HODO_32Lower",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-L4X, HODOY, L4Z),HODO42Logical,"HODO_42",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-L3X, HODOY, L3Z),HODO52Logical,"HODO_52",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-L2X, HODOY, L2Z),HODO62Logical,"HODO_62",world_log,0,0,fCheckOverlaps);	
	new G4PVPlacement(RotNeg,G4ThreeVector(-L1X, HODOY, L1Z),HODO72Logical,"HODO_72",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R1X, HODOY, R1Z),HODO82Logical,"HODO_82",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R2X, HODOY, R2Z),HODO92Logical,"HODO_92",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R3X, HODOY, R3Z),HODO102Logical,"HODO_102",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R4X, HODOY, R4Z),HODO112Logical,"HODO_112",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R5X, HODOSplitPosY, R5Z),HODO122Logical,"HODO_122",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R5X, HODOSplitNegY, R5Z),HODO32SplitLogical,"HODO_122Lower",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R6X, HODOY, R6Z),HODO132Logical,"HODO_132",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(RotNeg,G4ThreeVector(-R7X, HODOY, R7Z),HODO142Logical,"HODO_142",world_log,0,0,fCheckOverlaps);
	
   

 ///////////////////////////////////////////////////////////////
  //Beam Pipes after Y-vacuum can
     

  G4double BPRO1 = 12.70 * cm;
  G4double BPRI1 = 12.5476 * cm;
  G4double BPH1 = (HODO_ZPOS*cm-HODOZ_DIM-(RC_PosZ+RCH))/2;//written to cover all empty space between Y-vacuum can and hodoscope

  G4double BP1_Y = 0.0 * cm;
  G4double BP1_Z = (RC_PosZ/cm+RCH/cm+BPH1/cm)*cm;
  G4double BP1_X = (tan(theta)*(BPH1/cm+11.5+0.8))* cm - 49.0 * cm; //have to check

  G4VSolid*B_Pipe = new G4Tubs("B_Pipe1", BPRI1, BPRO1, BPH1, 0.0, 360.0*deg);
  //G4VSolid*BPipe2 = new G4Tubs("B_PipeIn", ZeroR, BPRI1, BPH1, 0.0, 360.0*deg);
  //G4SubtractionSolid* B_Pipe = new G4SubtractionSolid("B_Pipe", BPipe1, BPipe2, 0, G4ThreeVector(BP1_X, BP1_Y, BP1_Z) );
  G4LogicalVolume* BPipe1Logical = new G4LogicalVolume( B_Pipe, aluminum, "BPipeOutLogical", 0, 0, 0);
  // G4LogicalVolume* BPipe2Logical  = new G4LogicalVolume( BPipe2, Vacuum, "BPipeInLogical", 0, 0, 0);
  BPipe1Logical->SetVisAttributes(AlumVisAtt);
  // BPipe2Logical->SetVisAttributes(VacVisAtt);
  
  new G4PVPlacement(RotNeg,G4ThreeVector(-BP1_X,BP1_Y,BP1_Z), BPipe1Logical,"BPipe1Out",world_log, 0, 0, fCheckOverlaps);
  // new G4PVPlacement(RotNeg,G4ThreeVector(-BP1_X,BP1_Y,BP1_Z), BPipe2Logical,"BPipe1In",world_log, 0, 0, fCheckOverlaps);

  new G4PVPlacement(RotPos,G4ThreeVector(BP1_X,BP1_Y,BP1_Z), BPipe1Logical,"BPipe2Out",world_log, 0, 0, fCheckOverlaps);
  // new G4PVPlacement(RotPos,G4ThreeVector(BP1_X,BP1_Y,BP1_Z), BPipe2Logical,"BPipe2In",world_log, 0, 0, fCheckOverlaps);

      


  return world_phys;
}

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //
