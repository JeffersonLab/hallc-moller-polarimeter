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
  G4VisAttributes* LeadVisAtt = new G4VisAttributes( G4Colour(149./255.,149./255.,100./255.,alphaMatStd) );
  G4VisAttributes* SteelVisAtt= new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaMatStd) );
  G4VisAttributes* AlumVisAtt = new G4VisAttributes( G4Colour(  0./255.,237./255.,  0./255.,alphaMatStd) );
  G4VisAttributes* VacVisAtt  = new G4VisAttributes( G4Colour(255./255.,255./255.,255./255.,alphaVacuum) );
  G4VisAttributes* CuVisAtt   = new G4VisAttributes( G4Colour(178./255.,102./255., 26./255.,alphaMatStd) );
  G4VisAttributes* ScintVisAtt= new G4VisAttributes( G4Colour(  0./255.,100./255.,100./255.,alphaMatStd) );
  G4VisAttributes* DipVisAtt  = new G4VisAttributes( G4Colour(  0./255., 80./255.,225./255.,alphaVacuum) );


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //// Materials Setup
  G4double a, z, density, pressure, temperature;
  G4int nelements, natoms;

  G4Element* N  = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* H  = new G4Element("Hydrogen", "H", z=1 , a=1.01 *g/mole);
  G4Element* C  = new G4Element("Carbon",   "C", z=6, a=12.01*g/mole);
  G4Element* As  = new G4Element("Arsen","As", z=33, a=74.92*g/mole);

  G4Element* Al = new G4Element("Aluminum", "Al", z=13, a=26.98*g/mole);
  G4Element* Fe = new G4Element("Iron"   , "Fe", z=26, a=55.845*g/mole);
  G4Element* Si = new G4Element("Silicon", "Si", z=14, a=28.09 *g/mole);
  G4Element* Pb = new G4Element("Lead", "Pb", z=82., a=207.19*g/mole);

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
  G4Material* lead = new G4Material("lead", z=82, a, density);

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

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Field Setup
  mEMFieldSetup = new MolPolEMFieldSetup();  //setup the field,
  G4FieldManager* Q1FieldManager = mEMFieldSetup->GetFieldManagerFZB1();
  G4FieldManager* Q2FieldManager = mEMFieldSetup->GetFieldManagerFZB2();
  G4FieldManager* Q3FieldManager = mEMFieldSetup->GetFieldManagerFZB3();
  G4FieldManager* Q4FieldManager = mEMFieldSetup->GetFieldManagerFZB4();
  G4FieldManager* DFieldManager  = mEMFieldSetup->GetFieldManagerFZB5();
  G4FieldManager* Q6FieldManager = mEMFieldSetup->GetFieldManagerFZB6();
  G4bool allLocal = true;

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Target
  G4double pMTATRin   = 0.0 * cm; G4double pMTATRout  = 1.5 * cm;   G4double pMTATHLZ = 0.0062 * cm;
  G4double pMTATPos_X = 0.0 * cm; G4double pMTATPos_Y = 0.0 * cm; G4double pMTATPos_Z = 0.0 * cm;
  G4VSolid* MTATSolid = new G4Tubs( "MTATTube", pMTATRin, pMTATRout, pMTATHLZ, 0.0, 360.0 * deg );
  G4LogicalVolume* MTATLogical = new G4LogicalVolume(MTATSolid, iron, "Target", 0, 0, 0);
  MTATLogical->SetVisAttributes(IronVisAtt);
  new G4PVPlacement(0, G4ThreeVector(0,0,pMTATPos_Z), MTATLogical, "Target", world_log, 0, 0, fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Target BPIPE
   G4double pBPITRin = 0.0 * cm;   G4double pBPITRout = 5.08 * cm;   G4double pBPITHLZ = 100.0 * cm;
   G4double pBPVTRin = 0.0 * cm;   G4double pBPVTRout = 4.78 * cm;   G4double pBPVTHLZ = 100.0 * cm;

  G4double pBPITPos_X = 0.0 * cm;   G4double pBPITPos_Y = 0.0 * cm;   G4double pBPITPos_Z = 0.0 * cm;
  G4double pBPVTPos_X = 0.0 * cm;   G4double pBPVTPos_Y = 0.0 * cm;   G4double pBPVTPos_Z = 0.0 * cm;

  G4VSolid* BPITSolid = new G4Tubs( "BPITTube", pBPVTRout, pBPITRout, pBPITHLZ, 0.0, 360.0 * deg );
  G4VSolid* BPVTSolid = new G4Tubs( "BPVTTube", pBPVTRin , pBPVTRout, pBPITHLZ, 0.0, 360.0 * deg );

  G4LogicalVolume* BPITLogical = new G4LogicalVolume(BPITSolid, aluminum, "BPITLogical", 0, 0, 0);
  G4LogicalVolume* BPVTLogical = new G4LogicalVolume(BPVTSolid, Vacuum,   "BPVTLogical", 0, 0, 0);
  BPITLogical->SetVisAttributes(AlumVisAtt);
  BPVTLogical->SetVisAttributes(VacVisAtt);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // HELMHOLTZ COIL Magnetic Volume
  G4double pQ6Rin  =  0    * cm;  G4double pQ6Rout =  25.4 * cm;  G4double pQ6HL   = 38.1  * cm;  G4double pQ6Pos_z=   -500 * cm;
  G4VSolid* Q6MagSolid = new G4Tubs( "Q6MagTubs", pQ6Rin, pQ6Rout, pQ6HL, 0.0, 360.0 * deg);
  G4LogicalVolume* Q6MagLogical = new G4LogicalVolume(Q6MagSolid, Vacuum, "Q6Mag", 0,0,0);
  // Q6MagLogical->SetFieldManager(Q6FieldManager, allLocal);
  // Q6MagLogical->SetVisAttributes(VacVisAtt);
  //  new G4PVPlacement(0, G4ThreeVector(0, 0, pQ6Pos_z), Q6MagLogical, "Q6Mag", world_log, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Helmholtz Coil 'Physical' Volume
  G4double pHLMZRin = 5.10 * cm;   G4double pHLMZRout = 15.0 * cm;   G4double pHLMZHLZ = 5.0 * cm;
  G4double pHLMZ1Pos_X = 0.0 * cm;   G4double pHLMZ1Pos_Y = 0.0 * cm;   G4double pHLMZ1Pos_Z = -12.6 * cm;
  G4double pHLMZ2Pos_X = 0.0 * cm;   G4double pHLMZ2Pos_Y = 0.0 * cm;   G4double pHLMZ2Pos_Z = 12.6 * cm;

  G4VSolid* HLMZ1Solid = new G4Tubs( "HLMZ1Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );
  G4VSolid* HLMZ2Solid = new G4Tubs( "HLMZ2Tube", pHLMZRin, pHLMZRout, pHLMZHLZ, 0.0, 360.0 * deg );

  G4LogicalVolume* HLMZ1Logical = new G4LogicalVolume(HLMZ1Solid, Cu, "Helmholtz1", 0, 0, 0);
  G4LogicalVolume* HLMZ2Logical = new G4LogicalVolume(HLMZ2Solid, Cu, "Helmholtz2", 0, 0, 0);
  HLMZ1Logical->SetVisAttributes(CuVisAtt);
  HLMZ2Logical->SetVisAttributes(CuVisAtt);

  new G4PVPlacement(0, G4ThreeVector(0, 0, pHLMZ1Pos_Z), HLMZ1Logical, "Helmholtz1", world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0, 0, pHLMZ2Pos_Z), HLMZ2Logical, "Helmholtz2", world_log, 0, 0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Beam pipe
  G4double pBPI1Rin = 0.0   * cm; G4double pBPI1Rout = 5.08  * cm; G4double pBPI1HL = 100.0 * cm;
  G4double pBPV1Rin = 0.0   * cm; G4double pBPV1Rout = 4.78  * cm; G4double pBPV1HL = 100.0 * cm;
  G4double pBPI2Rin = 0.0   * cm; G4double pBPI2Rout = 5.08  * cm; G4double pBPI2HL = 12.15 * cm;
  G4double pBPV2Rin = 0.0   * cm; G4double pBPV2Rout = 4.78  * cm; G4double pBPV2HL = 12.15 * cm;
  G4double pBPI3Rin = 0.0   * cm; G4double pBPI3Rout = 3.177 * cm; G4double pBPI3HL = 158.5 * cm;
  G4double pBPV3Rin = 0.0   * cm; G4double pBPV3Rout = 2.877 * cm; G4double pBPV3HL = 158.5 * cm;
  G4double pBPF3Rin = 1.905 * cm; G4double pBPF3Rout = 4.78  * cm; G4double pBPF3HL = 1.5   * cm;
  G4double pBPF4Rin = 3.000 * cm; G4double pBPF4Rout = 4.78  * cm; G4double pBPF4HL = 1.5   * cm;
  G4double pBPI4Rin = 0.0   * cm; G4double pBPI4Rout = 1.905 * cm; G4double pBPI4HL = 200.0 * cm;
  G4double pBPV4Rin = 0.0   * cm; G4double pBPV4Rout = 1.800 * cm; G4double pBPV4HL = 20.0  * cm;
  G4double pBPE4Rin = 2.2   * cm; G4double pBPE4Rout = 8.0   * cm; G4double pBPE4HL = 20.0  * cm;
  G4double pBPI5Rin = 0.0   * cm; G4double pBPI5Rout = 3.175 * cm; G4double pBPI5HL = 22.5  * cm;
  G4double pBPV5Rin = 0.0   * cm; G4double pBPV5Rout = 3.000 * cm; G4double pBPV5HL = 22.5  * cm;

  G4double pBPI1Pos_X = 0.0 * cm; G4double pBPI1Pos_Y = 0.0 * cm; G4double pBPI1Pos_Z = 200.0  * cm;
  G4double pBPV1Pos_X = 0.0 * cm; G4double pBPV1Pos_Y = 0.0 * cm; G4double pBPV1Pos_Z = 0.0    * cm;
  G4double pBPI2Pos_X = 0.0 * cm; G4double pBPI2Pos_Y = 0.0 * cm; G4double pBPI2Pos_Z = 312.15 * cm;
  G4double pBPV2Pos_X = 0.0 * cm; G4double pBPV2Pos_Y = 0.0 * cm; G4double pBPV2Pos_Z = 0.0    * cm;
  G4double pBPI3Pos_X = 0.0 * cm; G4double pBPI3Pos_Y = 0.0 * cm; G4double pBPI3Pos_Z = 726.08 * cm;
  G4double pBPV3Pos_X = 0.0 * cm; G4double pBPV3Pos_Y = 0.0 * cm; G4double pBPV3Pos_Z = 0.0    * cm;
  G4double pBPF3Pos_X = 0.0 * cm; G4double pBPF3Pos_Y = 0.0 * cm; G4double pBPF3Pos_Z = 157.0  * cm;
  G4double pBPF4Pos_X = 0.0 * cm; G4double pBPF4Pos_Y = 0.0 * cm; G4double pBPF4Pos_Z = -157.0 * cm;
  G4double pBPI4Pos_X = 0.0 * cm; G4double pBPI4Pos_Y = 0.0 * cm; G4double pBPI4Pos_Z = 1084.58 * cm;
  G4double pBPV4Pos_X = 0.0 * cm; G4double pBPV4Pos_Y = 0.0 * cm; G4double pBPV4Pos_Z = 0.0    * cm;
  G4double pBPE4Pos_X = 0.0 * cm; G4double pBPE4Pos_Y = 0.0 * cm; G4double pBPE4Pos_Z = 940.0  * cm;
  G4double pBPI5Pos_X = 0.0 * cm; G4double pBPI5Pos_Y = 0.0 * cm; G4double pBPI5Pos_Z = 545.08 * cm;
  G4double pBPV5Pos_X = 0.0 * cm; G4double pBPV5Pos_Y = 0.0 * cm; G4double pBPV5Pos_Z = 0.0    * cm;

  G4VSolid* BPI1Solid = new G4Tubs( "BPI1Tube", pBPI1Rin, pBPI1Rout, pBPI1HL, 0.0, 360.0 * deg );
  G4VSolid* BPV1Solid = new G4Tubs( "BPV1Tube", pBPV1Rin, pBPV1Rout, pBPV1HL, 0.0, 360.0 * deg );
  G4VSolid* BPI2Solid = new G4Tubs( "BPI2Tube", pBPI2Rin, pBPI2Rout, pBPI2HL, 0.0, 360.0 * deg );
  G4VSolid* BPV2Solid = new G4Tubs( "BPV2Tube", pBPV2Rin, pBPV2Rout, pBPV2HL, 0.0, 360.0 * deg );
  G4VSolid* BPI3Solid = new G4Tubs( "BPI3Tube", pBPI3Rin, pBPI3Rout, pBPI3HL, 0.0, 360.0 * deg );
  G4VSolid* BPV3Solid = new G4Tubs( "BPV3Tube", pBPV3Rin, pBPV3Rout, pBPV3HL, 0.0, 360.0 * deg );
  G4VSolid* BPF3Solid = new G4Tubs( "BPF3Tube", pBPF3Rin, pBPF3Rout, pBPF3HL, 0.0, 360.0 * deg );
  G4VSolid* BPF4Solid = new G4Tubs( "BPF4Tube", pBPF4Rin, pBPF4Rout, pBPF4HL, 0.0, 360.0 * deg );
  G4VSolid* BPI4Solid = new G4Tubs( "BPI4Tube", pBPI4Rin, pBPI4Rout, pBPI4HL, 0.0, 360.0 * deg );
  G4VSolid* BPV4Solid = new G4Tubs( "BPV4Tube", pBPV4Rin, pBPV4Rout, pBPV4HL, 0.0, 360.0 * deg );
  G4VSolid* BPE4Solid = new G4Tubs( "BPE4Tube", pBPE4Rin, pBPE4Rout, pBPE4HL, 0.0, 360.0 * deg );
  G4VSolid* BPI5Solid = new G4Tubs( "BPI5Tube", pBPI5Rin, pBPI5Rout, pBPI5HL, 0.0, 360.0 * deg );
  G4VSolid* BPV5Solid = new G4Tubs( "BPV5Tube", pBPV5Rin, pBPV5Rout, pBPV5HL, 0.0, 360.0 * deg );

  G4SubtractionSolid* sub12 = new G4SubtractionSolid("sub12"  , BPI1Solid, BPV1Solid, 0, G4ThreeVector(pBPV1Pos_X, pBPV1Pos_Y, pBPV1Pos_Z) );
  G4SubtractionSolid* sub13 = new G4SubtractionSolid("sub13"  , BPI2Solid, BPV2Solid, 0, G4ThreeVector(pBPV2Pos_X, pBPV2Pos_Y, pBPV2Pos_Z) );
  G4SubtractionSolid* sub14 = new G4SubtractionSolid("sub14"  , BPI3Solid, BPV3Solid, 0, G4ThreeVector(pBPV3Pos_X, pBPV3Pos_Y, pBPV3Pos_Z) );
  G4SubtractionSolid* sub15 = new G4SubtractionSolid("sub15"  , BPI4Solid, BPV4Solid, 0, G4ThreeVector(pBPV4Pos_X, pBPV4Pos_Y, pBPV4Pos_Z) );
  G4SubtractionSolid* sub16 = new G4SubtractionSolid("sub16"  , BPI5Solid, BPV5Solid, 0, G4ThreeVector(pBPV5Pos_X, pBPV5Pos_Y, pBPV5Pos_Z) );

  G4LogicalVolume* sub12Logical = new G4LogicalVolume( sub12, aluminum, "sub12Logical", 0, 0, 0);
  G4LogicalVolume* sub13Logical = new G4LogicalVolume( sub13, aluminum, "sub13Logical", 0, 0, 0);
  G4LogicalVolume* sub14Logical = new G4LogicalVolume( sub14, aluminum, "sub14Logical", 0, 0, 0);
  G4LogicalVolume* sub15Logical = new G4LogicalVolume( sub15, aluminum, "sub15Logical", 0, 0, 0);
  G4LogicalVolume* sub16Logical = new G4LogicalVolume( sub16, aluminum, "sub16Logical", 0, 0, 0);
  //G4LogicalVolume* sub17Logical = new G4LogicalVolume( sub17, siliconsteel, "sub17Logical", 0, 0, 0);

  G4LogicalVolume* BPV1Logical  = new G4LogicalVolume( BPV1Solid, Vacuum, "BPV1SolidLogical", 0, 0, 0);
  G4LogicalVolume* BPV2Logical  = new G4LogicalVolume( BPV2Solid, Vacuum, "BPV2SolidLogical", 0, 0, 0);
  G4LogicalVolume* BPV3Logical  = new G4LogicalVolume( BPV3Solid, Vacuum, "BPV3SolidLogical", 0, 0, 0);
  G4LogicalVolume* BPV4Logical  = new G4LogicalVolume( BPV4Solid, Vacuum, "BPV4SolidLogical", 0, 0, 0);
  G4LogicalVolume* BPV5Logical  = new G4LogicalVolume( BPV5Solid, Vacuum, "BPV5SolidLogical", 0, 0, 0);

  sub12Logical->SetVisAttributes(AlumVisAtt);
  sub13Logical->SetVisAttributes(AlumVisAtt);
  sub14Logical->SetVisAttributes(AlumVisAtt);
  sub15Logical->SetVisAttributes(AlumVisAtt);
  sub16Logical->SetVisAttributes(AlumVisAtt);
  //sub17Logical->SetVisAttributes(LeadVisAtt);

  BPV1Logical->SetVisAttributes(VacVisAtt);
  BPV2Logical->SetVisAttributes(VacVisAtt);
  BPV3Logical->SetVisAttributes(VacVisAtt);
  BPV4Logical->SetVisAttributes(VacVisAtt);
  BPV5Logical->SetVisAttributes(VacVisAtt);

  //comment out all for test
  /*
  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI1Pos_Z), sub12Logical,"BPI1Phys",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI2Pos_Z), sub13Logical,"BPI2Phys",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI3Pos_Z), sub14Logical,"BPI3Phys",world_log, 0, 0, fCheckOverlaps);
  //  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI4Pos_Z), sub15Logical,"BPI4Phys",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI5Pos_Z), sub16Logical,"BPI5Phys",world_log, 0, 0, fCheckOverlaps);
  //  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI3Pos_Z), sub17Logical,"BPI7Phys",world_log, 0, 0, fCheckOverlaps);

  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI1Pos_Z), BPV1Logical, "BPV1Phys",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI2Pos_Z), BPV2Logical, "BPV2Phys",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI3Pos_Z), BPV3Logical, "BPV3Phys",world_log, 0, 0, fCheckOverlaps);
  //  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI4Pos_Z), BPV4Logical, "BPV4Phys",world_log, 0, 0, fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pBPI5Pos_Z), BPV5Logical, "BPV5Phys",world_log, 0, 0, fCheckOverlaps);
  */


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Dipole magnetic field volume
  G4double pDMagHLX   =  8.0   * cm;  G4double pDMagHLY   = 30.0 * cm;  G4double pDMagHLZ   = 82.25 * cm;
  G4double pDMagPos_X =  0.0   * cm;  G4double pDMagPos_Y =  0.0 * cm;  G4double pDMagPos_Z =  423.4 * cm;

  G4VSolid* DMagSolid = new G4Box ( "DMagBox" , pDMagHLX , pDMagHLY , pDMagHLZ );
  G4LogicalVolume* DLogical = new G4LogicalVolume ( DMagSolid, Vacuum, "DipoleMag", 0, 0, 0);
  // DLogical->SetFieldManager(DFieldManager,allLocal);
  // DLogical->SetVisAttributes(VacVisAtt);

  // new G4PVPlacement(0,G4ThreeVector(pDMagPos_X, pDMagPos_Y - 9*cm, pDMagPos_Z), DLogical,"DipoleMag",world_log,0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DIPOLE BOX
  G4double pDBI1HLX =  6.0   * cm;  G4double pDBI1HLY = 16.5 * cm;  G4double pDBI1HLZ = 98.5  * cm;
  G4double pDBV1HLX =  5.295 * cm;  G4double pDBV1HLY = 16.0 * cm;  G4double pDBV1HLZ = 98.5  * cm;
  G4double pDBW1HLX = 11.66  * cm;  G4double pDBW1HLY = 21.4 * cm;  G4double pDBW1HLZ =  0.64 * cm;
  G4double pDBW2Rin =  0.0   * cm;  G4double pDBW2Rout=  3.0 * cm;  G4double pDBW2HL  =  0.64 * cm;
  G4double pDBW3HLX =  1.18  * cm;  G4double pDBW3HLY =  8.0 * cm;  G4double pDBW3HLZ =  0.64 * cm;
  G4double pDBW4HLX =  1.18  * cm;  G4double pDBW4HLY =  8.0 * cm;  G4double pDBW4HLZ =  0.64 * cm;
  G4double pDBT3HLX =  1.18  * cm;  G4double pDBT3HLY =  8.0 * cm;  G4double pDBT3HLZ =  0.005* cm;
  G4double pDBT4HLX =  1.18  * cm;  G4double pDBT4HLY =  8.0 * cm;  G4double pDBT4HLZ =  0.005* cm;
  G4double pDMS1HLX =  3.0   * cm;  G4double pDMS1HLY = 15.0 * cm;  G4double pDMS1HLZ = 96.2  * cm;
  G4double pDMBHRin =  0.0   * cm;  G4double pDMBHRout=  2.0 * cm;  G4double pDMBHHL  = 96.2  * cm;
  G4double pDMFRRin =  1.27  * cm;  G4double pDMFRRout=  2.0 * cm;  G4double pDMFRHL  =  9.0  * cm;

  G4double pDBI1Pos_X =  0.0   * cm;  G4double pDBI1Pos_Y = -9.0 * cm;  G4double pDBI1Pos_Z =422.8  * cm;
  G4double pDBV1Pos_X =  0.0   * cm;  G4double pDBV1Pos_Y =  0.0 * cm;  G4double pDBV1Pos_Z =  0.0  * cm;
  G4double pDBW1Pos_X =  0.0   * cm;  G4double pDBW1Pos_Y = -9.0 * cm;  G4double pDBW1Pos_Z =521.94 * cm;
  G4double pDBW2Pos_X =  0.0   * cm;  G4double pDBW2Pos_Y =  9.0 * cm;  G4double pDBW2Pos_Z =  0.0  * cm;
  G4double pDBW3Pos_X = -4.13  * cm;  G4double pDBW3Pos_Y = -7.0 * cm;  G4double pDBW3Pos_Z =  0.0  * cm;
  G4double pDBW4Pos_X =  4.13  * cm;  G4double pDBW4Pos_Y = -7.0 * cm;  G4double pDBW4Pos_Z =  0.0  * cm;
  G4double pDBT3Pos_X =  0.0   * cm;  G4double pDBT3Pos_Y =  0.0 * cm;  G4double pDBT3Pos_Z =  0.0  * cm;
  G4double pDBT4Pos_X =  0.0   * cm;  G4double pDBT4Pos_Y =  0.0 * cm;  G4double pDBT4Pos_Z =  0.0  * cm;
  G4double pDMS1Pos_X =  0.0   * cm;  G4double pDMS1Pos_Y =  0.0 * cm;  G4double pDMS1Pos_Z =  0.0  * cm;
  G4double pDMBHPos_X =  0.0   * cm;  G4double pDMBHPos_Y =  9.0 * cm;  G4double pDMBHPos_Z =  0.0  * cm;
  G4double pDMFRPos_X =  0.0   * cm;  G4double pDMFRPos_Y =  0.0 * cm;  G4double pDMFRPos_Z =-87.2  * cm;

  G4VSolid* DBI1Solid = new G4Box ( "DBI1Box" , pDBI1HLX , pDBI1HLY  , pDBI1HLZ );
  G4VSolid* DBV1Solid = new G4Box ( "DBV1Box" , pDBV1HLX , pDBV1HLY  , pDBV1HLZ );
  G4VSolid* DBW1Solid = new G4Box ( "DBW1Box" , pDBW1HLX , pDBW1HLY  , pDBW1HLZ );
  G4VSolid* DBW2Solid = new G4Tubs( "DBW2Tubs", pDBW2Rin , pDBW2Rout , pDBW2HL  , 0.0, 360.0 * deg );
  G4VSolid* DBW3Solid = new G4Box ( "DBW3Box" , pDBW3HLX , pDBW3HLY  , pDBW3HLZ );
  G4VSolid* DBW4Solid = new G4Box ( "DBW4Box" , pDBW4HLX , pDBW4HLY  , pDBW4HLZ );
  G4VSolid* DBT3Solid = new G4Box ( "DBT3Box" , pDBT3HLX , pDBT3HLY  , pDBT3HLZ );
  G4VSolid* DBT4Solid = new G4Box ( "DBT4Box" , pDBT4HLX , pDBT4HLY  , pDBT4HLZ );
  G4VSolid* DMS1Solid = new G4Box ( "DMS1Box" , pDMS1HLX , pDMS1HLY  , pDMS1HLZ );
  G4VSolid* DMBHSolid = new G4Tubs( "DMBHTubs", pDMBHRin , pDMBHRout , pDMBHHL  , 0.0, 360.0 * deg );
  G4VSolid* DMFRSolid = new G4Tubs( "DMFRTubs", pDMFRRin , pDMFRRout , pDMFRHL  , 0.0, 360.0 * deg );

  G4SubtractionSolid* sub1 = new G4SubtractionSolid("sub1", DBI1Solid, DBV1Solid, 0,
                           G4ThreeVector(pDBV1Pos_X, pDBV1Pos_Y, pDBV1Pos_Z) );

  G4UnionSolid*       sub2 = new G4UnionSolid("sub2", sub1 , DMS1Solid, 0,
                           G4ThreeVector(pDMS1Pos_X + pDBV1Pos_X,
                           pDMS1Pos_Y + pDBV1Pos_Y, pDMS1Pos_Z + pDBV1Pos_Z) );

  G4SubtractionSolid* sub3 = new G4SubtractionSolid("sub3", sub2 , DMBHSolid, 0,
                           G4ThreeVector(pDMBHPos_X + pDMS1Pos_X + pDBV1Pos_X,
                           pDMBHPos_Y + pDMS1Pos_Y + pDBV1Pos_Y,
                           pDMBHPos_Z + pDMS1Pos_Z + pDBV1Pos_Z) );

  G4UnionSolid*       sub4 = new G4UnionSolid      ("sub4", sub3 , DMFRSolid, 0,
                           G4ThreeVector(pDMFRPos_X + pDMBHPos_X + pDMS1Pos_X + pDBV1Pos_X,
                           pDMFRPos_Y + pDMBHPos_Y + pDMS1Pos_Y + pDBV1Pos_Y,
                           pDMFRPos_Z + pDMBHPos_Z + pDMS1Pos_Z + pDBV1Pos_Z) );

  G4LogicalVolume* sub4Logical = new G4LogicalVolume ( sub4, siliconsteel, "DipoleVacuumBox", 0, 0, 0);
  //  sub4Logical->SetVisAttributes(SteelVisAtt);

  //  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),sub4Logical,"DipoleVacuumBox",DLogical,0,0,fCheckOverlaps);

  G4SubtractionSolid* sub9 = new G4SubtractionSolid("sub9"  , DBW1Solid, DBW2Solid, 0, G4ThreeVector(pDBW2Pos_X, pDBW2Pos_Y, pDBW2Pos_Z) );
  G4SubtractionSolid* sub10= new G4SubtractionSolid("sub10" , sub9     , DBW3Solid, 0, G4ThreeVector(pDBW3Pos_X, pDBW3Pos_Y, pDBW3Pos_Z) );
  G4SubtractionSolid* sub11= new G4SubtractionSolid("sub11" , sub10    , DBW4Solid, 0, G4ThreeVector(pDBW4Pos_X, pDBW4Pos_Y, pDBW4Pos_Z) );

  G4LogicalVolume* sub11Logical = new G4LogicalVolume ( sub11, siliconsteel, "sub11Logical", 0, 0, 0);
  //  sub11Logical->SetVisAttributes(LeadVisAtt);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // COLLIMATOR
  G4double pDCOLHLX = 1.00   * cm;  G4double pDCOLHLY = 4.00 * cm;  G4double pDCOLHLZ = 3.00  * cm;
  G4double pDSLOHLX = 1.00   * cm;  G4double pDSLOHLY = 1.50 * cm;  G4double pDSLOHLZ = 3.00  * cm;
  G4double pDCOLPos_X =-4.00   * cm;  G4double pDCOLPos_Y = 9.00 * cm;  G4double pDCOLPos_Z =-92.00  * cm;
  G4double pDSLOPos_X = 0.00   * cm;  G4double pDSLOPos_Y = 0.00 * cm;  G4double pDSLOPos_Z =  0.00  * cm;

  G4VSolid* DCOLSolid = new G4Box ( "DCOLBox" , pDCOLHLX , pDCOLHLY  , pDCOLHLZ );
  G4VSolid* DSLOSolid = new G4Box ( "DSLOBox" , pDSLOHLX , pDSLOHLY  , pDSLOHLZ );

  G4SubtractionSolid* subcol = new G4SubtractionSolid("subcol", DCOLSolid, DSLOSolid, 0, G4ThreeVector(pDSLOPos_X,pDSLOPos_Y,pDSLOPos_Z) );

  G4LogicalVolume* subcolLogical = new G4LogicalVolume ( subcol, lead, "Collimator", 0, 0, 0);
  subcolLogical ->SetVisAttributes(LeadVisAtt);

  // new G4PVPlacement(0,G4ThreeVector( pDCOLPos_X + pDBV1Pos_X + pDBI1Pos_X,
  // pDCOLPos_Y + pDBV1Pos_Y + pDBI1Pos_Y,
  //         pDCOLPos_Z + pDBV1Pos_Z + pDBI1Pos_Z),subcolLogical,"Collimator1",world_log,0,0,fCheckOverlaps);
 // new G4PVPlacement(0,G4ThreeVector(-pDCOLPos_X + pDBV1Pos_X + pDBI1Pos_X,
  //        pDCOLPos_Y + pDBV1Pos_Y + pDBI1Pos_Y,
  //        pDCOLPos_Z + pDBV1Pos_Z + pDBI1Pos_Z),subcolLogical,"Collimator2",world_log,0,0,fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Dipole Magnets Physical
  G4double pDHLX   = 10.00 * cm;  G4double pDHLY   = 40.00 * cm;  G4double pDHLZ   = 76.00 * cm;
  G4double pD1Pos_X=-16.00 * cm;  G4double pD2Pos_X= 16.00 * cm;  G4double pDPos_Y = -9.00 * cm;  G4double pDPos_Z =422.80 * cm;
  G4VSolid* DSolid  = new G4Box ( "DBox"  , pDHLX , pDHLY  , pDHLZ );
  G4LogicalVolume* DMagLogical  = new G4LogicalVolume(DSolid ,siliconsteel, "DMagLogical" ,0,0,0);
  // DMagLogical->SetVisAttributes(DipVisAtt);

  //  new G4PVPlacement(0,G4ThreeVector(pD1Pos_X,pDPos_Y,pDPos_Z),DMagLogical,"Dipole1",world_log,0,0,0);
  // new G4PVPlacement(0,G4ThreeVector(pD2Pos_X,pDPos_Y,pDPos_Z),DMagLogical,"Dipole2",world_log,0,0,0);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // MAGNETS
  G4double pQ1Rin  =  4.7625 * cm;  G4double pQ1Rout = 20.00 * cm;  G4double pQ1HL   = 18.29 * cm;  G4double pQ1Pos_Z= 85.0 * cm;
  G4double pQ2Rin  =  12.7 * cm;  G4double pQ2Rout = 86.36 * cm;  G4double pQ2HL   = 61.595 * cm;  G4double pQ2Pos_Z= 297.8 * cm;
  G4double pQ3Rin  =  12.7 * cm;  G4double pQ3Rout = 86.36 * cm;  G4double pQ3HL   = 61.595 * cm;  G4double pQ3Pos_Z= 431.1 * cm;

  G4VSolid* Q1Solid = new G4Tubs( "Q1Tubs", pQ1Rin, pQ1Rout, pQ1HL, 0.0, 360.0 * deg );
  G4VSolid* Q2Solid = new G4Tubs( "Q2Tubs", pQ2Rin, pQ2Rout, pQ2HL, 0.0, 360.0 * deg );
  G4VSolid* Q3Solid = new G4Tubs( "Q3Tubs", pQ3Rin, pQ3Rout, pQ3HL, 0.0, 360.0 * deg );
  G4VSolid* Q1MagSolid = new G4Tubs( "Q1MagTubs", 0., pQ1Rin, pQ1HL, 0.0, 360.0 * deg );
  G4VSolid* Q2MagSolid = new G4Tubs( "Q2MagTubs", 0., pQ2Rin, pQ2HL, 0.0, 360.0 * deg );
  G4VSolid* Q3MagSolid = new G4Tubs( "Q3MagTubs", 0., pQ3Rin, pQ3HL, 0.0, 360.0 * deg );

  G4LogicalVolume* Q1Logical = new G4LogicalVolume(Q1Solid,siliconsteel,"Q1Logical",0,0,0);
  G4LogicalVolume* Q2Logical = new G4LogicalVolume(Q2Solid,siliconsteel,"Q2Logical",0,0,0);
  G4LogicalVolume* Q3Logical = new G4LogicalVolume(Q3Solid,siliconsteel,"Q3Logical",0,0,0);

  G4LogicalVolume* Q1MagLogical = new G4LogicalVolume(Q1MagSolid,Vacuum,"Q1MagLogical",0,0,0);
  G4LogicalVolume* Q2MagLogical = new G4LogicalVolume(Q2MagSolid,Vacuum,"Q2MagLogical",0,0,0);
  G4LogicalVolume* Q3MagLogical = new G4LogicalVolume(Q3MagSolid,Vacuum,"Q3MagLogical",0,0,0);

  /*G4VisAttributes *Q1VisAtt(SteelVisAtt);
  Q1VisAtt->SetForceWireframe(true);
  Q1Logical->SetVisAttributes(Q1VisAtt);*/

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

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // Planes for Virtual Detectors
  G4VSolid*        VBSolid   = new G4Tubs("VBSolid",0,pQ1Rout, 0.00001 * mm ,0.0,360.0*deg);
  G4LogicalVolume* Q1ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q1ENLogical",0,0,0);
  G4LogicalVolume* Q1EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q1EXLogical",0,0,0);
  G4LogicalVolume* Q2ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q2ENLogical",0,0,0);
  G4LogicalVolume* Q2EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q2EXLogical",0,0,0);
  G4LogicalVolume* Q3ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q3ENLogical",0,0,0);
  G4LogicalVolume* Q3EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q3EXLogical",0,0,0);
  G4LogicalVolume* Q4ENLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q4ENLogical",0,0,0);
  G4LogicalVolume* Q4EXLogical = new G4LogicalVolume(VBSolid, Vacuum, "Q4EXLogical",0,0,0);

  /*G4VisAttributes *Q1ENVisAtt(VacVisAtt);
  Q1ENVisAtt->SetForceWireframe(true);
  Q1ENLogical->SetVisAttributes(Q1ENVisAtt);*/

  Q1ENLogical->SetVisAttributes(VacVisAtt);
  Q1EXLogical->SetVisAttributes(VacVisAtt);
  Q2ENLogical->SetVisAttributes(VacVisAtt);
  Q2EXLogical->SetVisAttributes(VacVisAtt);
  Q3ENLogical->SetVisAttributes(VacVisAtt);
  Q3EXLogical->SetVisAttributes(VacVisAtt);
  Q4ENLogical->SetVisAttributes(VacVisAtt);
  Q4EXLogical->SetVisAttributes(VacVisAtt);

  MolPolDetector* Q1ENSD = new MolPolDetector("q1en", 1);
  MolPolDetector* Q1EXSD = new MolPolDetector("q1ex", 2);
  MolPolDetector* Q2ENSD = new MolPolDetector("q2en", 3);
  MolPolDetector* Q2EXSD = new MolPolDetector("q2ex", 4);
  MolPolDetector* Q3ENSD = new MolPolDetector("q3en", 5);
  MolPolDetector* Q3EXSD = new MolPolDetector("q3ex", 6);
  MolPolDetector* Q4ENSD = new MolPolDetector("q4en", 7);
  MolPolDetector* Q4EXSD = new MolPolDetector("q4ex", 8);

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
  SDman->AddNewDetector(Q4ENSD);
  SDman->AddNewDetector(Q4EXSD);
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
  Q4ENLogical->SetSensitiveDetector(Q4ENSD);
  Q4EXLogical->SetSensitiveDetector(Q4EXSD);

  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z - pQ1HL ), Q1ENLogical,"VP.Q1.Entr",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ1Pos_Z + pQ1HL ), Q1EXLogical,"VP.Q1.Exit",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z - pQ2HL ), Q2ENLogical,"VP.Q2.Entr",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ2Pos_Z + pQ2HL ), Q2EXLogical,"VP.Q2.Exit",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z - pQ3HL ), Q3ENLogical,"VP.Q3.Entr",world_log,0,0,fCheckOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0,0,pQ3Pos_Z + pQ3HL ), Q3EXLogical,"VP.Q3.Exit",world_log,0,0,fCheckOverlaps);
  //  new G4PVPlacement(0,G4ThreeVector(0,0,pQ4Pos_Z - pQ4HL ), Q4ENLogical,"VP.Q4.Entr",world_log,0,0,fCheckOverlaps);
  //  new G4PVPlacement(0,G4ThreeVector(0,0,pQ4Pos_Z + pQ4HL ), Q4EXLogical,"VP.Q4.Exit",world_log,0,0,fCheckOverlaps);

  // Virtual Plane immediately after target :: these six lines can safely be removed when no longer needed. -Eric King
  G4LogicalVolume* TargVPLogical = new G4LogicalVolume(VBSolid, Vacuum, "TargVPLogical",0,0,0);
  TargVPLogical->SetVisAttributes(VacVisAtt);
  MolPolDetector* TARGVP = new MolPolDetector("targ",16);
  SDman->AddNewDetector(TARGVP);
  TargVPLogical->SetSensitiveDetector(TARGVP);
  new G4PVPlacement(0,G4ThreeVector(0,0,pMTATHLZ + 0.1 * mm ), TargVPLogical,"VP.Targ.Exit",Q6MagLogical,0,0,fCheckOverlaps);

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DIPOLE Virtual Planes
  G4double pVP2HLX    = 6.00 * cm;   G4double pVP2HLY   = 20.00 * cm;  G4double pVP2HLZ   = 0.1   * cm;
  //Note: Vertical positions adjusted by 9*cm.
  G4double pVP2Pos_X  = 0.00  * cm;  G4double pVP2Pos_Y = -9.00  * cm;  G4double pVP2Pos_Z = (422.8 - 98.5 - 1) * cm;
  G4double pVP3Pos_X  = 0.00  * cm;  G4double pVP3Pos_Y = -9.00  * cm;  G4double pVP3Pos_Z = (422.8 + 98.5 + 1) * cm;

  G4VSolid* VP2Solid  = new G4Box( "VP2BOX",  pVP2HLX, pVP2HLY, pVP2HLZ );
  G4VSolid* VP3Solid  = new G4Box( "VP3BOX",  pVP2HLX, pVP2HLY, pVP2HLZ );

  G4LogicalVolume* VP2Logical = new G4LogicalVolume(VP2Solid, Vacuum, "VP2Logical", 0,0,0);
  G4LogicalVolume* VP3Logical = new G4LogicalVolume(VP3Solid, Vacuum, "VP3Logical", 0,0,0);
  VP2Logical->SetSensitiveDetector( DPIN );
  VP3Logical->SetSensitiveDetector( DPOUT );
  VP2Logical->SetVisAttributes(VacVisAtt);
  VP3Logical->SetVisAttributes(VacVisAtt);

  //  new G4PVPlacement(0,G4ThreeVector(pVP2Pos_X, pVP2Pos_Y, pVP2Pos_Z), VP2Logical, "VP.Dipole.Entr", world_log, 0,0, fCheckOverlaps);
  //  new G4PVPlacement(0,G4ThreeVector(pVP3Pos_X, pVP3Pos_Y, pVP3Pos_Z), VP3Logical, "VP.Dipole.Exit", world_log, 0,0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DETECTOR Virtual Plane
  G4double pVP1HLX    = 37.00 * cm;  G4double pVP1HLY   = 50.00 * cm;  G4double pVP1HLZ   = 0.1   * cm;
  G4double pVP1Pos_X  = 0.00  * cm;  G4double pVP1Pos_Y = -46.86* cm;  G4double pVP1Pos_Z = 660.0 * cm;

  G4VSolid* VP1Solid  = new G4Box( "VP1BOX",  pVP1HLX, pVP1HLY, pVP1HLZ );
  G4LogicalVolume* VP1Logical = new G4LogicalVolume(VP1Solid, Vacuum, "VPDetector", 0,0,0);
  VP1Logical->SetSensitiveDetector( DETVP );
  VP1Logical->SetVisAttributes(VacVisAtt);

  //  new G4PVPlacement(0,G4ThreeVector(pVP1Pos_X, pVP1Pos_Y, pVP1Pos_Z), VP1Logical, "VP.Detector.Entr", world_log, 0,0, fCheckOverlaps);


   //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // DETECTOR AND BOX
  G4double pMDBXHLX   = 37.00 * cm;  G4double pMDBXHLY   = 47.10 * cm;  G4double pMDBXHLZ   = 61.20 * cm;
  G4double pMDBAHLX   = 18.50 * cm;  G4double pMDBAHLY   = 29.80 * cm;  G4double pMDBAHLZ   = 40.00 * cm;
  G4double pMDBWHLX   = 13.70 * cm;  G4double pMDBWHLY   = 11.50 * cm;  G4double pMDBWHLZ   =  6.45 * cm;
  G4double pMDBWalpha =  6.0  *deg;  G4double pMDBWtheta = 0.0  *deg;  G4double pMDBWphi   =  0.0  *deg;
  G4double pMDBLHLX   = 13.70 * cm;  G4double pMDBLHLY   = 11.50 * cm;  G4double pMDBLHLZ   =  6.45 * cm;
  G4double pMDBLalpha = 10.0  *deg;  G4double pMDBLtheta = 0.0  *deg;  G4double pMDBLphi   =  0.0  *deg;
  G4double pMDETHLX   =  9.20 * cm;  G4double pMDETHLY   = 17.00 * cm;  G4double pMDETHLZ   = 34.00 * cm;
  G4double pDLGBHLX   =  8.20 * cm;  G4double pDLGBHLY   = 16.20 * cm;  G4double pDLGBHLZ   = 20.10 * cm;

  G4double pMDBXPos_X   =  -49.00 * cm;  G4double pMDBXPos_Y   = 0 * cm;  G4double pMDBXPos_Z   =1109.3 * cm; //10cm up
  G4double pMDBAPos_X   =  0 * cm;  G4double pMDBAPos_Y   = 0 * cm;  G4double pMDBAPos_Z   =  0 * cm;
  G4double pMDBWPos_X   =  0 * cm;  G4double pMDBWPos_Y   = 0 * cm;  G4double pMDBWPos_Z   = 0 * cm;
  G4double pMDBLPos_X   =  0 * cm;  G4double pMDBLPos_Y   = 0 * cm;  G4double pMDBLPos_Z   = 0 * cm;
  G4double pMDETPos_X   =  0 * cm;  G4double pMDETPos_Y   =  0.00 * cm;  G4double pMDETPos_Z   =  0.00 * cm;
  G4double pDLGBPos_X   =  0 * cm;  G4double pDLGBPos_Y   =  0.00 * cm;  G4double pDLGBPos_Z   = 0 * cm;

  //Rotations
  //RMATR07  90.   0.  80.  90.  10.  270.     Moller detector
  //RMATR09  90. 270.   0.   0.  90.  180.     I is oppos Y,II along Z,III oppos X

  G4RotationMatrix* pRot9 = new G4RotationMatrix();
  pRot9->rotateY(90.*deg);
  pRot9->rotateZ(90.*deg);

  G4RotationMatrix* pRot7 = new G4RotationMatrix();
  pRot7->rotateX(-7.3*deg);

  G4VSolid* MDBXSolid  = new G4Box ( "MDBXBox"  , pMDBXHLX, pMDBXHLY, pMDBXHLZ );
  G4VSolid* MDBASolid  = new G4Box ( "MDBABox"  , pMDBAHLX, pMDBAHLY, pMDBAHLZ );
  G4VSolid* MDBWSolid  = new G4Para( "MDBWPara" , pMDBWHLX, pMDBWHLY, pMDBWHLZ, pMDBWalpha, pMDBWtheta, pMDBWphi);
  G4VSolid* MDBLSolid  = new G4Para( "MDBLPara" , pMDBLHLX, pMDBLHLY, pMDBLHLZ, pMDBLalpha, pMDBLtheta, pMDBLphi);
  G4VSolid* MDETSolid  = new G4Box ( "MDETBox"  , pMDETHLX, pMDETHLY, pMDETHLZ );
  G4VSolid* DLGBSolid  = new G4Box ( "DLGBBox"  , pDLGBHLX, pDLGBHLY, pDLGBHLZ );

  G4SubtractionSolid* sub5 = new G4SubtractionSolid("sub5", MDBXSolid, MDBASolid, 0    , G4ThreeVector(pMDBAPos_X, pMDBAPos_Y, pMDBAPos_Z) );
  G4SubtractionSolid* sub6 = new G4SubtractionSolid("sub6", sub5     , MDBWSolid, pRot9, G4ThreeVector(pMDBWPos_X, pMDBWPos_Y, pMDBWPos_Z) );
  G4SubtractionSolid* sub7 = new G4SubtractionSolid("sub7", sub6     , MDBLSolid, pRot9, G4ThreeVector(pMDBLPos_X, pMDBLPos_Y, pMDBLPos_Z) );
  G4SubtractionSolid* sub8 = new G4SubtractionSolid("sub8", sub7     , MDETSolid, pRot7, G4ThreeVector(pMDETPos_X, pMDETPos_Y, pMDETPos_Z) );

  //  G4LogicalVolume* DETLogical = new G4LogicalVolume(DLGBSolid, aluminum, "DETLogical",0,0,0);
  //  DETLogical->SetVisAttributes(AlumVisAtt);
  G4LogicalVolume* DETLogical = new G4LogicalVolume(DLGBSolid, Vacuum, "DETLogical",0,0,0);
  DETLogical->SetSensitiveDetector(DETSD);

  G4LogicalVolume* sub8Logical = new G4LogicalVolume ( sub8, Vacuum, "sub8Logical", 0, 0, 0);
  //  G4LogicalVolume* sub8Logical = new G4LogicalVolume ( sub8, siliconsteel, "sub8Logical", 0, 0, 0);
  G4VisAttributes *sub8VisAtt(AlumVisAtt);
  sub8VisAtt->SetForceWireframe(true);
  sub8Logical->SetVisAttributes(sub8VisAtt);
  sub8Logical->SetSensitiveDetector(DETSD2 );

  new G4PVPlacement(0,G4ThreeVector(pMDBXPos_X, pMDBXPos_Y, pMDBXPos_Z),sub8Logical,"sub8",world_log,0,0,fCheckOverlaps);

  //  new G4PVPlacement(pRot7,G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pDLGBPos_X,
  //        pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pDLGBPos_Y,
  //        pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pDLGBPos_Z),
  //      DETLogical,"virtualBoundaryPhys_det",world_log,0,0,fCheckOverlaps);

 //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    // DETECTOR 2 AND BOX
    G4double pMDBXHLX2   = 37.00 * cm;  G4double pMDBXHLY2   = 47.10 * cm;  G4double pMDBXHLZ2   = 61.20 * cm;
    G4double pMDBAHLX2   = 18.50 * cm;  G4double pMDBAHLY2   = 29.80 * cm;  G4double pMDBAHLZ2   = 40.00 * cm;
    G4double pMDBWHLX2   = 13.70 * cm;  G4double pMDBWHLY2   = 11.50 * cm;  G4double pMDBWHLZ2   = 6.45 * cm;
    G4double pMDBWalpha2 = 6.0  *deg;  G4double pMDBWtheta2 = 0.0  *deg;  G4double pMDBWphi2   =  0.0  *deg;
    G4double pMDBLHLX2   = 13.70 * cm;  G4double pMDBLHLY2   = 11.50 * cm;  G4double pMDBLHLZ2   =  6.45 * cm;
    G4double pMDBLalpha2 = 10.0  *deg;  G4double pMDBLtheta2 = 0.0  *deg;  G4double pMDBLphi2   =  0.0  *deg;
    G4double pMDETHLX2   =  9.20 * cm;  G4double pMDETHLY2   = 17.00 * cm;  G4double pMDETHLZ2   = 34.00 * cm;
    G4double pDLGBHLX2   =  8.20 * cm;  G4double pDLGBHLY2   = 16.20 * cm;  G4double pDLGBHLZ2   = 20.10 * cm;

    G4double pMDBXPos_X2   =  49.00 * cm;  G4double pMDBXPos_Y2   = 0 * cm;  G4double pMDBXPos_Z2   =1109.3 * cm; //10cm up
    G4double pMDBAPos_X2  =  0 * cm;  G4double pMDBAPos_Y2   =  0.0 * cm;  G4double pMDBAPos_Z2   =  0 * cm;
    G4double pMDBWPos_X2   = 0 * cm;  G4double pMDBWPos_Y2   =  0 * cm;  G4double pMDBWPos_Z2   = 0 * cm;
    G4double pMDBLPos_X2   = 0 * cm;  G4double pMDBLPos_Y2   =  0 * cm;  G4double pMDBLPos_Z2   = 0 * cm;
    G4double pMDETPos_X2   = 0 * cm;  G4double pMDETPos_Y2   =  0.00 * cm;  G4double pMDETPos_Z2   =  0.00 * cm;
    G4double pDLGBPos_X2   = 0 * cm;  G4double pDLGBPos_Y2   =  0.00 * cm;  G4double pDLGBPos_Z2   = 0 * cm;

    //Rotations
    //RMATR07  90.   0.  80.  90.  10.  270.     Moller detector
    //RMATR09  90. 270.   0.   0.  90.  180.     I is oppos Y,II along Z,III oppos X

    G4RotationMatrix* pRot92 = new G4RotationMatrix();
    pRot92->rotateY(90.*deg);
    pRot92->rotateZ(90.*deg);

    G4RotationMatrix* pRot72 = new G4RotationMatrix();
    pRot72->rotateX(-7.3*deg);

    G4VSolid* MDBXSolid2  = new G4Box ( "MDBXBox2 "  , pMDBXHLX2, pMDBXHLY2, pMDBXHLZ2 );
    G4VSolid* MDBASolid2  = new G4Box ( "MDBABox2"  , pMDBAHLX2, pMDBAHLY2, pMDBAHLZ2 );
    G4VSolid* MDBWSolid2  = new G4Para( "MDBWPara2" , pMDBWHLX2, pMDBWHLY2, pMDBWHLZ2, pMDBWalpha2, pMDBWtheta2, pMDBWphi2);
    G4VSolid* MDBLSolid2  = new G4Para( "MDBLPara2" , pMDBLHLX2, pMDBLHLY2, pMDBLHLZ2, pMDBLalpha2, pMDBLtheta2, pMDBLphi2);
    G4VSolid* MDETSolid2  = new G4Box ( "MDETBox2"  , pMDETHLX2, pMDETHLY2, pMDETHLZ2 );
    G4VSolid* DLGBSolid2  = new G4Box ( "DLGBBox2"  , pDLGBHLX2, pDLGBHLY2, pDLGBHLZ2 );

    G4SubtractionSolid* sub52 = new G4SubtractionSolid("sub52", MDBXSolid2, MDBASolid2, 0    , G4ThreeVector(pMDBAPos_X2, pMDBAPos_Y2, pMDBAPos_Z2) );
    G4SubtractionSolid* sub62 = new G4SubtractionSolid("sub62", sub52     , MDBWSolid2, pRot92, G4ThreeVector(pMDBWPos_X2, pMDBWPos_Y2, pMDBWPos_Z2) );
    G4SubtractionSolid* sub72 = new G4SubtractionSolid("sub72", sub62     , MDBLSolid2, pRot92, G4ThreeVector(pMDBLPos_X2, pMDBLPos_Y2, pMDBLPos_Z2) );
    G4SubtractionSolid* sub82 = new G4SubtractionSolid("sub82", sub72     , MDETSolid2, pRot72, G4ThreeVector(pMDETPos_X2, pMDETPos_Y2, pMDETPos_Z2) );

    //  G4LogicalVolume* DETLogical = new G4LogicalVolume(DLGBSolid, aluminum, "DETLogical",0,0,0);
    //  DETLogical->SetVisAttributes(AlumVisAtt);
    G4LogicalVolume* DETLogical2 = new G4LogicalVolume(DLGBSolid2, Vacuum, "DETLogical2",0,0,0);
    DETLogical2->SetSensitiveDetector(DETSD);

    G4LogicalVolume* sub8Logical2 = new G4LogicalVolume ( sub82, Vacuum, "sub8Logical2", 0, 0, 0);
    //  G4LogicalVolume* sub8Logical = new G4LogicalVolume ( sub8, siliconsteel, "sub8Logical", 0, 0, 0);
    G4VisAttributes *sub8VisAtt2(AlumVisAtt);
    sub8VisAtt2->SetForceWireframe(true);
    sub8Logical2->SetVisAttributes(sub8VisAtt2);
    sub8Logical2->SetSensitiveDetector(DETSD );

    new G4PVPlacement(0,G4ThreeVector(pMDBXPos_X2, pMDBXPos_Y2, pMDBXPos_Z2),sub8Logical2,"sub82",world_log,0,0,fCheckOverlaps);

    // new G4PVPlacement(pRot72,G4ThreeVector(pMDBXPos_X2 + pMDBAPos_X2 + pMDETPos_X2 + pDLGBPos_X2,
    //        pMDBXPos_Y2 + pMDBAPos_Y2 + pMDETPos_Y2 + pDLGBPos_Y2,
    //        pMDBXPos_Z2 + pMDBAPos_Z2 + pMDETPos_Z2 + pDLGBPos_Z2),
    //      DETLogical2,"virtualBoundaryPhys_det2",world_log,0,0,fCheckOverlaps);


 //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    //Pb DETECTOR AND BOX
    G4double pMDBXHLX3   = 20.00 * cm;  G4double pMDBXHLY3   = 14.00 * cm;  G4double pMDBXHLZ3   = 23.00 * cm;
    G4double pMDBXPos_X3   =  49.00 * cm;  G4double pMDBXPos_Y3   = 0 * cm;  G4double pMDBXPos_Z3   = 1109.3 * cm; 


    G4RotationMatrix* pRot93 = new G4RotationMatrix();
    pRot93->rotateY(90.*deg);
    pRot93->rotateZ(90.*deg);

    G4RotationMatrix* pRot73 = new G4RotationMatrix();
    pRot73->rotateX(-7.3*deg);

    G4VSolid* MDBXSolid3  = new G4Box ( "MDBXBox3 "  , pMDBXHLX3, pMDBXHLY3, pMDBXHLZ3 );
    G4LogicalVolume* DETLogical3 = new G4LogicalVolume(MDBXSolid3, LgTF1, "DETLogical3",0,0,0);
    DETLogical3->SetVisAttributes(ScintVisAtt);
    MolPolDetector* MDBXBox3 = new MolPolDetector("PbDetector",20);
    SDman->AddNewDetector(MDBXBox3);
    DETLogical3->SetSensitiveDetector(MDBXBox3);

    new G4PVPlacement(0,G4ThreeVector(pMDBXPos_X3, pMDBXPos_Y3, pMDBXPos_Z3),DETLogical3,"LG",world_log,0,0,fCheckOverlaps);

 //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    //Pb DETECTOR 2 AND BOX
    G4double pMDBXHLX4   = 20.00 * cm;  G4double pMDBXHLY4   = 14.00 * cm;  G4double pMDBXHLZ4   = 23.00 * cm;
    G4double pMDBXPos_X4   = -49.00 * cm;  G4double pMDBXPos_Y4   = 0 * cm;  G4double pMDBXPos_Z4   = 1109.3 * cm; 


    G4RotationMatrix* pRot94 = new G4RotationMatrix();
    pRot94->rotateY(90.*deg);
    pRot94->rotateZ(90.*deg);

    G4RotationMatrix* pRot74 = new G4RotationMatrix();
    pRot74->rotateX(-7.3*deg);

    G4VSolid* MDBXSolid4  = new G4Box ( "MDBXBox4 "  , pMDBXHLX4, pMDBXHLY4, pMDBXHLZ4 );
    G4LogicalVolume* DETLogical4 = new G4LogicalVolume(MDBXSolid4,LgTF1 , "DETLogical4",0,0,0);
    DETLogical4->SetVisAttributes(ScintVisAtt);
    MolPolDetector* MDBXBox4 = new MolPolDetector("PbDetector2",21);
    SDman->AddNewDetector(MDBXBox4);
    DETLogical4->SetSensitiveDetector(MDBXBox4);
    new G4PVPlacement(0,G4ThreeVector(pMDBXPos_X4, pMDBXPos_Y4, pMDBXPos_Z4),DETLogical4,"Pb2",world_log,0,0,fCheckOverlaps);
 
//////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    // HODOSCOPE 1
	
    //X positions of each hodoscope segment 
    G4double HODO1X   = 53.207 * cm;  G4double HODO4X = 51.404 * cm;  G4double HODO7X   = 49.601 * cm;
    G4double HODO2X   = 53.606 * cm;  G4double HODO5X   = 50.803 * cm;  G4double HODO8X   = 48.399 * cm;
    G4double HODO3X   = 52.005 * cm;  G4double HODO6X   = 50.202 * cm;  G4double HODO9X   = 47.798 * cm;
    G4double HODO10X   = 47.197 * cm;  G4double HODO11X   = 46.596 * cm;  G4double HODO12X   =  45.995 * cm;	
    G4double HODO13X   =  44.394 * cm;  G4double HODO14X   = 44.793 * cm;  
  
	
    //Y and Z positions
    G4double HODOY   = 0.0 * cm;
    G4double HODOZ   =  1000.0 * cm;
    G4double HODOSplitNegY = -2.0 * cm;
    G4double HODOSplitPosY = 2.0 * cm;
    
	
    //X, Y and Z Dimmensions
	
    G4double HODOX_DIM = 1.2 * cm;
    G4double HODOY_DIM = 100.0 * cm;
    G4double HODOZ_DIM = 0.8 * cm;
    G4double HODOSplitY = 50.0 * cm;

    //Rotations, not sure how important this is
    //RMATR07  90.   0.  80.  90.  10.  270.     Moller detector
    //RMATR09  90. 270.   0.   0.  90.  180.     I is oppos Y,II along Z,III oppos X

    G4RotationMatrix* pRot911 = new G4RotationMatrix();
    pRot911->rotateY(90.*deg);
    pRot911->rotateZ(90.*deg);

    G4RotationMatrix* pRot711 = new G4RotationMatrix();
    pRot711->rotateX(-7.3*deg);

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
    G4LogicalVolume* HODO1Logical = new G4LogicalVolume(HODOBOX1, Vacuum, "HODO1Logical",0,0,0);
    HODO1Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO1 = new MolPolDetector("HODO1",25);
    SDman -> AddNewDetector(HODO1);
    HODO1Logical -> SetSensitiveDetector(HODO1);

    G4LogicalVolume* HODO2Logical = new G4LogicalVolume(HODOBOX2, Vacuum, "HODO2",0,0,0);
    HODO2Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO2 = new MolPolDetector("HODO2",26);
    SDman -> AddNewDetector(HODO2);
    HODO2Logical -> SetSensitiveDetector(HODO2);

    G4LogicalVolume* HODO3Logical = new G4LogicalVolume(HODOBOX3, Vacuum, "HODO3",0,0,0);
    HODO3Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO3 = new MolPolDetector("HODO3",27);
    SDman -> AddNewDetector(HODO3);
    HODO3Logical -> SetSensitiveDetector(HODO3); 
 
    G4LogicalVolume*HODO3SplitLogical = new G4LogicalVolume(HODOBOX3Split,Vacuum, "HODO3Split",0,0,0);
    HODO3SplitLogical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO3Split = new MolPolDetector("HODO3Split",28);
    SDman -> AddNewDetector(HODO3Split);
    HODO3SplitLogical -> SetSensitiveDetector(HODO3Split);

    G4LogicalVolume* HODO4Logical = new G4LogicalVolume(HODOBOX4, Vacuum, "HODO4",0,0,0);
    HODO4Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO4 = new MolPolDetector("HODO4",29);
    SDman -> AddNewDetector(HODO4);
    HODO4Logical -> SetSensitiveDetector(HODO4);


    G4LogicalVolume* HODO5Logical = new G4LogicalVolume(HODOBOX5, Vacuum, "HODO5",0,0,0);
    HODO5Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO5 = new MolPolDetector("HODO5",30);
    SDman -> AddNewDetector(HODO5);
    HODO5Logical -> SetSensitiveDetector(HODO5);
	
    G4LogicalVolume* HODO6Logical = new G4LogicalVolume(HODOBOX6, Vacuum, "HODO6",0,0,0);
    HODO6Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO6 = new MolPolDetector("HODO6",31);
    SDman -> AddNewDetector(HODO6);
    HODO1Logical -> SetSensitiveDetector(HODO6);
	
    G4LogicalVolume* HODO7Logical = new G4LogicalVolume(HODOBOX7, Vacuum, "HODO7",0,0,0);
    HODO7Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO7 = new MolPolDetector("HODO7",32);
    SDman -> AddNewDetector(HODO7);
    HODO7Logical -> SetSensitiveDetector(HODO7);
	
    G4LogicalVolume* HODO8Logical = new G4LogicalVolume(HODOBOX8, Vacuum, "HODO8",0,0,0);
    HODO8Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO8 = new MolPolDetector("HODO8",33);
    SDman -> AddNewDetector(HODO8);
    HODO8Logical -> SetSensitiveDetector(HODO8);

    G4LogicalVolume* HODO9Logical = new G4LogicalVolume(HODOBOX9, Vacuum, "HODO9",0,0,0);
    HODO9Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO9 = new MolPolDetector("HODO9",34);
    SDman -> AddNewDetector(HODO9);
    HODO9Logical -> SetSensitiveDetector(HODO9);

    G4LogicalVolume* HODO10Logical = new G4LogicalVolume(HODOBOX10, Vacuum, "HODO10",0,0,0);
    HODO10Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO10 = new MolPolDetector("HODO10",35);
    SDman -> AddNewDetector(HODO10);
    HODO10Logical -> SetSensitiveDetector(HODO10);
	
    G4LogicalVolume* HODO11Logical = new G4LogicalVolume(HODOBOX11, Vacuum, "HODO11",0,0,0);
    HODO11Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO11 = new MolPolDetector("HODO11",36);
    SDman -> AddNewDetector(HODO11);
    HODO11Logical -> SetSensitiveDetector(HODO11);

    G4LogicalVolume* HODO12Logical = new G4LogicalVolume(HODOBOX12, Vacuum, "HODO12",0,0,0);
    HODO12Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO12 = new MolPolDetector("HODO12",37);
    SDman -> AddNewDetector(HODO12);
    HODO12Logical -> SetSensitiveDetector(HODO12);

    G4LogicalVolume*HODO12SplitLogical = new G4LogicalVolume(HODOBOX12Split,Vacuum, "HODO12Split",0,0,0);
    HODO12SplitLogical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO12Split = new MolPolDetector("HODO12Split",38);
    SDman -> AddNewDetector(HODO12Split);
    HODO12SplitLogical -> SetSensitiveDetector(HODO12Split);

    G4LogicalVolume* HODO13Logical = new G4LogicalVolume(HODOBOX13, Vacuum, "HODO13",0,0,0);
    HODO13Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO13 = new MolPolDetector("HODO13",39);
    SDman -> AddNewDetector(HODO13);
    HODO13Logical -> SetSensitiveDetector(HODO13);

    G4LogicalVolume* HODO14Logical = new G4LogicalVolume(HODOBOX14, Vacuum, "HODO14",0,0,0);
    HODO14Logical-> SetVisAttributes(VacVisAtt);
    MolPolDetector*HODO14 = new MolPolDetector("HODO14",40);
    SDman -> AddNewDetector(HODO14);
    HODO14Logical -> SetSensitiveDetector(HODO14);				
		
    	
    new G4PVPlacement(0,G4ThreeVector(HODO1X, HODOY, HODOZ),HODO1Logical,"HODO_1",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO2X, HODOY, HODOZ),HODO2Logical,"HODO_2",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO3X, HODOSplitPosY, HODOZ),HODO3Logical,"HODO_3",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO3X, HODOSplitNegY, HODOZ),HODO3SplitLogical,"HODO_3Lower",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO4X, HODOY, HODOZ),HODO4Logical,"HODO_4",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO5X, HODOY, HODOZ),HODO5Logical,"HODO_5",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO6X, HODOY, HODOZ),HODO6Logical,"HODO_6",world_log,0,0,fCheckOverlaps);	
    new G4PVPlacement(0,G4ThreeVector(HODO7X, HODOY, HODOZ),HODO7Logical,"HODO_7",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO8X, HODOY, HODOZ),HODO8Logical,"HODO_8",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO9X, HODOY, HODOZ),HODO9Logical,"HODO_9",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO10X, HODOY, HODOZ),HODO10Logical,"HODO_10",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO11X, HODOY, HODOZ),HODO11Logical,"HODO_11",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO12X, HODOSplitPosY, HODOZ),HODO12Logical,"HODO_12",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO12X, HODOSplitNegY, HODOZ),HODO12SplitLogical,"HODO_12Lower",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO13X, HODOY, HODOZ),HODO13Logical,"HODO_13",world_log,0,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(HODO14X, HODOY, HODOZ),HODO14Logical,"HODO_14",world_log,0,0,fCheckOverlaps);
   

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
    // HODOSCOPE 2
	
 //X positions of each hodoscope segment 
    G4double HODO1X2   = -53.207 * cm;  G4double HODO4X2 = -51.404 * cm;  G4double HODO7X2   = -49.601 * cm;
    G4double HODO2X2   = -53.606 * cm;  G4double HODO5X2   = -50.803 * cm;  G4double HODO8X2   = -48.399 * cm;
    G4double HODO3X2   = -52.005 * cm;  G4double HODO6X2   = -50.202 * cm;  G4double HODO9X2   = -47.798 * cm;
    G4double HODO10X2   = -47.197 * cm;  G4double HODO11X2   = -46.596 * cm;  G4double HODO12X2   =  -45.995 * cm;	
    G4double HODO13X2   =  -44.394 * cm;  G4double HODO14X2   = -44.793 * cm; 
	
	//Y and Z positions
    G4double HODOY2   = 0 * cm;
    G4double HODOZ2   =  1000.0 * cm;
    G4double HODOSplitNegY2 = -2.0 * cm;
    G4double HODOSplitPosY2 = 2.0 * cm;

	
	//X, Y and Z Dimmensions
	
	G4double HODOX_DIM2 = 1.2 * cm;
	G4double HODOY_DIM2 = 100.0 * cm;
	G4double HODOZ_DIM2 = 0.8 * cm;
	G4double HODOSplitY2 = 50.0 * cm;

    //Rotations, not sure how important this is
    //RMATR07  90.   0.  80.  90.  10.  270.     Moller detector
    //RMATR09  90. 270.   0.   0.  90.  180.     I is oppos Y,II along Z,III oppos X

    G4RotationMatrix* pRot9112 = new G4RotationMatrix();
    pRot9112->rotateY(90.*deg);
    pRot9112->rotateZ(90.*deg);

    G4RotationMatrix* pRot7112 = new G4RotationMatrix();
    pRot7112->rotateX(-7.3*deg);

	//Creating the 14 hodoscope boxes
        G4VSolid* HODOBOX012  = new G4Box ( "HODOBOX012 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX22  = new G4Box ( "HODOBOX22 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX32  = new G4Box ( "HODOBOX32 "  , HODOX_DIM2, HODOSplitY2, HODOZ_DIM2 );
	G4VSolid* HODOBOX32Split = new G4Box ( "HODOBOX32Split" , HODOX_DIM2, HODOSplitY2,  HODOZ_DIM2);
	G4VSolid* HODOBOX42  = new G4Box ( "HODOBOX42"  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX52  = new G4Box ( "HODOBOX52 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX62  = new G4Box ( "HODOBOX62 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX72  = new G4Box ( "HODOBOX72 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX82  = new G4Box ( "HODOBOX82 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX92  = new G4Box ( "HODOBOX92 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX102  = new G4Box ( "HODOBOX102"  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX112  = new G4Box ( "HODOBOX112 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX122  = new G4Box ( "HODOBOX122 "  , HODOX_DIM2, HODOSplitY2, HODOZ_DIM2 );
	G4VSolid* HODOBOX122Split = new G4Box ( "HODOBOX122Split" , HODOX_DIM2, HODOSplitY2,  HODOZ_DIM2);
	G4VSolid* HODOBOX132  = new G4Box ( "HODOBOX132 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	G4VSolid* HODOBOX142  = new G4Box ( "HODOBOX142 "  , HODOX_DIM2, HODOY_DIM2, HODOZ_DIM2 );
	
	
	G4LogicalVolume* HODO012Logical = new G4LogicalVolume(HODOBOX012, Vacuum, "HODO012",0,0,0);
	HODO012Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO012 = new MolPolDetector("HODO012",41);
	SDman -> AddNewDetector(HODO012);
	HODO012Logical -> SetSensitiveDetector(HODO012);        

	G4LogicalVolume* HODO22Logical = new G4LogicalVolume(HODOBOX22, Vacuum, "HODO22",0,0,0);
        HODO22Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO22 = new MolPolDetector("HODO22",42);
	SDman -> AddNewDetector(HODO22);
	HODO22Logical -> SetSensitiveDetector(HODO22); 

	G4LogicalVolume* HODO32Logical = new G4LogicalVolume(HODOBOX32, Vacuum, "HODO32",0,0,0);
        HODO32Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO32 = new MolPolDetector("HODO32",43);
	SDman -> AddNewDetector(HODO32);
	HODO32Logical -> SetSensitiveDetector(HODO32);

	G4LogicalVolume* HODO32SplitLogical = new G4LogicalVolume(HODOBOX32Split,Vacuum, "HODO32Split",0,0,0);
	HODO32SplitLogical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO32Split = new MolPolDetector("HODO32Split",44);
	SDman -> AddNewDetector(HODO32Split);
	HODO32SplitLogical -> SetSensitiveDetector(HODO32Split);

	G4LogicalVolume* HODO42Logical = new G4LogicalVolume(HODOBOX42, Vacuum, "HODO42",0,0,0);
	HODO42Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO42 = new MolPolDetector("HODO42",45);
	SDman -> AddNewDetector(HODO42);
	HODO42Logical -> SetSensitiveDetector(HODO42);        


	G4LogicalVolume* HODO52Logical = new G4LogicalVolume(HODOBOX52, Vacuum, "HODO52",0,0,0);
        HODO52Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO52 = new MolPolDetector("HODO52",46);
	SDman -> AddNewDetector(HODO52);
	HODO52Logical -> SetSensitiveDetector(HODO52);
	
	G4LogicalVolume* HODO62Logical = new G4LogicalVolume(HODOBOX62, Vacuum, "HODO62",0,0,0);
	HODO62Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO62 = new MolPolDetector("HODO62",47);
	SDman -> AddNewDetector(HODO62);
	HODO62Logical -> SetSensitiveDetector(HODO62);
	
	G4LogicalVolume* HODO72Logical = new G4LogicalVolume(HODOBOX72, Vacuum, "HODO72",0,0,0);
	HODO72Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO72 = new MolPolDetector("HODO72",48);
	SDman -> AddNewDetector(HODO72);
	HODO72Logical -> SetSensitiveDetector(HODO72);
	
	G4LogicalVolume* HODO82Logical = new G4LogicalVolume(HODOBOX82, Vacuum, "HODO82",0,0,0);
	HODO82Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO82 = new MolPolDetector("HODO82",50);
	SDman -> AddNewDetector(HODO82);
	HODO82Logical -> SetSensitiveDetector(HODO82);

	G4LogicalVolume* HODO92Logical = new G4LogicalVolume(HODOBOX92, Vacuum, "HODO92",0,0,0);
	HODO92Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO92 = new MolPolDetector("HODO92",51);
	SDman -> AddNewDetector(HODO92);
	HODO92Logical -> SetSensitiveDetector(HODO92);

	G4LogicalVolume* HODO102Logical = new G4LogicalVolume(HODOBOX102, Vacuum, "HODO102",0,0,0);
	HODO102Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO102 = new MolPolDetector("HODO102",52);
	SDman -> AddNewDetector(HODO102);
	HODO102Logical -> SetSensitiveDetector(HODO102);
		
	G4LogicalVolume* HODO112Logical = new G4LogicalVolume(HODOBOX112, Vacuum, "HODO112",0,0,0);
	HODO112Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO112 = new MolPolDetector("HODO112",53);
	SDman -> AddNewDetector(HODO112);
	HODO112Logical -> SetSensitiveDetector(HODO112);

	G4LogicalVolume* HODO122Logical = new G4LogicalVolume(HODOBOX122, Vacuum, "HODO122",0,0,0);
        HODO122Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO122 = new MolPolDetector("HODO122",54);
	SDman -> AddNewDetector(HODO122);
	HODO122Logical -> SetSensitiveDetector(HODO122);

	G4LogicalVolume*HODO122SplitLogical = new G4LogicalVolume(HODOBOX122Split,Vacuum, "HODO122Split",0,0,0);
	HODO122SplitLogical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO122Split = new MolPolDetector("HODO122Split",55);
	SDman -> AddNewDetector(HODO122Split);
	HODO122SplitLogical -> SetSensitiveDetector(HODO122Split);

	G4LogicalVolume* HODO132Logical = new G4LogicalVolume(HODOBOX132, Vacuum, "HODO132",0,0,0);
	HODO132Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO132 = new MolPolDetector("HODO132",56);
	SDman -> AddNewDetector(HODO132);
	HODO132Logical -> SetSensitiveDetector(HODO132);

	G4LogicalVolume* HODO142Logical = new G4LogicalVolume(HODOBOX142, Vacuum, "HODO142",0,0,0);
	HODO142Logical-> SetVisAttributes(VacVisAtt);
	MolPolDetector*HODO142 = new MolPolDetector("HODO142",57);
	SDman -> AddNewDetector(HODO142);
	HODO142Logical -> SetSensitiveDetector(HODO142);        
    
		
		
	new G4PVPlacement(0,G4ThreeVector(HODO1X2, HODOY, HODOZ),HODO012Logical,"HODO_012",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO2X2, HODOY, HODOZ),HODO22Logical,"HODO_22",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO3X2, HODOSplitPosY, HODOZ),HODO32Logical,"HODO_32",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO3X2, HODOSplitNegY, HODOZ),HODO32SplitLogical,"HODO_32Lower",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO4X2, HODOY, HODOZ),HODO42Logical,"HODO_42",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO5X2, HODOY, HODOZ),HODO52Logical,"HODO_52",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO6X2, HODOY, HODOZ),HODO62Logical,"HODO_62",world_log,0,0,fCheckOverlaps);	
	new G4PVPlacement(0,G4ThreeVector(HODO7X2, HODOY, HODOZ),HODO72Logical,"HODO_72",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO8X2, HODOY, HODOZ),HODO82Logical,"HODO_82",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO9X2, HODOY, HODOZ),HODO92Logical,"HODO_92",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO10X2, HODOY, HODOZ),HODO102Logical,"HODO_102",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO11X2, HODOY, HODOZ),HODO112Logical,"HODO_112",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO12X2, HODOSplitPosY, HODOZ),HODO122Logical,"HODO_122",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO12X2, HODOSplitNegY, HODOZ),HODO32SplitLogical,"HODO_122Lower",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO13X2, HODOY, HODOZ),HODO132Logical,"HODO_132",world_log,0,0,fCheckOverlaps);
	new G4PVPlacement(0,G4ThreeVector(HODO14X2, HODOY, HODOZ),HODO142Logical,"HODO_142",world_log,0,0,fCheckOverlaps);
	
       



  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // FIXME!!!!! This section of the geometry seems incomplete.
  // Hodoscope
  G4double pHOD1HLX   = 9.2 * cm;   G4double pHOD1HLY   = 15.5 * cm;   G4double pHOD1HLZ   = 0.75 * cm;
  G4double pHOD2HLX   = 9.2 * cm;   G4double pHOD2HLY   = 16.5 * cm;   G4double pHOD2HLZ   = 15.2 * cm;
  G4double pHOD3HLX   = 9.2 * cm;   G4double pHOD3HLY   = 17.5 * cm;   G4double pHOD3HLZ   = 0.1  * cm;

  G4double pHOD1Pos_X = 0.0 * cm;   G4double pHOD1Pos_Y = 0.0  * cm;   G4double pHOD1Pos_Z = -32.5* cm;
  G4double pHOD2Pos_X = 0.0 * cm;   G4double pHOD2Pos_Y = 0.0  * cm;   G4double pHOD2Pos_Z = -16.0* cm;
  G4double pHOD3Pos_X = 0.0 * cm;   G4double pHOD3Pos_Y = 0.0  * cm;   G4double pHOD3Pos_Z = -33.5* cm;

  G4VSolid* HOD1Solid = new G4Box( "HOD1Box", pHOD1HLX, pHOD1HLY, pHOD1HLZ );
  G4VSolid* HOD2Solid = new G4Box( "HOD2Box", pHOD2HLX, pHOD2HLY, pHOD2HLZ );
  G4VSolid* HOD3Solid = new G4Box( "HOD3Box", pHOD3HLX, pHOD3HLY, pHOD3HLZ );

  G4double offset = 3.5 * cm;
  // G4double offset = 0 * cm;

  G4double pAPP1HLX   = 2.0 * cm;   G4double pAPP1HLY   = 15.5 * cm;    G4double pAPP1HLZ   = 0.75 * cm;
  G4double pAPP1Pos_X = 4.4 * cm;   G4double pAPP1Pos_Y = 0.0 * cm + offset;    G4double pAPP1Pos_Z = 0.0 * cm;

  G4VSolid* APP1LSolid = new G4Box( "APP1LBOX", pAPP1HLX, pAPP1HLY, pAPP1HLZ );
  G4VSolid* APP1RSolid = new G4Box( "APP1RBOX", pAPP1HLX, pAPP1HLY, pAPP1HLZ );

  //  G4LogicalVolume* APP1LLogical = new G4LogicalVolume(APP1LSolid, scint, "APP1LLogical", 0,0,0);
  G4LogicalVolume* APP1LLogical = new G4LogicalVolume(APP1LSolid, Vacuum, "APP1LLogical", 0,0,0);
  APP1LLogical->SetVisAttributes(ScintVisAtt);
  //  APP1LLogical->SetSensitiveDetector(APPSD1);

  //  G4LogicalVolume* APP1RLogical = new G4LogicalVolume(APP1RSolid, scint, "APP1RLogical", 0,0,0);
  G4LogicalVolume* APP1RLogical = new G4LogicalVolume(APP1RSolid, Vacuum, "APP1RLogical", 0,0,0);
  APP1RLogical->SetVisAttributes(ScintVisAtt);
  //  APP1RLogical->SetSensitiveDetector(APPSD2);

  // new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X + pAPP1Pos_X,
  //           pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y,
  //         pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pAPP1Pos_Z),
  //    APP1LLogical, "APP1L", world_log, 0,0, fCheckOverlaps);

  // new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X - pAPP1Pos_X,
  //         pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y,
  //         pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pAPP1Pos_Z),
  //    APP1RLogical, "APP1R", world_log, 0,0, fCheckOverlaps);


  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  // HODOSCOPE Virtual Planes
  G4double pVPHODHLZ = 0.1 * cm;
  G4double pVPHODPos_Z = (0.0 - 0.85) * cm;

  G4VSolid* VPHOD1Solid = new G4Box( "VPHOD1", pAPP1HLX, pAPP1HLY, pVPHODHLZ);
  G4VSolid* VPHOD2Solid = new G4Box( "VPHOD2", pAPP1HLX, pAPP1HLY, pVPHODHLZ);

  G4LogicalVolume* VPHOD1Logical = new G4LogicalVolume(VPHOD1Solid, Vacuum, "VPHOD1Logical", 0,0,0);
  VPHOD1Logical->SetSensitiveDetector(APPSD1);
  VPHOD1Logical->SetVisAttributes(VacVisAtt);
  G4LogicalVolume* VPHOD2Logical = new G4LogicalVolume(VPHOD2Solid, Vacuum, "VPHOD2Logical", 0,0,0);
  VPHOD2Logical->SetSensitiveDetector(APPSD2);
  VPHOD2Logical->SetVisAttributes(VacVisAtt);

// new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X + pAPP1Pos_X,
//           pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y,
//           pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pVPHODPos_Z),
//      VPHOD1Logical, "VP.Hodo.Left", world_log, 0,0, fCheckOverlaps);

//  new G4PVPlacement(pRot7, G4ThreeVector(pMDBXPos_X + pMDBAPos_X + pMDETPos_X + pHOD1Pos_X - pAPP1Pos_X,
//           pMDBXPos_Y + pMDBAPos_Y + pMDETPos_Y + pHOD1Pos_Y + pAPP1Pos_Y,
//           pMDBXPos_Z + pMDBAPos_Z + pMDETPos_Z + pHOD1Pos_Z + pVPHODPos_Z),
//      VPHOD2Logical, "VP.Hodo.Right", world_log, 0,0, fCheckOverlaps);


  // // // // // // // // Sheilding // // // // // // // //

  /*
    GPARVOL120 'S1LD'  13  'HALL'    0.  -29.  537.   0  'BOX '  3  30.  25. 14.0
    GPARVOL121 'S1H1'  15  'S1LD'   -4.1   9.5   0.   9  'PARA'  6  10.5 14. 1.6 10. 0. 0.
    GPARVOL122 'S1H2'  15  'S1LD'    4.1   9.5   0.   9  'PARA'  6  10.5 14. 1.6 10. 0. 0.
    GPARVOL123 'S2LD'  13  'HALL'    0.   -7.5 610.   0  'BOX '  3  17.   3.5 50.0 !!!***
    GPARVOL124 'S21B'   9  'HALL'    0.  -13.3 620.   0  'BOX '  3  31.   1.9  3.8 !!!***
    GPARVOL125 'S21H'  15  'S21B'    0.   -0.5   0.   0  'BOX '  3  31.   1.4  3.3
    GPARVOL126 'S3LD'  13  'HALL'    0.   -7.5 723.   0  'BOX '  3  20.   2.0 57.0  !!!***
  */

  G4double pS1LDHLX   = 30.0 * cm;  G4double pS1LDHLY   = 25.0 * cm;  G4double pS1LDHLZ = 14.0 * cm;
  G4double pS1H1HLX   = 10.5 * cm;  G4double pS1H1HLY   = 14.0 * cm;  G4double pS1H1HLZ = 1.6  * cm;
  G4double pS1H1alpha = 10.0 * deg; G4double pS1H1theta = 0.0  * deg; G4double pS1H1phi = 0.0  * deg;
  G4double pS1H2HLX   = 10.5 * cm;  G4double pS1H2HLY   = 14.0 * cm;  G4double pS1H2HLZ = 1.6  * cm;
  G4double pS1H2alpha = 10.0 * deg; G4double pS1H2theta = 0.0  * deg; G4double pS1H2phi = 0.0  * deg;
  G4double pS2LDHLX   = 17.0 * cm;  G4double pS2LDHLY   = 3.5  * cm;  G4double pS2LDHLZ = 50.0 * cm;
  G4double pS21BHLX   = 31.0 * cm;  G4double pS21BHLY   = 1.9  * cm;  G4double pS21BHLZ = 3.8  * cm;
  G4double pS21HHLX   = 31.0 * cm;  G4double pS21HHLY   = 1.4  * cm;  G4double pS21HHLZ = 3.3  * cm;
  G4double pS3LDHLX   = 20.0 * cm;  G4double pS3LDHLY   = 2.0  * cm;  G4double pS3LDHLZ = 57.0 * cm;

  G4double pS1LDPos_X   = 0.0  * cm;  G4double pS1LDPos_Y   = -29.0 * cm;  G4double pS1LDPos_Z = 537.0 * cm;
  G4double pS1H1Pos_X   = -4.1 * cm;  G4double pS1H1Pos_Y   = 9.5   * cm;  G4double pS1H1Pos_Z = 0     * cm;
  G4double pS1H2Pos_X   = 4.1  * cm;  G4double pS1H2Pos_Y   = 9.5   * cm;  G4double pS1H2Pos_Z = 0     * cm;
  G4double pS2LDPos_X   = 0.0  * cm;  G4double pS2LDPos_Y   = -7.5  * cm;  G4double pS2LDPos_Z = 610.0 * cm;
  G4double pS21BPos_X   = 0.0  * cm;  G4double pS21BPos_Y   = -13.3 * cm;  G4double pS21BPos_Z = 620.0 * cm;
  G4double pS21HPos_X   = 0.0  * cm;  G4double pS21HPos_Y   = -0.5  * cm;  G4double pS21HPos_Z = 0.0   * cm;
  G4double pS3LDPos_X   = 0.0  * cm;  G4double pS3LDPos_Y   = -7.5  * cm;  G4double pS3LDPos_Z = 723.0 * cm;

  G4VSolid* S1LDSolid  = new G4Box ( "S1LDBox"  , pS1LDHLX, pS1LDHLY,  pS1LDHLZ );
  G4VSolid* S1H1Solid  = new G4Para( "S1H1Para" , pS1H1HLX, pS1H1HLY,  pS1H1HLZ, pS1H1alpha, pS1H1theta, pS1H1phi);
  G4VSolid* S1H2Solid  = new G4Para( "S1H2Para" , pS1H2HLX, pS1H2HLY,  pS1H2HLZ, pS1H2alpha, pS1H2theta, pS1H2phi);
  G4VSolid* S2LDSolid  = new G4Box ( "S2LDBox"  , pS2LDHLX, pS2LDHLY,  pS2LDHLZ );
  G4VSolid* S21BSolid  = new G4Box ( "S21BBox"  , pS21BHLX, pS21BHLY,  pS21BHLZ );
  G4VSolid* S21HSolid  = new G4Box ( "S21HBox"  , pS21HHLX, pS21HHLY,  pS21HHLZ );
  G4VSolid* S3LDSolid  = new G4Box ( "S3LDBox"  , pS3LDHLX, pS3LDHLY,  pS3LDHLZ );

  G4SubtractionSolid* subSHL1 = new G4SubtractionSolid("subSHL1", S1LDSolid, S1H1Solid, pRot9, G4ThreeVector(pS1H1Pos_X, pS1H1Pos_Y, pS1H1Pos_Z) );
  G4SubtractionSolid* subSHL2 = new G4SubtractionSolid("subSHL2", subSHL1, S1H2Solid, pRot9, G4ThreeVector(pS1H2Pos_X, pS1H2Pos_Y, pS1H2Pos_Z) );
  G4SubtractionSolid* subSHL3 = new G4SubtractionSolid("subSHL3", S21BSolid, S21HSolid, 0, G4ThreeVector(pS21HPos_X, pS21HPos_Y, pS21HPos_Z) );

  // FIXME!!!!!  Something here may require attention. Are materials correct? Does not match VisAttributes
  G4LogicalVolume* SHL1Logical = new G4LogicalVolume(subSHL2, lead, "Shield1", 0,0,0);
  G4LogicalVolume* SHL2Logical = new G4LogicalVolume(S2LDSolid, lead, "Shield2", 0,0,0);
  G4LogicalVolume* SHL3Logical = new G4LogicalVolume(subSHL3, aluminum, "Shield3", 0,0,0);
  G4LogicalVolume* SHL4Logical = new G4LogicalVolume(S3LDSolid, lead, "Shield4", 0,0,0);
  // FIXME!!!!!  Something here may require attention. Are VisAttributes correct? Does not match materials.
  // SHL1Logical->SetVisAttributes(LeadVisAtt);
  // SHL2Logical->SetVisAttributes(AlumVisAtt);
  // SHL3Logical->SetVisAttributes(LeadVisAtt);
  // SHL4Logical->SetVisAttributes(LeadVisAtt);

  // new G4PVPlacement(0, G4ThreeVector(pS1LDPos_X, pS1LDPos_Y, pS1LDPos_Z), SHL1Logical, "Shield1", world_log, 0,0,fCheckOverlaps);
  // new G4PVPlacement(0, G4ThreeVector(pS2LDPos_X, pS2LDPos_Y, pS2LDPos_Z), SHL2Logical, "Shield2", world_log, 0,0,fCheckOverlaps);
  // new G4PVPlacement(0, G4ThreeVector(pS21BPos_X, pS21BPos_Y, pS21BPos_Z), SHL3Logical, "Shield3", world_log, 0,0,fCheckOverlaps);
  // new G4PVPlacement(0, G4ThreeVector(pS3LDPos_X, pS3LDPos_Y, pS3LDPos_Z), SHL4Logical, "Shield4", world_log, 0,0,fCheckOverlaps);



  return world_phys;
}

  //////////////////////////////////////////////////////////////  (╯°□°）╯︵ ┻━┻
  //
