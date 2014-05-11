//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 18 00:04:13 2011 by ROOT version 5.28/00c
// from TTree vdMScanData/vdMScanData
// found on file: 2011MarScan1Raw-v8.root
//////////////////////////////////////////////////////////

#ifndef vdMScanData_h
#define vdMScanData_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class vdMScanData {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    ULong64_t       SCANDATA_StartTime;
    ULong64_t       SCANDATA_EndTime;
    UInt_t          SCANDATA_ScanRun;
    UInt_t          SCANDATA_ScanLB;
    Float_t         SCANDATA_StepProgress;
    Float_t         SCANDATA_ScanningIP;
    Float_t         SCANDATA_AcquisitionFlag;
    Float_t         SCANDATA_MovingBeam;
    Float_t         SCANDATA_NominalSeparation;
    Float_t         SCANDATA_ScanInPlane;
    Float_t         BEAMPOSITION_B1PositionH;
    Float_t         BEAMPOSITION_B1PositionV;
    Float_t         BEAMPOSITION_B2PositionH;
    Float_t         BEAMPOSITION_B2PositionV;
    Float_t         BEAMPOSITION_B1AngleH;
    Float_t         BEAMPOSITION_B1AngleV;
    Float_t         BEAMPOSITION_B2AngleH;
    Float_t         BEAMPOSITION_B2AngleV;
    Int_t           FILLPARAMSB1_B1Bunches;
    Int_t           FILLPARAMSB1_B1BCIDs[76];   //[B1Bunches]
    Int_t           FILLPARAMSB2_B2Bunches;
    Int_t           FILLPARAMSB2_B2BCIDs[76];   //[B2Bunches]
    Int_t           FILLPARAMS_LuminousBunches;
    Int_t           FILLPARAMS_LuminousBCIDs[68];   //[LuminousBunches]
    Float_t         BPTX_LBDATA_BPTX_B1Intensity;
    Float_t         BPTX_LBDATA_BPTX_B2Intensity;
    Float_t         BPTX_LBDATA_BPTX_B1IntensityAll;
    Float_t         BPTX_LBDATA_BPTX_B2IntensityAll;
    Float_t         BCT_LBDATA_BCT_B1Intensity;
    Float_t         BCT_LBDATA_BCT_B2Intensity;
    Float_t         BCT_LBDATA_BCT_B1IntensityAll;
    Float_t         BCT_LBDATA_BCT_B2IntensityAll;
    UInt_t          BCT_BUNCHDATA_BCT_Valid;
    Float_t         BCT_BUNCHDATA_BCT_B1BunchAverage;
    Float_t         BCT_BUNCHDATA_BCT_B2BunchAverage;
    Int_t           BCT_BUNCHDATA_BCT_B1Bunches;
    Int_t           BCT_BUNCHDATA_BCT_B2Bunches;
    Int_t           BCT_B1BCID[148];   //[BCT_B1Bunches]
    Int_t           BCT_B2BCID[76];   //[BCT_B2Bunches]
    Float_t         BCT_B1BunchIntensity[148];   //[BCT_B1Bunches]
    Float_t         BCT_B2BunchIntensity[76];   //[BCT_B2Bunches]
    UInt_t          BPTX_BUNCHDATA_BPTX_Valid;
    Float_t         BPTX_BUNCHDATA_BPTX_B1BunchAverage;
    Float_t         BPTX_BUNCHDATA_BPTX_B2BunchAverage;
    Int_t           BPTX_BUNCHDATA_BPTX_B1Bunches;
    Int_t           BPTX_BUNCHDATA_BPTX_B2Bunches;
    Int_t           BPTX_B1BCID[76];   //[BPTX_B1Bunches]
    Int_t           BPTX_B2BCID[76];   //[BPTX_B2Bunches]
    Float_t         BPTX_B1BunchIntensity[76];   //[BPTX_B1Bunches]
    Float_t         BPTX_B2BunchIntensity[76];   //[BPTX_B2Bunches]
    UInt_t          lucidEvtAND_BUNCHDATA_lucidEvtAND_Channel;
    UInt_t          lucidEvtAND_BUNCHDATA_lucidEvtAND_Valid;
    Float_t         lucidEvtAND_BUNCHDATA_lucidEvtAND_AverageRawInstLum;
    Int_t           lucidEvtAND_BUNCHDATA_lucidEvtAND_LuminousBunches;
    Int_t           lucidEvtAND_BCID[94];   //[lucidEvtAND_LuminousBunches]
    Float_t         lucidEvtAND_BunchRawInstLum[94];   //[lucidEvtAND_LuminousBunches]
    UInt_t          lucidEvtOR_BUNCHDATA_lucidEvtOR_Channel;
    UInt_t          lucidEvtOR_BUNCHDATA_lucidEvtOR_Valid;
    Float_t         lucidEvtOR_BUNCHDATA_lucidEvtOR_AverageRawInstLum;
    Int_t           lucidEvtOR_BUNCHDATA_lucidEvtOR_LuminousBunches;
    Int_t           lucidEvtOR_BCID[3564];   //[lucidEvtOR_LuminousBunches]
    Float_t         lucidEvtOR_BunchRawInstLum[3564];   //[lucidEvtOR_LuminousBunches]
    UInt_t          lucidHitOR_BUNCHDATA_lucidHitOR_Channel;
    UInt_t          lucidHitOR_BUNCHDATA_lucidHitOR_Valid;
    Float_t         lucidHitOR_BUNCHDATA_lucidHitOR_AverageRawInstLum;
    Int_t           lucidHitOR_BUNCHDATA_lucidHitOR_LuminousBunches;
    Int_t           lucidHitOR_BCID[3564];   //[lucidHitOR_LuminousBunches]
    Float_t         lucidHitOR_BunchRawInstLum[3564];   //[lucidHitOR_LuminousBunches]
    UInt_t          lucidHitAND_BUNCHDATA_lucidHitAND_Channel;
    UInt_t          lucidHitAND_BUNCHDATA_lucidHitAND_Valid;
    Float_t         lucidHitAND_BUNCHDATA_lucidHitAND_AverageRawInstLum;
    Int_t           lucidHitAND_BUNCHDATA_lucidHitAND_LuminousBunches;
    Int_t           lucidHitAND_BCID[94];   //[lucidHitAND_LuminousBunches]
    Float_t         lucidHitAND_BunchRawInstLum[94];   //[lucidHitAND_LuminousBunches]
    UInt_t          lucidEvtA_BUNCHDATA_lucidEvtA_Channel;
    UInt_t          lucidEvtA_BUNCHDATA_lucidEvtA_Valid;
    Float_t         lucidEvtA_BUNCHDATA_lucidEvtA_AverageRawInstLum;
    Int_t           lucidEvtA_BUNCHDATA_lucidEvtA_LuminousBunches;
    Int_t           lucidEvtA_BCID[3564];   //[lucidEvtA_LuminousBunches]
    Float_t         lucidEvtA_BunchRawInstLum[3564];   //[lucidEvtA_LuminousBunches]
    UInt_t          lucidEvtC_BUNCHDATA_lucidEvtC_Channel;
    UInt_t          lucidEvtC_BUNCHDATA_lucidEvtC_Valid;
    Float_t         lucidEvtC_BUNCHDATA_lucidEvtC_AverageRawInstLum;
    Int_t           lucidEvtC_BUNCHDATA_lucidEvtC_LuminousBunches;
    Int_t           lucidEvtC_BCID[3564];   //[lucidEvtC_LuminousBunches]
    Float_t         lucidEvtC_BunchRawInstLum[3564];   //[lucidEvtC_LuminousBunches]
    UInt_t          bcmHXORC_BUNCHDATA_bcmHXORC_Channel;
    UInt_t          bcmHXORC_BUNCHDATA_bcmHXORC_Valid;
    Float_t         bcmHXORC_BUNCHDATA_bcmHXORC_AverageRawInstLum;
    Int_t           bcmHXORC_BUNCHDATA_bcmHXORC_LuminousBunches;
    Int_t           bcmHXORC_BCID[954];   //[bcmHXORC_LuminousBunches]
    Float_t         bcmHXORC_BunchRawInstLum[954];   //[bcmHXORC_LuminousBunches]
    UInt_t          bcmVEvtOR_BUNCHDATA_bcmVEvtOR_Channel;
    UInt_t          bcmVEvtOR_BUNCHDATA_bcmVEvtOR_Valid;
    Float_t         bcmVEvtOR_BUNCHDATA_bcmVEvtOR_AverageRawInstLum;
    Int_t           bcmVEvtOR_BUNCHDATA_bcmVEvtOR_LuminousBunches;
    Int_t           bcmVEvtOR_BCID[1522];   //[bcmVEvtOR_LuminousBunches]
    Float_t         bcmVEvtOR_BunchRawInstLum[1522];   //[bcmVEvtOR_LuminousBunches]
    UInt_t          bcmVEvtAND_BUNCHDATA_bcmVEvtAND_Channel;
    UInt_t          bcmVEvtAND_BUNCHDATA_bcmVEvtAND_Valid;
    Float_t         bcmVEvtAND_BUNCHDATA_bcmVEvtAND_AverageRawInstLum;
    Int_t           bcmVEvtAND_BUNCHDATA_bcmVEvtAND_LuminousBunches;
    Int_t           bcmVEvtAND_BCID[81];   //[bcmVEvtAND_LuminousBunches]
    Float_t         bcmVEvtAND_BunchRawInstLum[81];   //[bcmVEvtAND_LuminousBunches]
    UInt_t          bcmVXORC_BUNCHDATA_bcmVXORC_Channel;
    UInt_t          bcmVXORC_BUNCHDATA_bcmVXORC_Valid;
    Float_t         bcmVXORC_BUNCHDATA_bcmVXORC_AverageRawInstLum;
    Int_t           bcmVXORC_BUNCHDATA_bcmVXORC_LuminousBunches;
    Int_t           bcmVXORC_BCID[926];   //[bcmVXORC_LuminousBunches]
    Float_t         bcmVXORC_BunchRawInstLum[926];   //[bcmVXORC_LuminousBunches]
    UInt_t          bcmHEvtOR_BUNCHDATA_bcmHEvtOR_Channel;
    UInt_t          bcmHEvtOR_BUNCHDATA_bcmHEvtOR_Valid;
    Float_t         bcmHEvtOR_BUNCHDATA_bcmHEvtOR_AverageRawInstLum;
    Int_t           bcmHEvtOR_BUNCHDATA_bcmHEvtOR_LuminousBunches;
    Int_t           bcmHEvtOR_BCID[1394];   //[bcmHEvtOR_LuminousBunches]
    Float_t         bcmHEvtOR_BunchRawInstLum[1394];   //[bcmHEvtOR_LuminousBunches]
    UInt_t          bcmHEvtAND_BUNCHDATA_bcmHEvtAND_Channel;
    UInt_t          bcmHEvtAND_BUNCHDATA_bcmHEvtAND_Valid;
    Float_t         bcmHEvtAND_BUNCHDATA_bcmHEvtAND_AverageRawInstLum;
    Int_t           bcmHEvtAND_BUNCHDATA_bcmHEvtAND_LuminousBunches;
    Int_t           bcmHEvtAND_BCID[81];   //[bcmHEvtAND_LuminousBunches]
    Float_t         bcmHEvtAND_BunchRawInstLum[81];   //[bcmHEvtAND_LuminousBunches]
    UInt_t          pref_LUMI_pref_LumiChannel;
    Float_t         pref_LUMI_pref_LBAvInstLumPhys;
    Float_t         pref_LUMI_pref_LBAvEvtsPerBXPhys;
    Float_t         pref_LUMI_pref_LBAvRawInstLumPhys;
    Float_t         pref_LUMI_pref_LBAvInstLumAll;
    Float_t         pref_LUMI_pref_LBAvRawInstLumAll;
    Float_t         pref_LUMI_pref_LBAvEvtsPerBXAll;
    Float_t         pref_LUMI_pref_LBAvOLCInstLum;
    Float_t         pref_LUMI_pref_LBAvOLCRawInstLum;
    Float_t         pref_LUMI_pref_LBAvOLCEvtsPerBX;
    UInt_t          pref_LUMI_pref_NOrbPhys;
    UInt_t          pref_LUMI_pref_NOrbAll;
    UInt_t          pref_LUMI_pref_NOrbOLC;
    UInt_t          pref_LUMI_pref_DetectorState;
    UInt_t          pref_LUMI_pref_LumiValid;
    UInt_t          lucidEvtAND_LUMI_lucidEvtAND_LumiChannel;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvInstLumPhys;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvEvtsPerBXPhys;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvRawInstLumPhys;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvInstLumAll;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvRawInstLumAll;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvEvtsPerBXAll;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvOLCInstLum;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvOLCRawInstLum;
    Float_t         lucidEvtAND_LUMI_lucidEvtAND_LBAvOLCEvtsPerBX;
    UInt_t          lucidEvtAND_LUMI_lucidEvtAND_NOrbPhys;
    UInt_t          lucidEvtAND_LUMI_lucidEvtAND_NOrbAll;
    UInt_t          lucidEvtAND_LUMI_lucidEvtAND_NOrbOLC;
    UInt_t          lucidEvtAND_LUMI_lucidEvtAND_DetectorState;
    UInt_t          lucidEvtAND_LUMI_lucidEvtAND_LumiValid;
    UInt_t          lucidEvtOR_LUMI_lucidEvtOR_LumiChannel;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvInstLumPhys;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvEvtsPerBXPhys;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvRawInstLumPhys;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvInstLumAll;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvRawInstLumAll;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvEvtsPerBXAll;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvOLCInstLum;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvOLCRawInstLum;
    Float_t         lucidEvtOR_LUMI_lucidEvtOR_LBAvOLCEvtsPerBX;
    UInt_t          lucidEvtOR_LUMI_lucidEvtOR_NOrbPhys;
    UInt_t          lucidEvtOR_LUMI_lucidEvtOR_NOrbAll;
    UInt_t          lucidEvtOR_LUMI_lucidEvtOR_NOrbOLC;
    UInt_t          lucidEvtOR_LUMI_lucidEvtOR_DetectorState;
    UInt_t          lucidEvtOR_LUMI_lucidEvtOR_LumiValid;
    UInt_t          lucidHitOR_LUMI_lucidHitOR_LumiChannel;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvInstLumPhys;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvEvtsPerBXPhys;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvRawInstLumPhys;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvInstLumAll;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvRawInstLumAll;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvEvtsPerBXAll;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvOLCInstLum;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvOLCRawInstLum;
    Float_t         lucidHitOR_LUMI_lucidHitOR_LBAvOLCEvtsPerBX;
    UInt_t          lucidHitOR_LUMI_lucidHitOR_NOrbPhys;
    UInt_t          lucidHitOR_LUMI_lucidHitOR_NOrbAll;
    UInt_t          lucidHitOR_LUMI_lucidHitOR_NOrbOLC;
    UInt_t          lucidHitOR_LUMI_lucidHitOR_DetectorState;
    UInt_t          lucidHitOR_LUMI_lucidHitOR_LumiValid;
    UInt_t          lucidHitAND_LUMI_lucidHitAND_LumiChannel;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvInstLumPhys;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvEvtsPerBXPhys;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvRawInstLumPhys;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvInstLumAll;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvRawInstLumAll;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvEvtsPerBXAll;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvOLCInstLum;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvOLCRawInstLum;
    Float_t         lucidHitAND_LUMI_lucidHitAND_LBAvOLCEvtsPerBX;
    UInt_t          lucidHitAND_LUMI_lucidHitAND_NOrbPhys;
    UInt_t          lucidHitAND_LUMI_lucidHitAND_NOrbAll;
    UInt_t          lucidHitAND_LUMI_lucidHitAND_NOrbOLC;
    UInt_t          lucidHitAND_LUMI_lucidHitAND_DetectorState;
    UInt_t          lucidHitAND_LUMI_lucidHitAND_LumiValid;
    UInt_t          lucidEvtA_LUMI_lucidEvtA_LumiChannel;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvInstLumPhys;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvEvtsPerBXPhys;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvRawInstLumPhys;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvInstLumAll;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvRawInstLumAll;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvEvtsPerBXAll;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvOLCInstLum;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvOLCRawInstLum;
    Float_t         lucidEvtA_LUMI_lucidEvtA_LBAvOLCEvtsPerBX;
    UInt_t          lucidEvtA_LUMI_lucidEvtA_NOrbPhys;
    UInt_t          lucidEvtA_LUMI_lucidEvtA_NOrbAll;
    UInt_t          lucidEvtA_LUMI_lucidEvtA_NOrbOLC;
    UInt_t          lucidEvtA_LUMI_lucidEvtA_DetectorState;
    UInt_t          lucidEvtA_LUMI_lucidEvtA_LumiValid;
    UInt_t          lucidEvtC_LUMI_lucidEvtC_LumiChannel;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvInstLumPhys;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvEvtsPerBXPhys;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvRawInstLumPhys;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvInstLumAll;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvRawInstLumAll;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvEvtsPerBXAll;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvOLCInstLum;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvOLCRawInstLum;
    Float_t         lucidEvtC_LUMI_lucidEvtC_LBAvOLCEvtsPerBX;
    UInt_t          lucidEvtC_LUMI_lucidEvtC_NOrbPhys;
    UInt_t          lucidEvtC_LUMI_lucidEvtC_NOrbAll;
    UInt_t          lucidEvtC_LUMI_lucidEvtC_NOrbOLC;
    UInt_t          lucidEvtC_LUMI_lucidEvtC_DetectorState;
    UInt_t          lucidEvtC_LUMI_lucidEvtC_LumiValid;
    UInt_t          bcmHXORC_LUMI_bcmHXORC_LumiChannel;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvInstLumPhys;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvEvtsPerBXPhys;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvRawInstLumPhys;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvInstLumAll;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvRawInstLumAll;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvEvtsPerBXAll;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvOLCInstLum;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvOLCRawInstLum;
    Float_t         bcmHXORC_LUMI_bcmHXORC_LBAvOLCEvtsPerBX;
    UInt_t          bcmHXORC_LUMI_bcmHXORC_NOrbPhys;
    UInt_t          bcmHXORC_LUMI_bcmHXORC_NOrbAll;
    UInt_t          bcmHXORC_LUMI_bcmHXORC_NOrbOLC;
    UInt_t          bcmHXORC_LUMI_bcmHXORC_DetectorState;
    UInt_t          bcmHXORC_LUMI_bcmHXORC_LumiValid;
    UInt_t          mbtsEvtOR_LUMI_mbtsEvtOR_LumiChannel;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvInstLumPhys;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvEvtsPerBXPhys;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvRawInstLumPhys;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvInstLumAll;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvRawInstLumAll;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvEvtsPerBXAll;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvOLCInstLum;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvOLCRawInstLum;
    Float_t         mbtsEvtOR_LUMI_mbtsEvtOR_LBAvOLCEvtsPerBX;
    UInt_t          mbtsEvtOR_LUMI_mbtsEvtOR_NOrbPhys;
    UInt_t          mbtsEvtOR_LUMI_mbtsEvtOR_NOrbAll;
    UInt_t          mbtsEvtOR_LUMI_mbtsEvtOR_NOrbOLC;
    UInt_t          mbtsEvtOR_LUMI_mbtsEvtOR_DetectorState;
    UInt_t          mbtsEvtOR_LUMI_mbtsEvtOR_LumiValid;
    UInt_t          mbtsEvtAND_LUMI_mbtsEvtAND_LumiChannel;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvInstLumPhys;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvEvtsPerBXPhys;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvRawInstLumPhys;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvInstLumAll;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvRawInstLumAll;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvEvtsPerBXAll;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvOLCInstLum;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvOLCRawInstLum;
    Float_t         mbtsEvtAND_LUMI_mbtsEvtAND_LBAvOLCEvtsPerBX;
    UInt_t          mbtsEvtAND_LUMI_mbtsEvtAND_NOrbPhys;
    UInt_t          mbtsEvtAND_LUMI_mbtsEvtAND_NOrbAll;
    UInt_t          mbtsEvtAND_LUMI_mbtsEvtAND_NOrbOLC;
    UInt_t          mbtsEvtAND_LUMI_mbtsEvtAND_DetectorState;
    UInt_t          mbtsEvtAND_LUMI_mbtsEvtAND_LumiValid;
    UInt_t          mbtsHitOR_LUMI_mbtsHitOR_LumiChannel;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvInstLumPhys;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvEvtsPerBXPhys;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvRawInstLumPhys;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvInstLumAll;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvRawInstLumAll;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvEvtsPerBXAll;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvOLCInstLum;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvOLCRawInstLum;
    Float_t         mbtsHitOR_LUMI_mbtsHitOR_LBAvOLCEvtsPerBX;
    UInt_t          mbtsHitOR_LUMI_mbtsHitOR_NOrbPhys;
    UInt_t          mbtsHitOR_LUMI_mbtsHitOR_NOrbAll;
    UInt_t          mbtsHitOR_LUMI_mbtsHitOR_NOrbOLC;
    UInt_t          mbtsHitOR_LUMI_mbtsHitOR_DetectorState;
    UInt_t          mbtsHitOR_LUMI_mbtsHitOR_LumiValid;
    UInt_t          zdcEvtAND_LUMI_zdcEvtAND_LumiChannel;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvInstLumPhys;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvEvtsPerBXPhys;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvRawInstLumPhys;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvInstLumAll;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvRawInstLumAll;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvEvtsPerBXAll;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvOLCInstLum;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvOLCRawInstLum;
    Float_t         zdcEvtAND_LUMI_zdcEvtAND_LBAvOLCEvtsPerBX;
    UInt_t          zdcEvtAND_LUMI_zdcEvtAND_NOrbPhys;
    UInt_t          zdcEvtAND_LUMI_zdcEvtAND_NOrbAll;
    UInt_t          zdcEvtAND_LUMI_zdcEvtAND_NOrbOLC;
    UInt_t          zdcEvtAND_LUMI_zdcEvtAND_DetectorState;
    UInt_t          zdcEvtAND_LUMI_zdcEvtAND_LumiValid;
    UInt_t          zdcEvtORA_LUMI_zdcEvtORA_LumiChannel;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvInstLumPhys;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvEvtsPerBXPhys;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvRawInstLumPhys;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvInstLumAll;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvRawInstLumAll;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvEvtsPerBXAll;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvOLCInstLum;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvOLCRawInstLum;
    Float_t         zdcEvtORA_LUMI_zdcEvtORA_LBAvOLCEvtsPerBX;
    UInt_t          zdcEvtORA_LUMI_zdcEvtORA_NOrbPhys;
    UInt_t          zdcEvtORA_LUMI_zdcEvtORA_NOrbAll;
    UInt_t          zdcEvtORA_LUMI_zdcEvtORA_NOrbOLC;
    UInt_t          zdcEvtORA_LUMI_zdcEvtORA_DetectorState;
    UInt_t          zdcEvtORA_LUMI_zdcEvtORA_LumiValid;
    UInt_t          bcmVEvtOR_LUMI_bcmVEvtOR_LumiChannel;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvInstLumPhys;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvEvtsPerBXPhys;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvRawInstLumPhys;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvInstLumAll;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvRawInstLumAll;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvEvtsPerBXAll;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvOLCInstLum;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvOLCRawInstLum;
    Float_t         bcmVEvtOR_LUMI_bcmVEvtOR_LBAvOLCEvtsPerBX;
    UInt_t          bcmVEvtOR_LUMI_bcmVEvtOR_NOrbPhys;
    UInt_t          bcmVEvtOR_LUMI_bcmVEvtOR_NOrbAll;
    UInt_t          bcmVEvtOR_LUMI_bcmVEvtOR_NOrbOLC;
    UInt_t          bcmVEvtOR_LUMI_bcmVEvtOR_DetectorState;
    UInt_t          bcmVEvtOR_LUMI_bcmVEvtOR_LumiValid;
    UInt_t          bcmVEvtAND_LUMI_bcmVEvtAND_LumiChannel;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvInstLumPhys;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvEvtsPerBXPhys;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvRawInstLumPhys;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvInstLumAll;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvRawInstLumAll;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvEvtsPerBXAll;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvOLCInstLum;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvOLCRawInstLum;
    Float_t         bcmVEvtAND_LUMI_bcmVEvtAND_LBAvOLCEvtsPerBX;
    UInt_t          bcmVEvtAND_LUMI_bcmVEvtAND_NOrbPhys;
    UInt_t          bcmVEvtAND_LUMI_bcmVEvtAND_NOrbAll;
    UInt_t          bcmVEvtAND_LUMI_bcmVEvtAND_NOrbOLC;
    UInt_t          bcmVEvtAND_LUMI_bcmVEvtAND_DetectorState;
    UInt_t          bcmVEvtAND_LUMI_bcmVEvtAND_LumiValid;
    UInt_t          bcmVXORC_LUMI_bcmVXORC_LumiChannel;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvInstLumPhys;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvEvtsPerBXPhys;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvRawInstLumPhys;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvInstLumAll;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvRawInstLumAll;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvEvtsPerBXAll;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvOLCInstLum;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvOLCRawInstLum;
    Float_t         bcmVXORC_LUMI_bcmVXORC_LBAvOLCEvtsPerBX;
    UInt_t          bcmVXORC_LUMI_bcmVXORC_NOrbPhys;
    UInt_t          bcmVXORC_LUMI_bcmVXORC_NOrbAll;
    UInt_t          bcmVXORC_LUMI_bcmVXORC_NOrbOLC;
    UInt_t          bcmVXORC_LUMI_bcmVXORC_DetectorState;
    UInt_t          bcmVXORC_LUMI_bcmVXORC_LumiValid;
    UInt_t          bcmHEvtOR_LUMI_bcmHEvtOR_LumiChannel;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvInstLumPhys;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvEvtsPerBXPhys;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvRawInstLumPhys;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvInstLumAll;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvRawInstLumAll;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvEvtsPerBXAll;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvOLCInstLum;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvOLCRawInstLum;
    Float_t         bcmHEvtOR_LUMI_bcmHEvtOR_LBAvOLCEvtsPerBX;
    UInt_t          bcmHEvtOR_LUMI_bcmHEvtOR_NOrbPhys;
    UInt_t          bcmHEvtOR_LUMI_bcmHEvtOR_NOrbAll;
    UInt_t          bcmHEvtOR_LUMI_bcmHEvtOR_NOrbOLC;
    UInt_t          bcmHEvtOR_LUMI_bcmHEvtOR_DetectorState;
    UInt_t          bcmHEvtOR_LUMI_bcmHEvtOR_LumiValid;
    UInt_t          zdcEvtORC_LUMI_zdcEvtORC_LumiChannel;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvInstLumPhys;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvEvtsPerBXPhys;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvRawInstLumPhys;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvInstLumAll;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvRawInstLumAll;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvEvtsPerBXAll;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvOLCInstLum;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvOLCRawInstLum;
    Float_t         zdcEvtORC_LUMI_zdcEvtORC_LBAvOLCEvtsPerBX;
    UInt_t          zdcEvtORC_LUMI_zdcEvtORC_NOrbPhys;
    UInt_t          zdcEvtORC_LUMI_zdcEvtORC_NOrbAll;
    UInt_t          zdcEvtORC_LUMI_zdcEvtORC_NOrbOLC;
    UInt_t          zdcEvtORC_LUMI_zdcEvtORC_DetectorState;
    UInt_t          zdcEvtORC_LUMI_zdcEvtORC_LumiValid;
    UInt_t          bcmHEvtAND_LUMI_bcmHEvtAND_LumiChannel;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvInstLumPhys;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvEvtsPerBXPhys;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvRawInstLumPhys;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvInstLumAll;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvRawInstLumAll;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvEvtsPerBXAll;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvOLCInstLum;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvOLCRawInstLum;
    Float_t         bcmHEvtAND_LUMI_bcmHEvtAND_LBAvOLCEvtsPerBX;
    UInt_t          bcmHEvtAND_LUMI_bcmHEvtAND_NOrbPhys;
    UInt_t          bcmHEvtAND_LUMI_bcmHEvtAND_NOrbAll;
    UInt_t          bcmHEvtAND_LUMI_bcmHEvtAND_NOrbOLC;
    UInt_t          bcmHEvtAND_LUMI_bcmHEvtAND_DetectorState;
    UInt_t          bcmHEvtAND_LUMI_bcmHEvtAND_LumiValid;
    Float_t         lucidEvtAND_TURN_lucidEvtAND_BXPhys;
    Float_t         lucidEvtAND_TURN_lucidEvtAND_CountsPhys;
    Float_t         lucidEvtAND_TURN_lucidEvtAND_BXAll;
    Float_t         lucidEvtAND_TURN_lucidEvtAND_CountsAll;
    Float_t         lucidEvtAND_TURN_lucidEvtAND_BXOLC;
    Float_t         lucidEvtAND_TURN_lucidEvtAND_CountsOLC;
    Float_t         lucidEvtOR_TURN_lucidEvtOR_BXPhys;
    Float_t         lucidEvtOR_TURN_lucidEvtOR_CountsPhys;
    Float_t         lucidEvtOR_TURN_lucidEvtOR_BXAll;
    Float_t         lucidEvtOR_TURN_lucidEvtOR_CountsAll;
    Float_t         lucidEvtOR_TURN_lucidEvtOR_BXOLC;
    Float_t         lucidEvtOR_TURN_lucidEvtOR_CountsOLC;
    Float_t         lucidHitOR_TURN_lucidHitOR_BXPhys;
    Float_t         lucidHitOR_TURN_lucidHitOR_CountsPhys;
    Float_t         lucidHitOR_TURN_lucidHitOR_BXAll;
    Float_t         lucidHitOR_TURN_lucidHitOR_CountsAll;
    Float_t         lucidHitOR_TURN_lucidHitOR_BXOLC;
    Float_t         lucidHitOR_TURN_lucidHitOR_CountsOLC;
    Float_t         lucidHitAND_TURN_lucidHitAND_BXPhys;
    Float_t         lucidHitAND_TURN_lucidHitAND_CountsPhys;
    Float_t         lucidHitAND_TURN_lucidHitAND_BXAll;
    Float_t         lucidHitAND_TURN_lucidHitAND_CountsAll;
    Float_t         lucidHitAND_TURN_lucidHitAND_BXOLC;
    Float_t         lucidHitAND_TURN_lucidHitAND_CountsOLC;
    Float_t         lucidEvtA_TURN_lucidEvtA_BXPhys;
    Float_t         lucidEvtA_TURN_lucidEvtA_CountsPhys;
    Float_t         lucidEvtA_TURN_lucidEvtA_BXAll;
    Float_t         lucidEvtA_TURN_lucidEvtA_CountsAll;
    Float_t         lucidEvtA_TURN_lucidEvtA_BXOLC;
    Float_t         lucidEvtA_TURN_lucidEvtA_CountsOLC;
    Float_t         lucidEvtC_TURN_lucidEvtC_BXPhys;
    Float_t         lucidEvtC_TURN_lucidEvtC_CountsPhys;
    Float_t         lucidEvtC_TURN_lucidEvtC_BXAll;
    Float_t         lucidEvtC_TURN_lucidEvtC_CountsAll;
    Float_t         lucidEvtC_TURN_lucidEvtC_BXOLC;
    Float_t         lucidEvtC_TURN_lucidEvtC_CountsOLC;
    Float_t         bcmHXORC_TURN_bcmHXORC_BXPhys;
    Float_t         bcmHXORC_TURN_bcmHXORC_CountsPhys;
    Float_t         bcmHXORC_TURN_bcmHXORC_BXAll;
    Float_t         bcmHXORC_TURN_bcmHXORC_CountsAll;
    Float_t         bcmHXORC_TURN_bcmHXORC_BXOLC;
    Float_t         bcmHXORC_TURN_bcmHXORC_CountsOLC;
    Float_t         mbtsEvtOR_TURN_mbtsEvtOR_BXPhys;
    Float_t         mbtsEvtOR_TURN_mbtsEvtOR_CountsPhys;
    Float_t         mbtsEvtOR_TURN_mbtsEvtOR_BXAll;
    Float_t         mbtsEvtOR_TURN_mbtsEvtOR_CountsAll;
    Float_t         mbtsEvtOR_TURN_mbtsEvtOR_BXOLC;
    Float_t         mbtsEvtOR_TURN_mbtsEvtOR_CountsOLC;
    Float_t         mbtsEvtAND_TURN_mbtsEvtAND_BXPhys;
    Float_t         mbtsEvtAND_TURN_mbtsEvtAND_CountsPhys;
    Float_t         mbtsEvtAND_TURN_mbtsEvtAND_BXAll;
    Float_t         mbtsEvtAND_TURN_mbtsEvtAND_CountsAll;
    Float_t         mbtsEvtAND_TURN_mbtsEvtAND_BXOLC;
    Float_t         mbtsEvtAND_TURN_mbtsEvtAND_CountsOLC;
    Float_t         mbtsHitOR_TURN_mbtsHitOR_BXPhys;
    Float_t         mbtsHitOR_TURN_mbtsHitOR_CountsPhys;
    Float_t         mbtsHitOR_TURN_mbtsHitOR_BXAll;
    Float_t         mbtsHitOR_TURN_mbtsHitOR_CountsAll;
    Float_t         mbtsHitOR_TURN_mbtsHitOR_BXOLC;
    Float_t         mbtsHitOR_TURN_mbtsHitOR_CountsOLC;
    Float_t         zdcEvtAND_TURN_zdcEvtAND_BXPhys;
    Float_t         zdcEvtAND_TURN_zdcEvtAND_CountsPhys;
    Float_t         zdcEvtAND_TURN_zdcEvtAND_BXAll;
    Float_t         zdcEvtAND_TURN_zdcEvtAND_CountsAll;
    Float_t         zdcEvtAND_TURN_zdcEvtAND_BXOLC;
    Float_t         zdcEvtAND_TURN_zdcEvtAND_CountsOLC;
    Float_t         zdcEvtORA_TURN_zdcEvtORA_BXPhys;
    Float_t         zdcEvtORA_TURN_zdcEvtORA_CountsPhys;
    Float_t         zdcEvtORA_TURN_zdcEvtORA_BXAll;
    Float_t         zdcEvtORA_TURN_zdcEvtORA_CountsAll;
    Float_t         zdcEvtORA_TURN_zdcEvtORA_BXOLC;
    Float_t         zdcEvtORA_TURN_zdcEvtORA_CountsOLC;
    Float_t         bcmVEvtOR_TURN_bcmVEvtOR_BXPhys;
    Float_t         bcmVEvtOR_TURN_bcmVEvtOR_CountsPhys;
    Float_t         bcmVEvtOR_TURN_bcmVEvtOR_BXAll;
    Float_t         bcmVEvtOR_TURN_bcmVEvtOR_CountsAll;
    Float_t         bcmVEvtOR_TURN_bcmVEvtOR_BXOLC;
    Float_t         bcmVEvtOR_TURN_bcmVEvtOR_CountsOLC;
    Float_t         bcmVEvtAND_TURN_bcmVEvtAND_BXPhys;
    Float_t         bcmVEvtAND_TURN_bcmVEvtAND_CountsPhys;
    Float_t         bcmVEvtAND_TURN_bcmVEvtAND_BXAll;
    Float_t         bcmVEvtAND_TURN_bcmVEvtAND_CountsAll;
    Float_t         bcmVEvtAND_TURN_bcmVEvtAND_BXOLC;
    Float_t         bcmVEvtAND_TURN_bcmVEvtAND_CountsOLC;
    Float_t         bcmVXORC_TURN_bcmVXORC_BXPhys;
    Float_t         bcmVXORC_TURN_bcmVXORC_CountsPhys;
    Float_t         bcmVXORC_TURN_bcmVXORC_BXAll;
    Float_t         bcmVXORC_TURN_bcmVXORC_CountsAll;
    Float_t         bcmVXORC_TURN_bcmVXORC_BXOLC;
    Float_t         bcmVXORC_TURN_bcmVXORC_CountsOLC;
    Float_t         bcmHEvtOR_TURN_bcmHEvtOR_BXPhys;
    Float_t         bcmHEvtOR_TURN_bcmHEvtOR_CountsPhys;
    Float_t         bcmHEvtOR_TURN_bcmHEvtOR_BXAll;
    Float_t         bcmHEvtOR_TURN_bcmHEvtOR_CountsAll;
    Float_t         bcmHEvtOR_TURN_bcmHEvtOR_BXOLC;
    Float_t         bcmHEvtOR_TURN_bcmHEvtOR_CountsOLC;
    Float_t         zdcEvtORC_TURN_zdcEvtORC_BXPhys;
    Float_t         zdcEvtORC_TURN_zdcEvtORC_CountsPhys;
    Float_t         zdcEvtORC_TURN_zdcEvtORC_BXAll;
    Float_t         zdcEvtORC_TURN_zdcEvtORC_CountsAll;
    Float_t         zdcEvtORC_TURN_zdcEvtORC_BXOLC;
    Float_t         zdcEvtORC_TURN_zdcEvtORC_CountsOLC;
    Float_t         bcmHEvtAND_TURN_bcmHEvtAND_BXPhys;
    Float_t         bcmHEvtAND_TURN_bcmHEvtAND_CountsPhys;
    Float_t         bcmHEvtAND_TURN_bcmHEvtAND_BXAll;
    Float_t         bcmHEvtAND_TURN_bcmHEvtAND_CountsAll;
    Float_t         bcmHEvtAND_TURN_bcmHEvtAND_BXOLC;
    Float_t         bcmHEvtAND_TURN_bcmHEvtAND_CountsOLC;

    // List of branches
    TBranch        *b_SCANDATA;   //!
    TBranch        *b_BEAMPOSITION;   //!
    TBranch        *b_FILLPARAMSB1;   //!
    TBranch        *b_FILLPARAMSB2;   //!
    TBranch        *b_FILLPARAMS;   //!
    TBranch        *b_BPTX_LBDATA;   //!
    TBranch        *b_BCT_LBDATA;   //!
    TBranch        *b_BCT_BUNCHDATA;   //!
    TBranch        *b_BCT_B1BCID;   //!
    TBranch        *b_BCT_B2BCID;   //!
    TBranch        *b_BCT_B1BunchIntensity;   //!
    TBranch        *b_BCT_B2BunchIntensity;   //!
    TBranch        *b_BPTX_BUNCHDATA;   //!
    TBranch        *b_BPTX_B1BCID;   //!
    TBranch        *b_BPTX_B2BCID;   //!
    TBranch        *b_BPTX_B1BunchIntensity;   //!
    TBranch        *b_BPTX_B2BunchIntensity;   //!
    TBranch        *b_lucidEvtAND_BUNCHDATA;   //!
    TBranch        *b_lucidEvtAND_BCID;   //!
    TBranch        *b_lucidEvtAND_BunchRawInstLum;   //!
    TBranch        *b_lucidEvtOR_BUNCHDATA;   //!
    TBranch        *b_lucidEvtOR_BCID;   //!
    TBranch        *b_lucidEvtOR_BunchRawInstLum;   //!
    TBranch        *b_lucidHitOR_BUNCHDATA;   //!
    TBranch        *b_lucidHitOR_BCID;   //!
    TBranch        *b_lucidHitOR_BunchRawInstLum;   //!
    TBranch        *b_lucidHitAND_BUNCHDATA;   //!
    TBranch        *b_lucidHitAND_BCID;   //!
    TBranch        *b_lucidHitAND_BunchRawInstLum;   //!
    TBranch        *b_lucidEvtA_BUNCHDATA;   //!
    TBranch        *b_lucidEvtA_BCID;   //!
    TBranch        *b_lucidEvtA_BunchRawInstLum;   //!
    TBranch        *b_lucidEvtC_BUNCHDATA;   //!
    TBranch        *b_lucidEvtC_BCID;   //!
    TBranch        *b_lucidEvtC_BunchRawInstLum;   //!
    TBranch        *b_bcmHXORC_BUNCHDATA;   //!
    TBranch        *b_bcmHXORC_BCID;   //!
    TBranch        *b_bcmHXORC_BunchRawInstLum;   //!
    TBranch        *b_bcmVEvtOR_BUNCHDATA;   //!
    TBranch        *b_bcmVEvtOR_BCID;   //!
    TBranch        *b_bcmVEvtOR_BunchRawInstLum;   //!
    TBranch        *b_bcmVEvtAND_BUNCHDATA;   //!
    TBranch        *b_bcmVEvtAND_BCID;   //!
    TBranch        *b_bcmVEvtAND_BunchRawInstLum;   //!
    TBranch        *b_bcmVXORC_BUNCHDATA;   //!
    TBranch        *b_bcmVXORC_BCID;   //!
    TBranch        *b_bcmVXORC_BunchRawInstLum;   //!
    TBranch        *b_bcmHEvtOR_BUNCHDATA;   //!
    TBranch        *b_bcmHEvtOR_BCID;   //!
    TBranch        *b_bcmHEvtOR_BunchRawInstLum;   //!
    TBranch        *b_bcmHEvtAND_BUNCHDATA;   //!
    TBranch        *b_bcmHEvtAND_BCID;   //!
    TBranch        *b_bcmHEvtAND_BunchRawInstLum;   //!
    TBranch        *b_pref_LUMI;   //!
    TBranch        *b_lucidEvtAND_LUMI;   //!
    TBranch        *b_lucidEvtOR_LUMI;   //!
    TBranch        *b_lucidHitOR_LUMI;   //!
    TBranch        *b_lucidHitAND_LUMI;   //!
    TBranch        *b_lucidEvtA_LUMI;   //!
    TBranch        *b_lucidEvtC_LUMI;   //!
    TBranch        *b_bcmHXORC_LUMI;   //!
    TBranch        *b_mbtsEvtOR_LUMI;   //!
    TBranch        *b_mbtsEvtAND_LUMI;   //!
    TBranch        *b_mbtsHitOR_LUMI;   //!
    TBranch        *b_zdcEvtAND_LUMI;   //!
    TBranch        *b_zdcEvtORA_LUMI;   //!
    TBranch        *b_bcmVEvtOR_LUMI;   //!
    TBranch        *b_bcmVEvtAND_LUMI;   //!
    TBranch        *b_bcmVXORC_LUMI;   //!
    TBranch        *b_bcmHEvtOR_LUMI;   //!
    TBranch        *b_zdcEvtORC_LUMI;   //!
    TBranch        *b_bcmHEvtAND_LUMI;   //!
    TBranch        *b_lucidEvtAND_TURN;   //!
    TBranch        *b_lucidEvtOR_TURN;   //!
    TBranch        *b_lucidHitOR_TURN;   //!
    TBranch        *b_lucidHitAND_TURN;   //!
    TBranch        *b_lucidEvtA_TURN;   //!
    TBranch        *b_lucidEvtC_TURN;   //!
    TBranch        *b_bcmHXORC_TURN;   //!
    TBranch        *b_mbtsEvtOR_TURN;   //!
    TBranch        *b_mbtsEvtAND_TURN;   //!
    TBranch        *b_mbtsHitOR_TURN;   //!
    TBranch        *b_zdcEvtAND_TURN;   //!
    TBranch        *b_zdcEvtORA_TURN;   //!
    TBranch        *b_bcmVEvtOR_TURN;   //!
    TBranch        *b_bcmVEvtAND_TURN;   //!
    TBranch        *b_bcmVXORC_TURN;   //!
    TBranch        *b_bcmHEvtOR_TURN;   //!
    TBranch        *b_zdcEvtORC_TURN;   //!
    TBranch        *b_bcmHEvtAND_TURN;   //!

    vdMScanData(TTree *tree=0);
    virtual ~vdMScanData();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef vdMScanData_cxx
vdMScanData::vdMScanData(TTree *tree) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("2011MarScan1Raw-v8.root");
    if (!f) {
      f = new TFile("2011MarScan1Raw-v8.root");
    }
    tree = (TTree*)gDirectory->Get("vdMScanData");

  }
  Init(tree);
}

vdMScanData::~vdMScanData() {
  if (!fChain) {
    return;
  }
  delete fChain->GetCurrentFile();
}

Int_t vdMScanData::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) {
    return 0;
  }
  return fChain->GetEntry(entry);
}
Long64_t vdMScanData::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) {
    return -5;
  }
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) {
    return centry;
  }
  if (!fChain->InheritsFrom(TChain::Class())) {
    return centry;
  }
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void vdMScanData::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) {
    return;
  }
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("SCANDATA", &SCANDATA_StartTime, &b_SCANDATA);
  fChain->SetBranchAddress("BEAMPOSITION", &BEAMPOSITION_B1PositionH, &b_BEAMPOSITION);
  fChain->SetBranchAddress("FILLPARAMSB1", &FILLPARAMSB1_B1Bunches, &b_FILLPARAMSB1);
  fChain->SetBranchAddress("FILLPARAMSB2", &FILLPARAMSB2_B2Bunches, &b_FILLPARAMSB2);
  fChain->SetBranchAddress("FILLPARAMS", &FILLPARAMS_LuminousBunches, &b_FILLPARAMS);
  fChain->SetBranchAddress("BPTX_LBDATA", &BPTX_LBDATA_BPTX_B1Intensity, &b_BPTX_LBDATA);
  fChain->SetBranchAddress("BCT_LBDATA", &BCT_LBDATA_BCT_B1Intensity, &b_BCT_LBDATA);
  fChain->SetBranchAddress("BCT_BUNCHDATA", &BCT_BUNCHDATA_BCT_Valid, &b_BCT_BUNCHDATA);
  fChain->SetBranchAddress("BCT_B1BCID", BCT_B1BCID, &b_BCT_B1BCID);
  fChain->SetBranchAddress("BCT_B2BCID", BCT_B2BCID, &b_BCT_B2BCID);
  fChain->SetBranchAddress("BCT_B1BunchIntensity", BCT_B1BunchIntensity, &b_BCT_B1BunchIntensity);
  fChain->SetBranchAddress("BCT_B2BunchIntensity", BCT_B2BunchIntensity, &b_BCT_B2BunchIntensity);
  fChain->SetBranchAddress("BPTX_BUNCHDATA", &BPTX_BUNCHDATA_BPTX_Valid, &b_BPTX_BUNCHDATA);
  fChain->SetBranchAddress("BPTX_B1BCID", BPTX_B1BCID, &b_BPTX_B1BCID);
  fChain->SetBranchAddress("BPTX_B2BCID", BPTX_B2BCID, &b_BPTX_B2BCID);
  fChain->SetBranchAddress("BPTX_B1BunchIntensity", BPTX_B1BunchIntensity, &b_BPTX_B1BunchIntensity);
  fChain->SetBranchAddress("BPTX_B2BunchIntensity", BPTX_B2BunchIntensity, &b_BPTX_B2BunchIntensity);
  fChain->SetBranchAddress("lucidEvtAND_BUNCHDATA", &lucidEvtAND_BUNCHDATA_lucidEvtAND_Channel, &b_lucidEvtAND_BUNCHDATA);
  fChain->SetBranchAddress("lucidEvtAND_BCID", lucidEvtAND_BCID, &b_lucidEvtAND_BCID);
  fChain->SetBranchAddress("lucidEvtAND_BunchRawInstLum", lucidEvtAND_BunchRawInstLum, &b_lucidEvtAND_BunchRawInstLum);
  fChain->SetBranchAddress("lucidEvtOR_BUNCHDATA", &lucidEvtOR_BUNCHDATA_lucidEvtOR_Channel, &b_lucidEvtOR_BUNCHDATA);
  fChain->SetBranchAddress("lucidEvtOR_BCID", lucidEvtOR_BCID, &b_lucidEvtOR_BCID);
  fChain->SetBranchAddress("lucidEvtOR_BunchRawInstLum", lucidEvtOR_BunchRawInstLum, &b_lucidEvtOR_BunchRawInstLum);
  fChain->SetBranchAddress("lucidHitOR_BUNCHDATA", &lucidHitOR_BUNCHDATA_lucidHitOR_Channel, &b_lucidHitOR_BUNCHDATA);
  fChain->SetBranchAddress("lucidHitOR_BCID", lucidHitOR_BCID, &b_lucidHitOR_BCID);
  fChain->SetBranchAddress("lucidHitOR_BunchRawInstLum", lucidHitOR_BunchRawInstLum, &b_lucidHitOR_BunchRawInstLum);
  fChain->SetBranchAddress("lucidHitAND_BUNCHDATA", &lucidHitAND_BUNCHDATA_lucidHitAND_Channel, &b_lucidHitAND_BUNCHDATA);
  fChain->SetBranchAddress("lucidHitAND_BCID", lucidHitAND_BCID, &b_lucidHitAND_BCID);
  fChain->SetBranchAddress("lucidHitAND_BunchRawInstLum", lucidHitAND_BunchRawInstLum, &b_lucidHitAND_BunchRawInstLum);
  fChain->SetBranchAddress("lucidEvtA_BUNCHDATA", &lucidEvtA_BUNCHDATA_lucidEvtA_Channel, &b_lucidEvtA_BUNCHDATA);
  fChain->SetBranchAddress("lucidEvtA_BCID", lucidEvtA_BCID, &b_lucidEvtA_BCID);
  fChain->SetBranchAddress("lucidEvtA_BunchRawInstLum", lucidEvtA_BunchRawInstLum, &b_lucidEvtA_BunchRawInstLum);
  fChain->SetBranchAddress("lucidEvtC_BUNCHDATA", &lucidEvtC_BUNCHDATA_lucidEvtC_Channel, &b_lucidEvtC_BUNCHDATA);
  fChain->SetBranchAddress("lucidEvtC_BCID", lucidEvtC_BCID, &b_lucidEvtC_BCID);
  fChain->SetBranchAddress("lucidEvtC_BunchRawInstLum", lucidEvtC_BunchRawInstLum, &b_lucidEvtC_BunchRawInstLum);
  fChain->SetBranchAddress("bcmHXORC_BUNCHDATA", &bcmHXORC_BUNCHDATA_bcmHXORC_Channel, &b_bcmHXORC_BUNCHDATA);
  fChain->SetBranchAddress("bcmHXORC_BCID", bcmHXORC_BCID, &b_bcmHXORC_BCID);
  fChain->SetBranchAddress("bcmHXORC_BunchRawInstLum", bcmHXORC_BunchRawInstLum, &b_bcmHXORC_BunchRawInstLum);
  fChain->SetBranchAddress("bcmVEvtOR_BUNCHDATA", &bcmVEvtOR_BUNCHDATA_bcmVEvtOR_Channel, &b_bcmVEvtOR_BUNCHDATA);
  fChain->SetBranchAddress("bcmVEvtOR_BCID", bcmVEvtOR_BCID, &b_bcmVEvtOR_BCID);
  fChain->SetBranchAddress("bcmVEvtOR_BunchRawInstLum", bcmVEvtOR_BunchRawInstLum, &b_bcmVEvtOR_BunchRawInstLum);
  fChain->SetBranchAddress("bcmVEvtAND_BUNCHDATA", &bcmVEvtAND_BUNCHDATA_bcmVEvtAND_Channel, &b_bcmVEvtAND_BUNCHDATA);
  fChain->SetBranchAddress("bcmVEvtAND_BCID", bcmVEvtAND_BCID, &b_bcmVEvtAND_BCID);
  fChain->SetBranchAddress("bcmVEvtAND_BunchRawInstLum", bcmVEvtAND_BunchRawInstLum, &b_bcmVEvtAND_BunchRawInstLum);
  fChain->SetBranchAddress("bcmVXORC_BUNCHDATA", &bcmVXORC_BUNCHDATA_bcmVXORC_Channel, &b_bcmVXORC_BUNCHDATA);
  fChain->SetBranchAddress("bcmVXORC_BCID", bcmVXORC_BCID, &b_bcmVXORC_BCID);
  fChain->SetBranchAddress("bcmVXORC_BunchRawInstLum", bcmVXORC_BunchRawInstLum, &b_bcmVXORC_BunchRawInstLum);
  fChain->SetBranchAddress("bcmHEvtOR_BUNCHDATA", &bcmHEvtOR_BUNCHDATA_bcmHEvtOR_Channel, &b_bcmHEvtOR_BUNCHDATA);
  fChain->SetBranchAddress("bcmHEvtOR_BCID", bcmHEvtOR_BCID, &b_bcmHEvtOR_BCID);
  fChain->SetBranchAddress("bcmHEvtOR_BunchRawInstLum", bcmHEvtOR_BunchRawInstLum, &b_bcmHEvtOR_BunchRawInstLum);
  fChain->SetBranchAddress("bcmHEvtAND_BUNCHDATA", &bcmHEvtAND_BUNCHDATA_bcmHEvtAND_Channel, &b_bcmHEvtAND_BUNCHDATA);
  fChain->SetBranchAddress("bcmHEvtAND_BCID", bcmHEvtAND_BCID, &b_bcmHEvtAND_BCID);
  fChain->SetBranchAddress("bcmHEvtAND_BunchRawInstLum", bcmHEvtAND_BunchRawInstLum, &b_bcmHEvtAND_BunchRawInstLum);
  fChain->SetBranchAddress("pref_LUMI", &pref_LUMI_pref_LumiChannel, &b_pref_LUMI);
  fChain->SetBranchAddress("lucidEvtAND_LUMI", &lucidEvtAND_LUMI_lucidEvtAND_LumiChannel, &b_lucidEvtAND_LUMI);
  fChain->SetBranchAddress("lucidEvtOR_LUMI", &lucidEvtOR_LUMI_lucidEvtOR_LumiChannel, &b_lucidEvtOR_LUMI);
  fChain->SetBranchAddress("lucidHitOR_LUMI", &lucidHitOR_LUMI_lucidHitOR_LumiChannel, &b_lucidHitOR_LUMI);
  fChain->SetBranchAddress("lucidHitAND_LUMI", &lucidHitAND_LUMI_lucidHitAND_LumiChannel, &b_lucidHitAND_LUMI);
  fChain->SetBranchAddress("lucidEvtA_LUMI", &lucidEvtA_LUMI_lucidEvtA_LumiChannel, &b_lucidEvtA_LUMI);
  fChain->SetBranchAddress("lucidEvtC_LUMI", &lucidEvtC_LUMI_lucidEvtC_LumiChannel, &b_lucidEvtC_LUMI);
  fChain->SetBranchAddress("bcmHXORC_LUMI", &bcmHXORC_LUMI_bcmHXORC_LumiChannel, &b_bcmHXORC_LUMI);
  fChain->SetBranchAddress("mbtsEvtOR_LUMI", &mbtsEvtOR_LUMI_mbtsEvtOR_LumiChannel, &b_mbtsEvtOR_LUMI);
  fChain->SetBranchAddress("mbtsEvtAND_LUMI", &mbtsEvtAND_LUMI_mbtsEvtAND_LumiChannel, &b_mbtsEvtAND_LUMI);
  fChain->SetBranchAddress("mbtsHitOR_LUMI", &mbtsHitOR_LUMI_mbtsHitOR_LumiChannel, &b_mbtsHitOR_LUMI);
  fChain->SetBranchAddress("zdcEvtAND_LUMI", &zdcEvtAND_LUMI_zdcEvtAND_LumiChannel, &b_zdcEvtAND_LUMI);
  fChain->SetBranchAddress("zdcEvtORA_LUMI", &zdcEvtORA_LUMI_zdcEvtORA_LumiChannel, &b_zdcEvtORA_LUMI);
  fChain->SetBranchAddress("bcmVEvtOR_LUMI", &bcmVEvtOR_LUMI_bcmVEvtOR_LumiChannel, &b_bcmVEvtOR_LUMI);
  fChain->SetBranchAddress("bcmVEvtAND_LUMI", &bcmVEvtAND_LUMI_bcmVEvtAND_LumiChannel, &b_bcmVEvtAND_LUMI);
  fChain->SetBranchAddress("bcmVXORC_LUMI", &bcmVXORC_LUMI_bcmVXORC_LumiChannel, &b_bcmVXORC_LUMI);
  fChain->SetBranchAddress("bcmHEvtOR_LUMI", &bcmHEvtOR_LUMI_bcmHEvtOR_LumiChannel, &b_bcmHEvtOR_LUMI);
  fChain->SetBranchAddress("zdcEvtORC_LUMI", &zdcEvtORC_LUMI_zdcEvtORC_LumiChannel, &b_zdcEvtORC_LUMI);
  fChain->SetBranchAddress("bcmHEvtAND_LUMI", &bcmHEvtAND_LUMI_bcmHEvtAND_LumiChannel, &b_bcmHEvtAND_LUMI);
  fChain->SetBranchAddress("lucidEvtAND_TURN", &lucidEvtAND_TURN_lucidEvtAND_BXPhys, &b_lucidEvtAND_TURN);
  fChain->SetBranchAddress("lucidEvtOR_TURN", &lucidEvtOR_TURN_lucidEvtOR_BXPhys, &b_lucidEvtOR_TURN);
  fChain->SetBranchAddress("lucidHitOR_TURN", &lucidHitOR_TURN_lucidHitOR_BXPhys, &b_lucidHitOR_TURN);
  fChain->SetBranchAddress("lucidHitAND_TURN", &lucidHitAND_TURN_lucidHitAND_BXPhys, &b_lucidHitAND_TURN);
  fChain->SetBranchAddress("lucidEvtA_TURN", &lucidEvtA_TURN_lucidEvtA_BXPhys, &b_lucidEvtA_TURN);
  fChain->SetBranchAddress("lucidEvtC_TURN", &lucidEvtC_TURN_lucidEvtC_BXPhys, &b_lucidEvtC_TURN);
  fChain->SetBranchAddress("bcmHXORC_TURN", &bcmHXORC_TURN_bcmHXORC_BXPhys, &b_bcmHXORC_TURN);
  fChain->SetBranchAddress("mbtsEvtOR_TURN", &mbtsEvtOR_TURN_mbtsEvtOR_BXPhys, &b_mbtsEvtOR_TURN);
  fChain->SetBranchAddress("mbtsEvtAND_TURN", &mbtsEvtAND_TURN_mbtsEvtAND_BXPhys, &b_mbtsEvtAND_TURN);
  fChain->SetBranchAddress("mbtsHitOR_TURN", &mbtsHitOR_TURN_mbtsHitOR_BXPhys, &b_mbtsHitOR_TURN);
  fChain->SetBranchAddress("zdcEvtAND_TURN", &zdcEvtAND_TURN_zdcEvtAND_BXPhys, &b_zdcEvtAND_TURN);
  fChain->SetBranchAddress("zdcEvtORA_TURN", &zdcEvtORA_TURN_zdcEvtORA_BXPhys, &b_zdcEvtORA_TURN);
  fChain->SetBranchAddress("bcmVEvtOR_TURN", &bcmVEvtOR_TURN_bcmVEvtOR_BXPhys, &b_bcmVEvtOR_TURN);
  fChain->SetBranchAddress("bcmVEvtAND_TURN", &bcmVEvtAND_TURN_bcmVEvtAND_BXPhys, &b_bcmVEvtAND_TURN);
  fChain->SetBranchAddress("bcmVXORC_TURN", &bcmVXORC_TURN_bcmVXORC_BXPhys, &b_bcmVXORC_TURN);
  fChain->SetBranchAddress("bcmHEvtOR_TURN", &bcmHEvtOR_TURN_bcmHEvtOR_BXPhys, &b_bcmHEvtOR_TURN);
  fChain->SetBranchAddress("zdcEvtORC_TURN", &zdcEvtORC_TURN_zdcEvtORC_BXPhys, &b_zdcEvtORC_TURN);
  fChain->SetBranchAddress("bcmHEvtAND_TURN", &bcmHEvtAND_TURN_bcmHEvtAND_BXPhys, &b_bcmHEvtAND_TURN);
  Notify();
}

Bool_t vdMScanData::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void vdMScanData::Show(Long64_t entry) {
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) {
    return;
  }
  fChain->Show(entry);
}
Int_t vdMScanData::Cut(Long64_t entry) {
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef vdMScanData_cxx
