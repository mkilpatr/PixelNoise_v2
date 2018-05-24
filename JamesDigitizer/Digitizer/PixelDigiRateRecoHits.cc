// system include files
#include <memory>
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

// To use root histos
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" 
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

// data formats
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

// For ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TTree.h>

using namespace std;
using namespace edm;

class PixelDigiRateRecoHits : public EDAnalyzer {
  public:
    explicit PixelDigiRateRecoHits(const ParameterSet&);
    ~PixelDigiRateRecoHits();
    virtual void beginJob(); // obsolete const EventSetup& iSetup);
    virtual void analyze(const Event&, const EventSetup&);
    virtual void endJob();

  private:
    ParameterSet iConfig_;
    unsigned int numevt;
    TH1D *numEvents;
    TH2F *hpixMap1, *hpixMap2, *hpixMap3, *hpixMap4; 
    TH2F *hpixMapF1[2],*hpixMapF2[2],*hpixMapF3[2];
    TFile *myTFile_;

    struct eventInfo {
      unsigned int aRunNum;
      unsigned int aLumiBlock;
      unsigned int aEvtNum;
    };
    eventInfo aEvent;

    TTree *digiTree;  //full out tree of all digis
    struct myDigi {  //fully specified hit
      int aLayer;  // layer 1-4 barrel, 5-7 -z disks, 8-10 +z disks
      int aSector; // azimuthal index (ladder/sector/blade) 1-64 b / 1-56 f
      int aModule; // z index = module in barrel (1-8); panel/side in fpix (1-2)
      int aRoc;    // readout chip (0-15)
      int aDcol;   // double column (0-25)
      int aPix;    // pixel in double column (0-159)
    };
    myDigi aDigi;

    TTree *recHitTree;
    struct myRecHit {
      int aLayer;  // layer 1-4 barrel, 5-7 -z disks, 8-10 +z disks
      int aSector; // azimuthal index (ladder/sector/blade) 1-64 b / 1-56 f
      int aModule; // z index = module in barrel (1-8); panel/side in fpix (1-2)
//      int aDcol;   // double column (0-25)
//      int aPix;    // pixel in double column (0-159)
      float x;
      float y;
//      float xx;
//      float xy;
//      float yy;
      float gx;
      float gy;
      float gz;
      int nsimhit;
      float hx, hy;
    };
    myRecHit aRecHit_;

    InputTag src1_;
    InputTag src2_;
};

PixelDigiRateRecoHits::PixelDigiRateRecoHits(const ParameterSet& iConfig) {
  iConfig_ = iConfig;
  src1_ = iConfig_.getParameter <InputTag> ("src1");
  src2_ = iConfig_.getParameter <InputTag> ("src2");
  numevt = 0;
  cout << " Construct PixelDigiRateRecoHits " << endl;
}

PixelDigiRateRecoHits::~PixelDigiRateRecoHits() {
  cout << " Destroy PixelDigiRateRecoHits " << endl;
}

void PixelDigiRateRecoHits::beginJob() {
  cout << "Initialize PixelDigiRateRecoHits " <<endl;
  string outputFile = iConfig_.getParameter<string>("OutputFile");
  myTFile_ = new TFile (outputFile.c_str() , "RECREATE");
  digiTree = new TTree("pixdigis", "SiPixelDigis");
  digiTree->Branch("eventInfo", &aEvent, "runNum/i:lumiBlock:evtNum");
  digiTree->Branch("digi", &aDigi, "layer/I:sector:module:roc:dcol:pix");
  recHitTree = new TTree("WhatExactlyAmIRecording", "YouWillNeedToEditThisLater");
  recHitTree->Branch("eventInfo", &aEvent, "runNum/i:lumiBlock:evtNum");
  recHitTree->Branch("recHit", &aRecHit_, "layer/I:sector:module:x/F:y:gx:gy:gz:nsimhit/I:hx/F:hy");
  numEvents = new TH1D("numEvents"," ",1,0,1);
  hpixMap1 = new TH2F("hpixMap1"," ",416,0.,416.,160,0.,160.);
  hpixMap1->SetOption("colz");
  hpixMap2 = new TH2F("hpixMap2"," ",416,0.,416.,160,0.,160.);
  hpixMap2->SetOption("colz");
  hpixMap3 = new TH2F("hpixMap3"," ",416,0.,416.,160,0.,160.);
  hpixMap3->SetOption("colz");
  hpixMap4 = new TH2F("hpixMap4"," ",416,0.,416.,160,0.,160.);
  hpixMap4->SetOption("colz");
  hpixMapF1[0] = new TH2F("hpixMapF1in","Fwd Disk1 inner",416,0.,416.,160,0.,160.);
  hpixMapF1[0]->SetOption("colz");
  hpixMapF1[1] = new TH2F("hpixMapF1out","Fwd Disk1 outer",416,0.,416.,160,0.,160.);
  hpixMapF1[1]->SetOption("colz");
  hpixMapF2[0] = new TH2F("hpixMapF2in","Fwd Disk2 inner",416,0.,416.,160,0.,160.);
  hpixMapF2[0]->SetOption("colz");
  hpixMapF2[1] = new TH2F("hpixMapF2out","Fwd Disk2 outer",416,0.,416.,160,0.,160.);
  hpixMapF2[1]->SetOption("colz");
  hpixMapF3[0] = new TH2F("hpixMapF3in","Fwd Disk3 inner",416,0.,416.,160,0.,160.);
  hpixMapF3[0]->SetOption("colz");
  hpixMapF3[1] = new TH2F("hpixMapF3out","Fwd Disk3 outer",416,0.,416.,160,0.,160.);
  hpixMapF3[1]->SetOption("colz");
}

void PixelDigiRateRecoHits::analyze(const Event& iEvent, const EventSetup& iSetup) {
  numEvents->Fill(0);
// Get digis
  Handle <DetSetVector <PixelDigi> > pixelDigis;
  iEvent.getByLabel(src1_, pixelDigis);
// Get event setup (to get global transformation)
  ESHandle <TrackerGeometry> geom;
  iSetup.get<TrackerDigiGeometryRecord>().get(geom);
  const TrackerGeometry& theTracker(*geom);

  aEvent.aRunNum = iEvent.id().run();
  aEvent.aLumiBlock = iEvent.id().luminosityBlock();
  aEvent.aEvtNum = iEvent.id().event();

// Iterate on detector units
  bool thisOneHadHits = false;
  DetSetVector <PixelDigi>::const_iterator DSViter;
  if(pixelDigis->size() != 0){
    numevt++;
    thisOneHadHits = true;
  } else {
    cout << "run " << iEvent.id().run() << ", lumi block " << iEvent.id().luminosityBlock() << ", event " << iEvent.id().event() << " did not have any pixel digis!!!" << endl;
  }
  for(DSViter = pixelDigis->begin(); DSViter != pixelDigis->end(); DSViter++) {
    unsigned int detid = DSViter->id; // = rawid
    DetId detId(detid);
    unsigned int detType = detId.det(); // det type, tracker=1
    unsigned int subid = detId.subdetId(); //subdetector type, barrel=1
    if(detType != 1) continue; // look only at tracker
// Get the geom-detector 
    const PixelGeomDetUnit *theGeomDet = dynamic_cast <const PixelGeomDetUnit*> (theTracker.idToDet(detId));
    double detR = theGeomDet->surface().position().perp();
    unsigned int layer = 0, ladder = 0, zindex = 0, disk = 0; //1,2,3
    int sdisk = 0; // -3,-2,-1, 1,2,3
    unsigned int blade = 0; //1-24
    unsigned int side = 0; //size=1 for -z, 2 for +z
    unsigned int panel = 0; //panel=1
    unsigned int inout = 0; // innerdisk =0, outerdisk=1 
// Subdet it, pix barrel=1, forward=2
    if(subid == 2) {   // forward
      if(detR > 10.) {
        inout = 1; //set for outerdisk
      }
      PXFDetId pdetId = PXFDetId(detid);
      disk = pdetId.disk(); //1,2,3
      blade = pdetId.blade(); //1-56  (1-22 inner; 23-36 outer)
      zindex = pdetId.module(); //1 ; not used
      side = pdetId.side(); //size=1 for -z, 2 for +z
      sdisk = disk + 4;  //layer = 5,6,7 for disk 1,2,3 (-z) 
      if(side == 2) {
        sdisk += 3;// layer=8,9,10 for disk 1,2,3 (+z)
      }
      panel = pdetId.panel(); //panel=1,2
    } else if(subid == 1) { // Barrel 
      PXBDetId pdetId = PXBDetId(detid);
      layer = pdetId.layer();  // Barell layer = 1,2,3,4
      ladder = pdetId.ladder();  // Barrel ladder id 1-20,32,44.
      zindex = pdetId.module();  // Barrel Z-index=1,8
    } // end fb-bar

// Look at digis now
    DetSet<PixelDigi>::const_iterator di;
    for(di = DSViter->data.begin(); di != DSViter->data.end(); di++) {
      int col = di->column(); // column 
      int row = di->row();    // row
// warning: detailed geometry not considered carefully, e.g. (ROC ID,col,row) is unique but might not match reality
      if(layer == 1) {
        hpixMap1->Fill(float(col),float(row));
      } else if(layer == 2) {
        hpixMap2->Fill(float(col),float(row));
      } else if(layer == 3) {
        hpixMap3->Fill(float(col),float(row));
      } else if(layer == 4){
        hpixMap4->Fill(float(col),float(row));
      } else if(disk == 1) {
        hpixMapF1[inout]->Fill(float(col),float(row));
      } else if(disk == 2) {
        hpixMapF2[inout]->Fill(float(col),float(row));
      } else if(disk == 3) {
        hpixMapF3[inout]->Fill(float(col),float(row));
      }
      if(1 == subid) {
        aDigi.aLayer = layer;
        aDigi.aSector = ladder;
        aDigi.aModule = zindex;
      }
      else {
        aDigi.aLayer = sdisk;
        aDigi.aSector = blade;
        aDigi.aModule = panel;
      }
//      aDigi.aRoc = col/52 + 8*(row/80);
      aDigi.aDcol = col;
      aDigi.aPix = row;
      digiTree->Fill();
      if(!thisOneHadHits) {
        cout << "even though i didn't have pixel digis, i still found something to fill in the digiTree!!  wtf?!?!?!" << endl;
      }
    } // end for digis
  } // end for det-units

  Handle <SiPixelRecHitCollection> recHitColl;
  iEvent.getByLabel(src2_, recHitColl);
// for finding matched simhit
  TrackerHitAssociator associate(iEvent, iConfig_);
  if((recHitColl.product())->dataSize() > 0) {
    vector<PSimHit> matched;
    vector<PSimHit>::const_iterator closest_simhit;
// Loop over Detector IDs
    for(SiPixelRecHitCollection::const_iterator recHitIdIterator = recHitColl.product()->begin(); recHitIdIterator != recHitColl.product()->end(); recHitIdIterator++) {
      SiPixelRecHitCollection::DetSet detset = *recHitIdIterator;
      if(detset.empty()) continue;
      DetId detId = DetId(detset.detId()); // Get the Detid object
      unsigned int detType = detId.det(); // det type, tracker=1
      if(detType != 1) continue; // look only at tracker
      const PixelGeomDetUnit *theGeomDet = dynamic_cast <const PixelGeomDetUnit*> (theTracker.idToDet(detId));
      //const GeomDet* geomDet(theTracker->idToDet(detId));
// Loop over rechits for this detid
      for(SiPixelRecHitCollection::DetSet::const_iterator iterRecHit = detset.begin(); iterRecHit != detset.end(); ++iterRecHit) {
        // get matched simhit
        matched.clear();
        matched = associate.associateHit(*iterRecHit);
        if(!matched.empty()) {
          float closest = 9999.9;
          vector<PSimHit>::const_iterator closestit = matched.begin();
          LocalPoint lp = iterRecHit->localPosition();
          float rechit_x = lp.x();
          float rechit_y = lp.y();
//loop over simhits and find closest
          for(vector<PSimHit>::const_iterator m = matched.begin(); m < matched.end(); m++) {
            float sim_x1 = (*m).entryPoint().x();
            float sim_x2 = (*m).exitPoint().x();
            float sim_xpos = 0.5 * (sim_x1 + sim_x2);
            float sim_y1 = (*m).entryPoint().y();
            float sim_y2 = (*m).exitPoint().y();
            float sim_ypos = 0.5 * (sim_y1 + sim_y2);
            float x_res = fabs(sim_xpos - rechit_x);
            float y_res = fabs(sim_ypos - rechit_y);
            float dist = sqrt(x_res * x_res + y_res * y_res);
            if(dist < closest) {
              closest = dist;
              closestit = m;
            }
          } // end of simhit loop
          closest_simhit = closestit;
        } // end matched emtpy
        unsigned int subid = detId.subdetId(); //subdetector type, barrel=1
        int sdisk = 0; // -3,-2,-1, 1,2,3
        unsigned int disk = 0; //1,2,3
        unsigned int side = 0; //size=1 for -z, 2 for +z
        // Subdet it, pix barrel=1, forward=2
        if(subid == 1) { // Barrel 
          PXBDetId pdetId = PXBDetId(detset.detId());
          aRecHit_.aLayer = pdetId.layer();  // Barell layer = 1,2,3,4
          aRecHit_.aSector = pdetId.ladder();  // Barrel ladder id 1-20,32,44.
          aRecHit_.aModule = pdetId.module();  // Barrel Z-index=1,8
        } else if(subid == 2) {   // forward
          PXFDetId pdetId = PXFDetId(detset.detId());
          disk = pdetId.disk(); //1,2,3
          aRecHit_.aSector = pdetId.blade(); //1-56  (1-22 inner; 23-36 outer)
          side = pdetId.side(); //size=1 for -z, 2 for +z
          sdisk = disk + 4;  //layer = 5,6,7 for disk 1,2,3 (-z) 
          if(side == 2) {
            sdisk += 3;// layer=8,9,10 for disk 1,2,3 (+z)
          }
          aRecHit_.aModule = pdetId.panel(); //panel=1,2
          aRecHit_.aLayer = sdisk;
        }
        LocalPoint lp = iterRecHit->localPosition();
        aRecHit_.x = lp.x();
        aRecHit_.y = lp.y();
        //LocalError le = iterRecHit->localPositionError();
        //aRecHit_.xx = le.xx();
        //aRecHit_.xy = le.xy();
        //aRecHit_.yy = le.yy();
        GlobalPoint GP = theGeomDet->surface().toGlobal(iterRecHit->localPosition());
        aRecHit_.gx = GP.x();
        aRecHit_.gy = GP.y();
        aRecHit_.gz = GP.z();
        aRecHit_.nsimhit = matched.size();
        if(matched.size() > 0) {
          float sim_x1 = (*closest_simhit).entryPoint().x();
          float sim_x2 = (*closest_simhit).exitPoint().x();
          aRecHit_.hx = 0.5 * (sim_x1 + sim_x2);
          float sim_y1 = (*closest_simhit).entryPoint().y();
          float sim_y2 = (*closest_simhit).exitPoint().y();
          aRecHit_.hy = 0.5 * (sim_y1 + sim_y2);
        }
        recHitTree->Fill();
        if(!thisOneHadHits) {
          cout << "even though i didn't have pixel digis, i still found something to fill in the recHitTree!!  wtf?!?!?!" << endl;
        }
      } // end of rechit loop
    } // end of detid loop
  } // end of loop test on recHitColl size
}

void PixelDigiRateRecoHits::endJob(){
  cout << " PROCESSED EVENTS:  " << numevt << endl;
  cout << " End PixelDigiRateRecoHits " << endl;
  myTFile_->Write();
  myTFile_->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(PixelDigiRateRecoHits);
