#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TLorentzVector.h>

//#include "AliPythia8.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

float RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo) {
  TLorentzVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = PartOne;
  TLorentzVector PartTwoCMS = PartTwo;

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}

void SimLambda1520Kinem(int seed=42, int nEvents=20000) {

    // Script to obtain pT and eta distributions for Lambda(1520) from Xi(1530)
    
    Pythia pythia;
    
    // set seed for simulation
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:all = on");
    pythia.readString("Beams:eCM = 13000");

    // init
    pythia.init();

    // save output
    std::string outFileName = "./SimLambda1520Kinem_" + std::to_string(seed) + ".root";  
    TFile oFile(outFileName.data(), "recreate");
    
    TH1F* hAcc = new TH1F("hAcc", ";#eta;Counts", 2000, -2, 2);
    TH1F* hPt = new TH1F("hPt", ";#p_T;Counts", 3000, 0, 6);

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if(iEvent%10000 == 0) {
            std::cout << Form("Processing event... ", iEvent) << endl;
        }
        if (!pythia.next()) {
            continue;
        } 
        for (int iPart = 2; iPart < pythia.event.size(); iPart++) {
            Particle part = pythia.event.at(iPart);
            if(part.id() == 3314) {
                hAcc->Fill(part.eta());
                hPt->Fill(part.pT());
            }
        }
    }

    oFile.cd();
    hAcc->Write();
    hPt->Write();
    oFile.Close();
}
