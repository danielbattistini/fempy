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

void SimMisidLambdas(int nEvents=2000, int seed=42) {
    
    std::map<TString, TH1F*> hSEPairs;
    std::map<TString, TH1F*> hSEPairsSmeared;
    std::vector<TString> keys = {"3122211", "3122-211", "-3122211", "-3122-211"};
    for(int iKey=0; iKey<keys.size(); iKey++) {
      TString smearHistoName = keys[iKey] + "_smeared";
      hSEPairs.insert({keys[iKey], new TH1F(keys[iKey].Data(), ";#it{k*};Counts", 1500, 0., 6.)});
      hSEPairsSmeared.insert({keys[iKey], new TH1F(smearHistoName.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
    }
    
    Pythia pythia;
    
    // set seed for simulation
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:all = on");

    pythia.readString("3122:onMode = off");
    pythia.readString("211:onMode = off");
    pythia.readString("3122:onIfMatch = 2212 -211");

    // init
    pythia.init();

    // perform the simulation

    // define histograms
    std::map<TString, TH1F*> hSEMothersPairs;
    std::vector<TString> histoKeys;
    
    // define vectors to store particles for every event
    std::vector<Particle> pions;
    std::vector<Particle> pionsFromLambdas;
    std::vector<Particle> protons;

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
      if(iEvent%10000 == 0) {
            std::cout << Form("Processing event... ", iEvent) << endl;
        }
        if (!pythia.next()) {
            continue;
        } 
        for (int iPart = 2; iPart < pythia.event.size(); iPart++) {
            Particle part = pythia.event.at(iPart);
            if(abs(part.id()) == 211) {
                int motherIdx = part.mother1();
                Particle pionMother = pythia.event.at(motherIdx);
                if(abs(pionMother.id()) == 3122) {
                    int daugProton = pionMother.daughter1();
                    Particle lambdaProton = pythia.event.at(daugProton);
                    int daugPion = pionMother.daughter2();
                    Particle lambdaPion = pythia.event.at(daugPion);
                    cout << "------------------------" << endl;
                    cout << "Lambda mother id: " << pionMother.id() << endl;
                    cout << "Proton daughter of lambda: " << lambdaProton.id() << endl;
                    cout << "Pion daughter of lambda: " << lambdaPion.id() << endl;
                    cout << "Stack idxs: " << iPart << " " << daugPion << endl;
                    cout << "------------------------" << endl;
                } else {
                    cout << "------------------------" << endl;
                    cout << "Particle id: " << part.id << endl;
                    cout << "Mother id: " << pionMother.id() << endl;
                    cout << "------------------------" << endl;
                    pions.push_back(part);
                }
            }
        }    
    // for(int iPart=2; iPart<pythia.event.size(); iPart++) {
    //     if(pythia.event[iPart].id()==3122) {
    //       Particle lambda = pythia.event[iPart];
    //       TLorentzVector momLambda = TLorentzVector(lambda.px(), lambda.py(), lambda.pz(), lambda.e());
          
    //       for(int jPart=3; jPart<pythia.event.size(); jPart++) {
    //         if(abs(pythia.event[jPart].id())==211) {
    //           Particle pion = pythia.event[jPart];
    //           TLorentzVector momPion = TLorentzVector(pion.px(), pion.py(), pion.pz(), pion.e());
          
    //           if(lambda.pT()>=0.3 && abs(lambda.eta())<=0.8 && pion.pT()>=0.3 && abs(pion.eta())<=0.8) {
    //             TString histoKey =  std::to_string(pythia.event[iPart].id()) + std::to_string(pythia.event[jPart].id());
    //             float kStar = RelativePairMomentum(momPion, momLambda);
    //             hSEPairs[histoKey]->Fill(kStar); 
    //           }
          
    //         }
    //       }
        
    //     }
    //   }
    }

    // std::vector<TH2F*> smearMatrices;
    // cout << "CIAO1" << endl;
    // TFile *MCdata = TFile::Open("/home/mdicostanzo/an/LPi/Trains/02_allpc/mc/data/AnalysisResultsAllPC.root");
    // cout << "CIAO2" << endl;
    // TDirectoryFile *folder = static_cast<TDirectoryFile*>(MCdata->Get("HMResultsQA1001"));
    // TList *toplist = static_cast<TList*>(folder->Get("HMResultsQA1001"));
    // TList *QAList = dynamic_cast<TList*>(toplist->FindObject("PairQA"));

    // TList *pairList02 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle0_Particle2"));
    // TH2F *smearMomentumMatrix02 = static_cast<TH2F*>(pairList02->FindObject("MomentumResolutionSE_Particle0_Particle2"));
    // TString smearMatrix02 = "2113122_smear_matr"; 
    // smearMatrices.push_back(smearMomentumMatrix02);

    // TList *pairList12 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle1_Particle2"));
    // TH2F *smearMomentumMatrix12 = static_cast<TH2F*>(pairList12->FindObject("MomentumResolutionSE_Particle1_Particle2"));
    // TString smearMatrix12 = "-2113122_smear_matr"; 
    // smearMatrices.push_back(smearMomentumMatrix12);

    // TList *pairList03 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle0_Particle3"));
    // TH2F *smearMomentumMatrix03 = static_cast<TH2F*>(pairList03->FindObject("MomentumResolutionSE_Particle0_Particle3"));
    // TString smearMatrix03 = "211-3122_smear_matr"; 
    // smearMatrices.push_back(smearMomentumMatrix03);

    // TList *pairList13 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle1_Particle3"));
    // TH2F *smearMomentumMatrix13 = static_cast<TH2F*>(pairList13->FindObject("MomentumResolutionSE_Particle1_Particle3"));
    // TString smearMatrix13 = "-211-3122_smear_matr"; 
    // smearMatrices.push_back(smearMomentumMatrix13);
    // cout << "CIAO3" << endl;

    // MCdata->Close();

    // for(int iKey=0; iKey<keys.size(); iKey++) {
    //   for(int iBin=0; iBin<hSEPairsSmeared[keys[iKey]]->GetNbinsX(); iBin++) {
    //     TH1D *hProjY = smearMatrices[iKey]->ProjectionY(Form("iBin_%i", iBin+1), iBin+1, iBin+1);
    //     if(hProjY->Integral() > 0) {    
    //         hProjY->Scale(hSEPairs[keys[iKey]]->GetBinContent(iBin+1)/hProjY->Integral());
    //         hSEPairsSmeared[keys[iKey]]->Add(hProjY);
    //     }
    //   }
    //   for(int iBin=0; iBin<hSEPairsSmeared[keys[iKey]]->GetNbinsX(); iBin++) {
    //     hSEPairsSmeared[keys[iKey]]->SetBinError(iBin+1, TMath::Sqrt(hSEPairsSmeared[keys[iKey]]->GetBinContent(iBin+1)));
    //   }
    // }
    // cout << "CIAO3" << endl;

    // cout << "CIAO1" << endl;
    // std::string outFileName = "/home/mdicostanzo/an/LPi/Simulation/outputs/SimLambda1520_KinemXi1530_" + std::to_string(seed) + ".root";  
    // TFile oFile(outFileName.data(), "recreate");
    // oFile.cd();
    // for(const auto &hSEPair : hSEPairs) {
    //     hSEPairs[hSEPair.first]->Write();
    //     hSEPairsSmeared[hSEPair.first]->Write();
    // }

    // TH1F * resumLPiPl = new TH1F("3122211_-3122-211", ";#it{k*};Counts", 1500, 0., 6.);
    // resumLPiPl->Add(hSEPairs["3122211"]);
    // resumLPiPl->Add(hSEPairs["-3122-211"]);
    // resumLPiPl->Write();
    // TH1F * resumLPiMin = new TH1F("3122-211_-3122211", ";#it{k*};Counts", 1500, 0., 6.);
    // resumLPiMin->Add(hSEPairs["3122-211"]);
    // resumLPiMin->Add(hSEPairs["3122-211"]);
    // resumLPiMin->Write();
    // TH1F * resumLPiPlSmear = new TH1F("3122211_-3122-211_smeared", ";#it{k*};Counts", 1500, 0., 6.);
    // resumLPiPlSmear->Add(hSEPairs["3122211"]);
    // resumLPiPlSmear->Add(hSEPairs["-3122-211"]);
    // resumLPiPlSmear->Write();
    // TH1F * resumLPiMinSmear = new TH1F("3122-211_-3122211_smeared", ";#it{k*};Counts", 1500, 0., 6.);
    // resumLPiMinSmear->Add(hSEPairs["3122-211"]);
    // resumLPiMinSmear->Add(hSEPairs["3122-211"]);
    // resumLPiMinSmear->Write();

    // cout << "CIAO2" << endl;
    // hAcc1530->Write();
    // cout << "CIAO3" << endl;
    // hPt1530->Write();
    // cout << "CIAO4" << endl;
    // fDecay->Write();
    // cout << "CIAO5" << endl;
    // oFile.Close();
    // cout << "CIAO6" << endl;

}