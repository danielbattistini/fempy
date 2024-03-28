#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//#include "yaml-cpp/yaml.h"

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

void SimLambda1520(int nEvents=200000000, int seed=42, 
                   std::string kinemFilePath="/home/mdicostanzo/an/LPi/Simulation/outputs/SimLambda1520Kinem_MERGED.root", // kinematic information
                   std::string MCFilePath="/home/mdicostanzo/an/LPi/Trains/02_allpc/mc/data/AnalysisResultsAllPC.root",    // momentum smearing
                   std::string outFilePath= "/home/mdicostanzo/an/LPi/Simulation/outputs/SimLambda1520_KinemXi1530_new") {
    
    // map to fill the histograms
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

    // decay mode of interest, from PYTHIA config file

    // <particle id="3124" name="Lambda(1520)0" antiName="Lambda(1520)bar0" spinType="4" chargeType="0" colType="0" 
    //           m0="1.51950" mWidth="0.01560" mMin="1.40000" mMax="1.65000">
    //  <channel onMode="1" bRatio="0.0667000" products="3122 211 -211"/>
    //  <channel onMode="1" bRatio="0.0020000" products="3212 211 -211"/>
    // </particle>

    // switch off lambda(1520) decays and switch on only the relevant ones
    pythia.readString("3124:onMode = off");
    pythia.readString("3122:onMode = off");
    pythia.readString("211:onMode = off");
    pythia.readString("3124:onIfMatch = 3122 211 -211");
    pythia.readString("3124:onIfMatch = 3212 211 -211");
    pythia.readString("3212:onIfMatch = 3122 111");
    pythia.readString("3212:onIfMatch = 3122 22");

    // init
    pythia.init();

    // define vectors to store particles for every event
    std::vector<std::tuple<int, int, int, int>> pions;
    std::vector<std::tuple<int, int, int, int>> lambdas;

    int pdgLambda1520 = 3124;
    double massLambda1520 = TDatabasePDG::Instance()->GetParticle(pdgLambda1520)->Mass();
    // Lifetime of TDatabasePDG in s, in Pythia the unity [mm/c] are used, so a conversion factor of 3*10^11 is needed
    double lifetimeLambda1520 = 3*pow(10, 11)*TDatabasePDG::Instance()->GetParticle(pdgLambda1520)->Lifetime(); 
    
    auto fDecay = new TF1("fDecay", "exp(-x/[0] )", 0., pow(10,-10)); // mm/c
    fDecay->FixParameter(0, lifetimeLambda1520);
    fDecay->SetNpx(1500); 
  
    // Read histograms from Xi(1530) to simulate acceptance and pT of the particle
    TFile *kinemXi1530 = TFile::Open(kinemFilePath.data(), "r");
    TH1F * hAcc1530 = (TH1F*)kinemXi1530->Get("hAcc");
    TH1F * hPt1530 = (TH1F*)kinemXi1530->Get("hPt");

    // Lambda 1520
    Particle lambda1520Feature;
    lambda1520Feature.id(pdgLambda1520);
    lambda1520Feature.m(massLambda1520);

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

        if(iEvent%5000 == 0)
          std::cout << "Processing event " << iEvent << endl;
        
        // initialize particle properties of Lambda(1520)
        double tauLambda1520 = fDecay->GetRandom();
        double phiLambda1520 = gRandom->Rndm() * 2 * TMath::Pi();
        double etaLambda1520 = hAcc1530->GetRandom();
        double thetaLambda1520 = 2 * TMath::ATan( TMath::Exp(-etaLambda1520) );
        double ptLambda1520 = hPt1530->GetRandom();
        double pxLambda1520 = ptLambda1520 * TMath::Cos(phiLambda1520);
        double pyLambda1520 = ptLambda1520 * TMath::Sin(phiLambda1520);
        double pzLambda1520 = ptLambda1520 * TMath::Cos(thetaLambda1520) * (1 / TMath::Sin(thetaLambda1520));
        double pLambda1520 = TMath::Sqrt(ptLambda1520 * ptLambda1520 + pzLambda1520 * pzLambda1520);
        double ELambda1520 = TMath::Sqrt(massLambda1520 * massLambda1520 + pLambda1520 * pLambda1520);
        
        Particle lambda1520;
        lambda1520.id(pdgLambda1520);
        lambda1520.status(81);
        lambda1520.m(massLambda1520);
        lambda1520.xProd(0.);
        lambda1520.yProd(0.);
        lambda1520.zProd(0.);
        lambda1520.tProd(0.);
        lambda1520.e(ELambda1520);
        lambda1520.px(pxLambda1520);
        lambda1520.py(pyLambda1520);
        lambda1520.pz(pzLambda1520);
        lambda1520.tau(tauLambda1520);
        
        // remove all particles generated in the event and append the Lambda(1520)
        pythia.event.reset();
        pythia.event.append(lambda1520);
        int idPartLambda = pythia.event[1].id();
        pythia.particleData.mayDecay(idPartLambda, true);
        
        // force the decay of the Lambda
        pythia.moreDecays();


        for(int iPart=2; iPart<pythia.event.size(); iPart++) {
            if(pythia.event[iPart].id()==3122) {
                  Particle lambda = pythia.event[iPart];
                  TLorentzVector momLambda = TLorentzVector(lambda.px(), lambda.py(), lambda.pz(), lambda.e());

                for(int jPart=3; jPart<pythia.event.size(); jPart++) {
                    if(abs(pythia.event[jPart].id())==211) {
                        Particle pion = pythia.event[jPart];
                        TLorentzVector momPion = TLorentzVector(pion.px(), pion.py(), pion.pz(), pion.e());

                        // acceptance and momentum cuts on the decay products
                        if(lambda.pT()>=0.3 && abs(lambda.eta())<=0.8 && pion.pT()>=0.3 && abs(pion.eta())<=0.8) {
                            TString histoKey =  std::to_string(pythia.event[iPart].id()) + std::to_string(pythia.event[jPart].id());
                            float kStar = RelativePairMomentum(momPion, momLambda);
                            hSEPairs[histoKey]->Fill(kStar); 
                        }
                    }
                }
            }
        }

        // repeat the same procedure for the AntiLambda(1520), using same initialization as the Lambda(1520)
        pythia.event.reset();

        // AntiLambda 1520
        Particle antiLambda1520;
        antiLambda1520.id(-pdgLambda1520);
        antiLambda1520.status(81);
        antiLambda1520.m(massLambda1520);
        antiLambda1520.xProd(0.);
        antiLambda1520.yProd(0.);
        antiLambda1520.zProd(0.);
        antiLambda1520.tProd(0.);
        antiLambda1520.e(ELambda1520);
        antiLambda1520.px(pxLambda1520);
        antiLambda1520.py(pyLambda1520);
        antiLambda1520.pz(pzLambda1520);
        antiLambda1520.tau(tauLambda1520);
        
        pythia.event.append(antiLambda1520);
        int idPartAntiLambda = pythia.event[1].id();
        pythia.particleData.mayDecay(idPartAntiLambda, true);
        pythia.moreDecays();

        for(int iPart=2; iPart<pythia.event.size(); iPart++) {
            if(pythia.event[iPart].id()==-3122) {
                Particle lambda = pythia.event[iPart];
                TLorentzVector momLambda = TLorentzVector(lambda.px(), lambda.py(), lambda.pz(), lambda.e());
            
                for(int jPart=3; jPart<pythia.event.size(); jPart++) {  
                    if(abs(pythia.event[jPart].id())==211) {
                      Particle pion = pythia.event[jPart];
                      TLorentzVector momPion = TLorentzVector(pion.px(), pion.py(), pion.pz(), pion.e());

                        if(lambda.pT()>=0.3 && abs(lambda.eta())<=0.8 && pion.pT()>=0.3 && abs(pion.eta())<=0.8) {
                            TString histoKey =  std::to_string(pythia.event[iPart].id()) + std::to_string(pythia.event[jPart].id());
                            float kStar = RelativePairMomentum(momPion, momLambda);
                            hSEPairs[histoKey]->Fill(kStar); 
                        }
                    }
                }
            }
        }
    }

    // load smearing matrices
    std::vector<TH2F*> smearMatrices;
    TFile *MCdata = TFile::Open(MCFilePath.data());
    TDirectoryFile *folder = static_cast<TDirectoryFile*>(MCdata->Get("HMResultsQA1001"));
    TList *toplist = static_cast<TList*>(folder->Get("HMResultsQA1001"));
    TList *QAList = dynamic_cast<TList*>(toplist->FindObject("PairQA"));

    TList *pairList02 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle0_Particle2"));
    TH2F *smearMomentumMatrix02 = static_cast<TH2F*>(pairList02->FindObject("MomentumResolutionSE_Particle0_Particle2"));
    TString smearMatrix02 = "2113122_smear_matr"; 
    smearMatrices.push_back(smearMomentumMatrix02);

    TList *pairList12 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle1_Particle2"));
    TH2F *smearMomentumMatrix12 = static_cast<TH2F*>(pairList12->FindObject("MomentumResolutionSE_Particle1_Particle2"));
    TString smearMatrix12 = "-2113122_smear_matr"; 
    smearMatrices.push_back(smearMomentumMatrix12);

    TList *pairList03 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle0_Particle3"));
    TH2F *smearMomentumMatrix03 = static_cast<TH2F*>(pairList03->FindObject("MomentumResolutionSE_Particle0_Particle3"));
    TString smearMatrix03 = "211-3122_smear_matr"; 
    smearMatrices.push_back(smearMomentumMatrix03);

    TList *pairList13 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle1_Particle3"));
    TH2F *smearMomentumMatrix13 = static_cast<TH2F*>(pairList13->FindObject("MomentumResolutionSE_Particle1_Particle3"));
    TString smearMatrix13 = "-211-3122_smear_matr"; 
    smearMatrices.push_back(smearMomentumMatrix13);

    MCdata->Close();

    // compute the smeared histograms
    for(int iKey=0; iKey<keys.size(); iKey++) {
      for(int iBin=0; iBin<hSEPairsSmeared[keys[iKey]]->GetNbinsX(); iBin++) {
        TH1D *hProjY = smearMatrices[iKey]->ProjectionY(Form("iBin_%i", iBin+1), iBin+1, iBin+1);
        if(hProjY->Integral() > 0) {    
            hProjY->Scale(hSEPairs[keys[iKey]]->GetBinContent(iBin+1)/hProjY->Integral());
            hSEPairsSmeared[keys[iKey]]->Add(hProjY);
        }
      }
      for(int iBin=0; iBin<hSEPairsSmeared[keys[iKey]]->GetNbinsX(); iBin++) {
        hSEPairsSmeared[keys[iKey]]->SetBinError(iBin+1, TMath::Sqrt(hSEPairsSmeared[keys[iKey]]->GetBinContent(iBin+1)));
      }
    }

    // Write histograms to file
    std::string outFileName = outFilePath + "_" + std::to_string(seed) + ".root";  
    TFile oFile(outFileName.data(), "recreate");
    oFile.cd();
    for(const auto &hSEPair : hSEPairs) {
        hSEPairs[hSEPair.first]->Write();
        hSEPairsSmeared[hSEPair.first]->Write();
    }

    // Resum pair and antipair
    TH1F * resumLPiPl = new TH1F("3122211_-3122-211", ";#it{k*};Counts", 1500, 0., 6.);
    resumLPiPl->Add(hSEPairs["3122211"]);
    resumLPiPl->Add(hSEPairs["-3122-211"]);
    resumLPiPl->Write();
    TH1F * resumLPiMin = new TH1F("3122-211_-3122211", ";#it{k*};Counts", 1500, 0., 6.);
    resumLPiMin->Add(hSEPairs["3122-211"]);
    resumLPiMin->Add(hSEPairs["3122-211"]);
    resumLPiMin->Write();
    TH1F * resumLPiPlSmear = new TH1F("3122211_-3122-211_smeared", ";#it{k*};Counts", 1500, 0., 6.);
    resumLPiPlSmear->Add(hSEPairs["3122211"]);
    resumLPiPlSmear->Add(hSEPairs["-3122-211"]);
    resumLPiPlSmear->Write();
    TH1F * resumLPiMinSmear = new TH1F("3122-211_-3122211_smeared", ";#it{k*};Counts", 1500, 0., 6.);
    resumLPiMinSmear->Add(hSEPairs["3122-211"]);
    resumLPiMinSmear->Add(hSEPairs["3122-211"]);
    resumLPiMinSmear->Write();

    hAcc1530->Write();
    hPt1530->Write();
    fDecay->Write();
    oFile.Close();
 
}