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

void SimLambda1520(int nEvents=10000, int seed=42) {
    
    std::map<TString, TH1F*> hSEPairs;
    TString lambdaPiPlus("3122211");
    hSEPairs.insert({lambdaPiPlus, new TH1F(lambdaPiPlus.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
    TString lambdaPiMinus("3122-211");
    hSEPairs.insert({lambdaPiMinus, new TH1F(lambdaPiMinus.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
    TString antiLambdaPiPlus("-3122211");
    hSEPairs.insert({antiLambdaPiPlus, new TH1F(antiLambdaPiPlus.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
    TString antiLambdaPiMinus("-3122-211");
    hSEPairs.insert({antiLambdaPiMinus, new TH1F(antiLambdaPiMinus.Data(), ";#it{k*};Counts", 1500, 0., 6.)});

    Pythia pythia;
    
    // set seed for simulation
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:all = on");

    // keep only interesting decays, to be reweighted a posteriori

    // <particle id="3124" name="Lambda(1520)0" antiName="Lambda(1520)bar0" spinType="4" chargeType="0" colType="0" 
    //           m0="1.51950" mWidth="0.01560" mMin="1.40000" mMax="1.65000">
    //  <channel onMode="1" bRatio="0.0667000" products="3122 211 -211"/>
    //  <channel onMode="1" bRatio="0.0020000" products="3212 211 -211"/>
    // </particle>

    pythia.readString("3124:onMode = off");
    pythia.readString("3122:onMode = off");
    pythia.readString("211:onMode = off");
    pythia.readString("3124:onIfMatch = 3122 211 -211");
    pythia.readString("3124:onIfMatch = 3212 211 -211");
    pythia.readString("3212:onIfMatch = 3122 111");
    pythia.readString("3212:onIfMatch = 3122 22");

    //pythia.readString("3124:tau0=4.41000e-01"); // bRatio="0.0800000"

    //pythia.readString("102134:oneChannel = 1 1 0 3122 211 -211");

    // init
    pythia.init();

    // perform the simulation

    // define histograms
    std::map<TString, TH1F*> hSEMothersPairs;
    std::vector<TString> histoKeys;
    
    // define vectors to store particles for every event
    std::vector<std::tuple<int, int, int, int>> pions;
    std::vector<std::tuple<int, int, int, int>> lambdas;
    //
    int countLambda1520 = 0;
    int pdgLambda1520 = 3124;
    double massLambda1520 = TDatabasePDG::Instance()->GetParticle(pdgLambda1520)->Mass();
    cout << "Lambda1520 mass " << massLambda1520 << endl;
    double hBar = 1.054571817*pow(10,-34);
    cout << "Lambda(1520) lifetime: " << TDatabasePDG::Instance()->GetParticle(pdgLambda1520)->Lifetime() << endl; 
    cout << "Lambda(1520) lifetime multiplied: " << pow(10, 20)*TDatabasePDG::Instance()->GetParticle(pdgLambda1520)->Lifetime() << endl; 
    double lifetimeLambda1520 = 3*pow(10, 11)*TDatabasePDG::Instance()->GetParticle(pdgLambda1520)->Lifetime(); 
    // Lifetime of TDatabasePDG in s, in Pythia the unity [mm/c] are used, so a conversion factor of 3*10^11 is needed
    //double lifetimeLambda1520 = 3*pow(10, 11)*TDatabasePDG::Instance()->GetParticle(pdgLambda1520)->Lifetime(); 
    cout << "Lambda1520 lifetime: " << lifetimeLambda1520 << endl;
    auto fDecay = new TF1("fDecay", "exp(-x/[0] )", 0., pow(10,-10)); // seconds
    fDecay->FixParameter(0, lifetimeLambda1520);
    fDecay->SetNpx(1500); // f1->GetNpx() > f1->GetXmax() / f1->GetXmin()
    auto fPt = new TF1("fPt", "exp(-x/1)", 0., 5.); // GeV
    //auto fDecay = new TF1("fDecay", "exp(-x/100 )", 0., 1000); // seconds

    // Lambda 1520
    Particle lambda1520Feature;
    lambda1520Feature.id(pdgLambda1520);
    lambda1520Feature.m(massLambda1520);
    cout << "Lambda(1520) mass: " << lambda1520Feature.m() << endl;

    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        cout << endl;
        cout << "-------------------------" << endl;
        // if(iEvent%10000 == 0)
        std::cout << "Processing event " << iEvent << endl;
        double tauLambda1520 = fDecay->GetRandom();
        double phiLambda1520 = gRandom->Rndm() * 2 * TMath::Pi();
        double etaLambda1520 = gRandom->Uniform(-2., 2.);
        double thetaLambda1520 = 2 * TMath::ATan( TMath::Exp(-etaLambda1520) );
        double ptLambda1520 = fPt->GetRandom();
        double pxLambda1520 = ptLambda1520 * TMath::Cos(phiLambda1520);
        double pyLambda1520 = ptLambda1520 * TMath::Sin(phiLambda1520);
        double pzLambda1520 = ptLambda1520 * TMath::Cos(thetaLambda1520) * (1 / TMath::Sin(thetaLambda1520));
        double pLambda1520 = TMath::Sqrt(ptLambda1520 * ptLambda1520 + pzLambda1520 * pzLambda1520);
        double ELambda1520 = TMath::Sqrt(massLambda1520 * massLambda1520 + pLambda1520 * pLambda1520);
        
        // Lambda 1520
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
        
        pythia.event.reset();
        pythia.event.append(lambda1520);
        int idPartLambda = pythia.event[1].id();
        cout << idPartLambda << endl;
        pythia.particleData.mayDecay(idPartLambda, true);
        pythia.moreDecays();

        for(int iPart=2; iPart<pythia.event.size(); iPart++) {
          if(pythia.event[iPart].id()==3122) {
            Particle lambda = pythia.event[iPart];
            cout << "Quadrimomentum of Lambda: " << lambda.px() << " " << lambda.py() << " " << lambda.pz() << " " << lambda.e() << " " << endl;
            TLorentzVector momLambda = TLorentzVector(lambda.px(), lambda.py(), lambda.pz(), lambda.e());
            for(int jPart=3; jPart<pythia.event.size(); jPart++) {
              if(abs(pythia.event[jPart].id())==211) {
                Particle pion = pythia.event[jPart];
                cout << "Quadrimomentum of pion: " << pion.px() << " " << pion.py() << " " << pion.pz() << " " << pion.e() << " " << endl;
                TLorentzVector momPion = TLorentzVector(pion.px(), pion.py(), pion.pz(), pion.e());
                TString histoKey =  std::to_string(pythia.event[iPart].id()) + std::to_string(pythia.event[jPart].id());
                float kStar = RelativePairMomentum(momPion, momLambda);
                hSEPairs[histoKey]->Fill(kStar); 
              }
            }
          }
        }

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
        cout << idPartAntiLambda << endl;
        pythia.particleData.mayDecay(idPartAntiLambda, true);
        pythia.moreDecays();

        for(int iPart=2; iPart<pythia.event.size(); iPart++) {
          if(pythia.event[iPart].id()==-3122) {
            Particle lambda = pythia.event[iPart];
            cout << "Quadrimomentum of Lambda: " << lambda.px() << " " << lambda.py() << " " << lambda.pz() << " " << lambda.e() << " " << endl;
            TLorentzVector momLambda = TLorentzVector(lambda.px(), lambda.py(), lambda.pz(), lambda.e());
            for(int jPart=3; jPart<pythia.event.size(); jPart++) {
              cout << pythia.event[jPart].id() << endl;
              if(abs(pythia.event[jPart].id())==211) {
                Particle pion = pythia.event[jPart];
                cout << "Quadrimomentum of pion: " << pion.px() << " " << pion.py() << " " << pion.pz() << " " << pion.e() << " " << endl;
                TLorentzVector momPion = TLorentzVector(pion.px(), pion.py(), pion.pz(), pion.e());
                TString histoKey =  std::to_string(pythia.event[iPart].id()) + std::to_string(pythia.event[jPart].id());
                float kStar = RelativePairMomentum(momPion, momLambda);
                cout << "Histo key: " << histoKey << endl;
                hSEPairs[histoKey]->Fill(kStar); 
              }
            }
          }
        }
      }

    std::string outFileName = "./SimLambda1520_" + std::to_string(seed) + ".root";  
    TFile oFile(outFileName.data(), "recreate");
    oFile.cd();
    for(const auto &hSEPair : hSEPairs) {
        hSEPairs[hSEPair.first]->Write();
    }
    fPt->Write();
    fDecay->Write();
    oFile.Close();
}