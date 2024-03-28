#include "TH1F.h"
#include "TClonesArray.h"
#include "AliPythia8.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"
#include "yaml-cpp/yaml.h"
//#include "Pythia8/Pythia8.h"

#include <array>
#include <map>
#include <string>
#include <vector>


#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>

enum tunes {kMonash, kCRMode0, kCRMode2, kCRMode3};

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

void CheckDecayChannels(

/*  
    Script that computes the SE distribution for a specific mother and
    specific decay chain, to be specified in a .yml file
*/
    int seed = 42,
    std::string cfgFilePath = "/home/mdicostanzo/an/LPi/Simulation/outputs/Sigma1385_recovered.yml",
    std::string MCFilePath="/home/mdicostanzo/an/LPi/Trains/02_allpc/mc/data/AnalysisResultsAllPC.root",
    int nEvents = 100000,
    tunes tune = kMonash) {

    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

    // after this line one can select the configuration for 
    // his simulation before calling the Initialize method
    AliPythia8 pythia8; 

    // select physics processes
    if (tune == kMonash) {
        pythia8.ReadString(Form("Tune:pp = 14"));
        pythia8.ReadString("SoftQCD:all = on");
    } else if (tune == kCRMode0) {
        pythia8.ReadString(Form("Tune:pp = 14"));
        pythia8.ReadString("ColourReconnection:mode = 1");
        pythia8.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia8.ReadString("ColourReconnection:m0 = 2.9");
        pythia8.ReadString("ColourReconnection:allowJunctions = on");
        pythia8.ReadString("ColourReconnection:junctionCorrection = 1.43");
        pythia8.ReadString("ColourReconnection:timeDilationMode = 0");
        pythia8.ReadString("StringPT:sigma = 0.335");
        pythia8.ReadString("StringZ:aLund = 0.36");
        pythia8.ReadString("StringZ:bLund = 0.56");
        pythia8.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia8.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia8.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia8.ReadString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia8.ReadString("BeamRemnants:remnantMode = 1");
        pythia8.ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode2) {
        pythia8.ReadString(Form("Tune:pp = 14"));
        pythia8.ReadString("ColourReconnection:mode = 1");
        pythia8.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia8.ReadString("ColourReconnection:m0 = 0.3");
        pythia8.ReadString("ColourReconnection:allowJunctions = on");
        pythia8.ReadString("ColourReconnection:junctionCorrection = 1.20");
        pythia8.ReadString("ColourReconnection:timeDilationMode = 2");
        pythia8.ReadString("ColourReconnection:timeDilationPar = 0.18");
        pythia8.ReadString("StringPT:sigma = 0.335");
        pythia8.ReadString("StringZ:aLund = 0.36");
        pythia8.ReadString("StringZ:bLund = 0.56");
        pythia8.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia8.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia8.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia8.ReadString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia8.ReadString("BeamRemnants:remnantMode = 1");
        pythia8.ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode3) {
        pythia8.ReadString(Form("Tune:pp = 14"));
        pythia8.ReadString("ColourReconnection:mode = 1");
        pythia8.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia8.ReadString("ColourReconnection:m0 = 0.3");
        pythia8.ReadString("ColourReconnection:allowJunctions = on");
        pythia8.ReadString("ColourReconnection:junctionCorrection = 1.15");
        pythia8.ReadString("ColourReconnection:timeDilationMode = 3");
        pythia8.ReadString("ColourReconnection:timeDilationPar = 0.073");
        pythia8.ReadString("StringPT:sigma = 0.335");
        pythia8.ReadString("StringZ:aLund = 0.36");
        pythia8.ReadString("StringZ:bLund = 0.56");
        pythia8.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia8.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia8.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia8.ReadString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia8.ReadString("BeamRemnants:remnantMode = 1");
        pythia8.ReadString("BeamRemnants:saturation = 5");
    }


    // set seed for simulation
    pythia8.ReadString(Form("Random:seed %d", seed));
    pythia8.ReadString("Random:setSeed = on");
    //pythia8.ReadString("Init:showProcesses = off");
    //pythia8.ReadString("Init:showMultipartonInteractions = off");

    // Open yaml configuration file       
    YAML::Node cfgFile = YAML::LoadFile(cfgFilePath.data());

    // switch on specific decay processes
    const YAML::Node& processes = cfgFile["switchon"];
    std::vector<std::string> initialStates;
    std::vector<std::vector<string>> finalStates;
    for(int proc = 0; proc < processes.size(); proc++){
        std::string initpart = processes[proc]["ch"][0].as<std::string>();
        std::string switchoffpart = initpart + ":onMode = off";
        //cout << switchoffpart << endl;
        pythia8.ReadString(switchoffpart.data());
        std::string finalpart = ""; 
        for(int finpart = 1; finpart < processes[proc]["ch"].size(); finpart++){
            finalpart += processes[proc]["ch"][finpart].as<std::string>();
            finalpart += " ";
        }
        finalpart.pop_back();
        std::string switchonchn = initpart + ":onIfMatch = " + finalpart;
        pythia8.ReadString(switchonchn.data());
    }

    // node including the specific chains to be checked
    const YAML::Node& channels = cfgFile["channels"];
    std::map<std::string, TH1F*> hSEMothers;           
    std::map<std::string, TH1F*> hSEPairs;
    std::map<std::string, TH1F*> hSEPairsSmeared;           
    std::vector<int> mothers;
    std::vector<std::string> hMotherKeys;
    std::vector<std::string> hPairKeys;

    // specific decay histories for the first and second particles of the pair,
    // the decay chains are given as vector so an arbitrary number of layers can
    // be handled
    std::vector<std::vector<int>> firstPartChains;
    std::vector<std::vector<int>> secondPartChains;       
    for(int ch = 0; ch < channels.size(); ch++){
        std::string hHistoMotherName = "hSE_" + channels[ch]["mothername"].as<std::string>();   
        std::string hHistoSmearedMotherName = hHistoMotherName + "_smeared";   
        hMotherKeys.push_back(channels[ch]["mothername"].as<std::string>());
        hSEMothers.insert({hMotherKeys.back(), new TH1F(hHistoMotherName.data(), ";#it{k*};Counts", 1500, 0., 6.)});
        
        mothers.push_back(channels[ch]["motherPDG"].as<int>());
        firstPartChains.push_back(channels[ch]["part1chain"].as<std::vector<int>>());
        secondPartChains.push_back(channels[ch]["part2chain"].as<std::vector<int>>());
        
        std::string hHistoPairName = "hSE_" + channels[ch]["part1"].as<std::string>() + channels[ch]["part2"].as<std::string>();   
        hPairKeys.push_back(std::to_string(secondPartChains[ch].back()) + std::to_string(firstPartChains[ch].back()));
        hSEPairs.insert({hPairKeys.back(), new TH1F(hHistoPairName.data(), ";#it{k*};Counts", 1500, 0., 6.)});
        std::string hHistoSmearedPairName = "hSE_" + channels[ch]["part1"].as<std::string>() + 
                                            channels[ch]["part2"].as<std::string>() + "_smeared";   
        hSEPairsSmeared.insert({hPairKeys.back(), new TH1F(hHistoSmearedPairName.data(), ";#it{k*};Counts", 1500, 0., 6.)});

    }

    pythia8.Initialize(2212, 2212, 13000);

    const std::string oFileName = cfgFile["oFileDir"].as<std::string>() + cfgFile["oFileName"].as<std::string>() 
                                  + "_" + std::to_string(seed) + ".root";
    TFile oFile(oFileName.data(), "recreate");       

    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia8.GenerateEvent();
        pythia8.ImportParticles(particles, "All");
        if(iEvent % 100 == 0){
            cout << "Processing " << iEvent << " ..." << endl;
        }
        
        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
            for(unsigned ch = 0; ch < mothers.size(); ch++){

                // check if particle is one of the mothers we are interested in
                if(mothers[ch] == particle->GetPdgCode()){
                    
                    // First particle
                    int firstDecChain = iPart;
                    int chainLayer = 0;
                    while(chainLayer < firstPartChains[ch].size()-1){
                        TParticle* chainLayerPart = dynamic_cast<TParticle*>(particles->At(firstDecChain));
                        
                        // pick the daughter having the PDG code specified in the decay chain in the .yml file
                        for(int i = 0; i < chainLayerPart->GetNDaughters(); i++){
                            TParticle* chainLayerDaug = (TParticle*) particles->At(chainLayerPart->GetDaughter(i));
                            // in the next iteration, the daughters of this particle will be checked
                            if(chainLayerDaug->GetPdgCode() == firstPartChains[ch][chainLayer+1]){
                                firstDecChain = chainLayerPart->GetDaughter(i);
                            }
                        }
                        chainLayer++;
                    }

                    // Same logic for the second particle
                    int secondDecChain = iPart;
                    chainLayer = 0;
                    while(chainLayer < secondPartChains[ch].size()-1){
                        TParticle* chainLayerPart = dynamic_cast<TParticle*>(particles->At(secondDecChain));
                        for(int i = 0; i < chainLayerPart->GetNDaughters(); i++){
                            TParticle* chainLayerDaug = (TParticle*) particles->At(chainLayerPart->GetDaughter(i));
                            if(chainLayerDaug->GetPdgCode() == secondPartChains[ch][chainLayer+1]){
                                secondDecChain = chainLayerPart->GetDaughter(i);
                            }
                        }
                        chainLayer++;
                    }
                    
                    TParticle* pair1 = dynamic_cast<TParticle*>(particles->At(firstDecChain));
                    TParticle* pair2 = dynamic_cast<TParticle*>(particles->At(secondDecChain));
                    int PDG1 = pair1->GetPdgCode();
                    int PDG2 = pair2->GetPdgCode();
                    if(PDG1 == firstPartChains[ch].back() && PDG2 == secondPartChains[ch].back()){
                        TLorentzVector mom1 = TLorentzVector(pair1->Px(), pair1->Py(), pair1->Pz(), pair1->Energy());
                        TLorentzVector mom2 = TLorentzVector(pair2->Px(), pair2->Py(), pair2->Pz(), pair2->Energy());
                        hSEMothers[hMotherKeys[ch]]->Fill(RelativePairMomentum(mom1, mom2));
                        hSEPairs[std::to_string(PDG2) + std::to_string(PDG1)]->Fill(RelativePairMomentum(mom1, mom2));
                    }
                }
            }
        }
    }           
    
    // Load smearing matrices
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
    oFile.cd();

    // Smearing appliable only to histos relative to pairs, not to specific mother particles
    for(int iPairKey=0; iPairKey<hPairKeys.size(); iPairKey++) {
      for(int iBin=0; iBin<hSEPairsSmeared[hPairKeys[iPairKey]]->GetNbinsX(); iBin++) {
        TH1D *hProjY = smearMatrices[iPairKey]->ProjectionY(Form("iBin_%i", iBin+1), iBin+1, iBin+1);
        if(hProjY->Integral() > 0) {    
            hProjY->Scale(hSEPairs[hPairKeys[iPairKey]]->GetBinContent(iBin+1)/hProjY->Integral());
            hSEPairsSmeared[hPairKeys[iPairKey]]->Add(hProjY);
        }
      }
      for(int iBin=0; iBin<hSEPairsSmeared[hPairKeys[iPairKey]]->GetNbinsX(); iBin++) {
        hSEPairsSmeared[hPairKeys[iPairKey]]->SetBinError(iBin+1, TMath::Sqrt(hSEPairsSmeared[hPairKeys[iPairKey]]->GetBinContent(iBin+1)));
      }
    }
    
    oFile.cd();
    for(const auto &hSEPair : hSEPairs) {
        hSEPairs[hSEPair.first]->Write();
        hSEPairsSmeared[hSEPair.first]->Write();
    }
    for(const auto &hSEMother : hSEMothers) {
        hSEMothers[hSEMother.first]->Write();
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
    resumLPiPlSmear->Add(hSEPairsSmeared["3122211"]);
    resumLPiPlSmear->Add(hSEPairsSmeared["-3122-211"]);
    resumLPiPlSmear->Write();
    TH1F * resumLPiMinSmear = new TH1F("3122-211_-3122211_smeared", ";#it{k*};Counts", 1500, 0., 6.);
    resumLPiMinSmear->Add(hSEPairsSmeared["3122-211"]);
    resumLPiMinSmear->Add(hSEPairsSmeared["3122-211"]);
    resumLPiMinSmear->Write();

}