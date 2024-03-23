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

#include "Pythia8/Pythia.h"

enum tunes {kMonash = 0, kCRMode0, kCRMode2, kCRMode3};

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

// Macro for files NewSimStandalonePythia
void SimLPiBkg(int seed=42, int nEvents=250000, tunes = kCRMode2,
               std::string MCFilePath="/home/mdicostanzo/an/LPi/Trains/02_allpc/mc/data/AnalysisResultsAllPC.root",
               std::string outFilePath= "/home/mdicostanzo/an/LPi/Simulation/outputs/SimLPi_allbkg_new") {
    
    Pythia pythia;
    
    if (tune == kMonash) {
        pythia.ReadString(Form("Tune:pp = 14"));
    } else if (tune == kCRMode0) {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 2.9");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.43");
        pythia.ReadString("ColourReconnection:timeDilationMode = 0");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode2) {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 0.3");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.20");
        pythia.ReadString("ColourReconnection:timeDilationMode = 2");
        pythia.ReadString("ColourReconnection:timeDilationPar = 0.18");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode3) {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 0.3");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.15");
        pythia.ReadString("ColourReconnection:timeDilationMode = 3");
        pythia.ReadString("ColourReconnection:timeDilationPar = 0.073");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    }

    // set seed for simulation
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    
    //pythia.readString("SoftQCD:all = on");
    //pythia.readString("Beams:eCM = 13000");

    // init
    pythia.init();

    // save output
    std::string outFileName = outFilePath + "_" + std::to_string(seed) + ".root";   
    TFile oFile(outFileName.data(), "recreate");

    // define histograms
    std::map<TString, TH1D*> hSEMothersPairs;
    std::map<TString, TH1D*> hSEMothersPairsSmeared;
    std::vector<TString> histoKeys;
    
    // define vectors to store particles for every event
    std::vector<std::tuple<int, int, int, int>> pions;
    std::vector<std::tuple<int, int, int, int>> lambdas;
    std::vector<std::tuple<int, int, int, int>> pionsDirectSigma;
    std::vector<std::tuple<int, int, int, int>> lambdasDirectSigma;

    int countLambda1520 = 0;
    for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
        if(iEvent%10000 == 0) {
            std::cout << Form("Processing event... ", iEvent) << endl;
        }
        if (!pythia.next()) {
            continue;
        } 
        for (int iPart = 2; iPart < pythia.event.size(); iPart++) {
            Particle part = pythia.event.at(iPart);
            if(abs(part.id()) == 3124) {
                countLambda1520++;
            }
            if(abs(part.id()) == 211 || abs(part.id()) == 3122) {
                
                if(abs(part.eta())>0.8) {
                    continue;
                } 
                if(abs(part.pT())<0.3) {
                    continue;
                } 

                if(part.daughter1() == part.daughter2() && part.daughter1()>0 && part.daughter2()>0) {
                    continue;
                }  

                if(pythia.event.at(part.mother1()).isHadron()) {

                    int moth1Idx = part.mother1();
                    int moth2Idx = part.mother2();
                    Particle intermPart = pythia.event.at(moth1Idx);
                
                    int partPrimMothIdx = moth1Idx;
                    int checkDirectDecay = 0;
                    bool foundMother = false;
                    if(abs(intermPart.status())>=81 && abs(intermPart.status())<=89) {
                        foundMother = true;
                     }
                    while(!foundMother) {
                        partPrimMothIdx = moth1Idx;
                        Particle intermPart = pythia.event.at(moth1Idx);
                        moth1Idx = intermPart.mother1();
                        moth2Idx = intermPart.mother2();
                        Particle moth1 = pythia.event.at(moth1Idx);
                        Particle moth2 = pythia.event.at(moth2Idx);
                        
                        if(abs(intermPart.status())>=81 && abs(intermPart.status())<=89 || !moth1.isHadron()) {
                            foundMother = true;
                        }
                        
                        checkDirectDecay++;
                        if(checkDirectDecay>=3)  {
                            if(intermPart.id() == 90) return;
                        }
                    }
                    
                    Particle primMoth = pythia.event.at(partPrimMothIdx);
                    int finalPartPdg = part.id(); 
                    int primMothPdg = primMoth.id(); 
                    std::tuple<int, int, int, int> partInfo = {finalPartPdg, iPart, primMothPdg, partPrimMothIdx};
                    if(finalPartPdg == -211 && primMothPdg == 3122 && checkDirectDecay == 1) continue;  // pair cleaner
                    if(finalPartPdg == 211 && primMothPdg == -3122 && checkDirectDecay == 1) continue;  // pair cleaner
                    if( (abs(primMothPdg) == 3114 || abs(primMothPdg) == 3224 ) && checkDirectDecay == 0) {
                        if(abs(finalPartPdg) == 211) {
                            pionsDirectSigma.push_back(partInfo);
                        } else if(abs(finalPartPdg) == 3122) {
                            lambdasDirectSigma.push_back(partInfo);
                        }    
                    } else {
                        if(abs(finalPartPdg) == 211) {
                            pions.push_back(partInfo);
                        } else if(abs(finalPartPdg) == 3122) {
                            lambdas.push_back(partInfo);
                        }
                    }
                }
            }
        }
                
        for(int iPion = 0; iPion < pions.size(); iPion++) {
            Particle pion = pythia.event.at(std::get<1>(pions[iPion]));
            TLorentzVector momPion = TLorentzVector(pion.px(), pion.py(), pion.pz(), pion.e());
            for(int iLambda = 0; iLambda < lambdas.size(); iLambda++) {
                if( std::get<2>(pions[iPion]) == std::get<2>(lambdas[iLambda]) && std::get<3>(pions[iPion]) == std::get<3>(lambdas[iLambda])) {
                    
                    TString pairKey = std::to_string(std::get<0>(pions[iPion])) + std::to_string(std::get<0>(lambdas[iLambda]));
                    TString motherKey = std::to_string(std::get<2>(pions[iPion]));
                    TString histoKey = "hSE_" + pairKey + "_" + motherKey; 
                    
                    if(std::find(histoKeys.begin(), histoKeys.end(), histoKey.Data()) == histoKeys.end()) { 
                        histoKeys.push_back(histoKey);
                        TString histoKeySmeared = "hSE_" + pairKey + "_" + motherKey + "_SMEARED"; 
                        hSEMothersPairs.insert({histoKey, new TH1D(histoKey.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
                        hSEMothersPairsSmeared.insert({histoKey, new TH1D(histoKeySmeared.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
                    }
                    Particle lambda = pythia.event.at(std::get<1>(lambdas[iLambda]));
                    TLorentzVector momLambda = TLorentzVector(lambda.px(), lambda.py(), lambda.pz(), lambda.e());
                    float kStar = RelativePairMomentum(momPion, momLambda);
                    hSEMothersPairs[histoKey]->Fill(kStar);
                }
            }
        }
        for(int iDirPion = 0; iDirPion < pionsDirectSigma.size(); iDirPion++) {
            Particle pion = pythia.event.at(std::get<1>(pionsDirectSigma[iDirPion]));
            TLorentzVector momPion = TLorentzVector(pion.px(), pion.py(), pion.pz(), pion.e());
            for(int iDirLambda = 0; iDirLambda < lambdasDirectSigma.size(); iDirLambda++) {
                if(std::get<2>(pionsDirectSigma[iDirPion]) == std::get<2>(lambdasDirectSigma[iDirLambda]) && 
                   std::get<3>(pionsDirectSigma[iDirPion]) == std::get<3>(lambdasDirectSigma[iDirLambda]) ) {
                    
                    TString pairKey = std::to_string(std::get<0>(pionsDirectSigma[iDirPion])) + std::to_string(std::get<0>(lambdasDirectSigma[iDirLambda]));
                    TString motherKey = std::to_string(std::get<2>(pionsDirectSigma[iDirPion]));
                    TString histoKey = "hSE_" + pairKey + "_" + motherKey + "_" + "DIRECT_SIGMA"; 
                    if(std::find(histoKeys.begin(), histoKeys.end(), histoKey.Data()) == histoKeys.end()) { 
                        histoKeys.push_back(histoKey);
                        TString histoKeySmeared = "hSE_" + pairKey + "_" + motherKey + "_" + "DIRECT_SIGMA" + "_SMEARED"; 
                        hSEMothersPairs.insert({histoKey, new TH1D(histoKey.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
                        hSEMothersPairsSmeared.insert({histoKey, new TH1D(histoKeySmeared.Data(), ";#it{k*};Counts", 1500, 0., 6.)});
                    }
                    Particle lambda = pythia.event.at(std::get<1>(lambdasDirectSigma[iDirLambda]));
                    TLorentzVector momLambda = TLorentzVector(lambda.px(), lambda.py(), lambda.pz(), lambda.e());
                    float kStar = RelativePairMomentum(momPion, momLambda);
                    hSEMothersPairs[histoKey]->Fill(kStar);
                }
            }
        }
        
        pions.clear();
        lambdas.clear();
        pionsDirectSigma.clear();
        lambdasDirectSigma.clear();

    }

    std::map<TString, TH2D*> smearMatrices;
    TFile *MCdata = TFile::Open(MCFilePath.data());
    TDirectoryFile *folder = static_cast<TDirectoryFile*>(MCdata->Get("HMResultsQA1001"));
    TList *toplist = static_cast<TList*>(folder->Get("HMResultsQA1001"));
    TList *QAList = dynamic_cast<TList*>(toplist->FindObject("PairQA"));

    TList *pairList02 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle0_Particle2"));
    TH2D *smearMomentumMatrix02 = static_cast<TH2D*>(pairList02->FindObject("MomentumResolutionSE_Particle0_Particle2"));
    TString smearMatrix02 = "SMEARMATRIX_2113122_"; 
    smearMatrices.insert({smearMatrix02, smearMomentumMatrix02});

    TList *pairList03 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle0_Particle3"));
    TH2D *smearMomentumMatrix03 = static_cast<TH2D*>(pairList03->FindObject("MomentumResolutionSE_Particle0_Particle3"));
    TString smearMatrix03 = "SMEARMATRIX_211-3122_"; 
    smearMatrices.insert({smearMatrix03, smearMomentumMatrix03});

    TList *pairList12 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle1_Particle2"));
    TH2D *smearMomentumMatrix12 = static_cast<TH2D*>(pairList12->FindObject("MomentumResolutionSE_Particle1_Particle2"));
    TString smearMatrix12 = "SMEARMATRIX_-2113122_"; 
    smearMatrices.insert({smearMatrix12, smearMomentumMatrix12});

    TList *pairList13 = dynamic_cast<TList*>(QAList->FindObject("QA_Particle1_Particle3"));
    TH2D *smearMomentumMatrix13 = static_cast<TH2D*>(pairList13->FindObject("MomentumResolutionSE_Particle1_Particle3"));
    TString smearMatrix13 = "SMEARMATRIX_-211-3122_"; 
    smearMatrices.insert({smearMatrix13, smearMomentumMatrix13});

    MCdata->Close();

    oFile.cd();
    std::vector<std::string> directories = {"_2113122", "_211-3122", "_-2113122", "_-211-3122"};
    std::vector<std::string> smearedDirectories = {"_2113122_SMEARED", "_211-3122_SMEARED", "_-2113122_SMEARED", "_-211-3122_SMEARED"};
    std::vector<std::string> mergedHistoNames = {"_2113122_MERGED", "_211-3122_MERGED", "_-2113122_MERGED", "_-211-3122_MERGED"};
    std::vector<TH1D*> mergedHistos;
    std::vector<std::string> mergedHistoNoDirectNames = {"_2113122_MERGED_NODIRECT", "_211-3122_MERGED_NODIRECT", "_-2113122_MERGED_NODIRECT", "_-211-3122_MERGED_NODIRECT"};
    std::vector<TH1D*> mergedHistosNoDirect;
    std::vector<std::string> mergedHistoSmearedNames = {"_2113122_MERGED_SMEARED", "_211-3122_MERGED_SMEARED", "_-2113122_MERGED_SMEARED", "_-211-3122_MERGED_SMEARED"};
    std::vector<TH1D*> mergedHistosSmeared;
    std::vector<std::string> mergedHistoSmearedNoDirectNames = {"_2113122_MERGED_NODIRECT_SMEARED", "_211-3122_MERGED_NODIRECT_SMEARED", "_-2113122_MERGED_NODIRECT_SMEARED", "-211-3122_MERGED_NODIRECT_SMEARED"};
    std::vector<TH1D*> mergedHistosSmearedNoDirect;    

    // cout << "CREATING DIRECTORIES" << endl;
    for(int iDir=0; iDir<directories.size(); iDir++) {
       oFile.mkdir(directories[iDir].data());
       oFile.mkdir(smearedDirectories[iDir].data());
       mergedHistos.push_back(new TH1D(mergedHistoNames[iDir].data(), ";#it{k*};Counts", 1500, 0., 6.));
       mergedHistosSmeared.push_back(new TH1D(mergedHistoSmearedNames[iDir].data(), ";#it{k*};Counts", 1500, 0., 6.));
       mergedHistosNoDirect.push_back(new TH1D(mergedHistoNoDirectNames[iDir].data(), ";#it{k*};Counts", 1500, 0., 6.));
       mergedHistosSmearedNoDirect.push_back(new TH1D(mergedHistoSmearedNoDirectNames[iDir].data(), ";#it{k*};Counts", 1500, 0., 6.));
    }

    // cout << "CREATED DIRECTORIES MERGED" << endl;

    for(const auto &hSEMothersPair : hSEMothersPairs) {
       for(int iDir=0; iDir<directories.size(); iDir++) {
           std::string checkPair = directories[iDir] + "_";
         
            if(hSEMothersPair.first.Contains(checkPair)) {
                
                oFile.cd(directories[iDir].data());
                hSEMothersPairs[hSEMothersPair.first]->Write();
                
                mergedHistos[iDir]->Add(hSEMothersPairs[hSEMothersPair.first]);
                
                // cout << "CIAO1" << endl;
                TString smearMatrixKey = "SMEARMATRIX" + checkPair;
                // cout << smearMatrixKey << endl;
                // cout << "CIAO2" << endl;
                for(int iBin=0; iBin<hSEMothersPairs[hSEMothersPair.first]->GetNbinsX(); iBin++) {
                // cout << "CIAO3" << endl;
                    TH1D *hProjY = smearMatrices[smearMatrixKey]->ProjectionY(Form("iBin_%i", iBin+1), iBin+1, iBin+1);
                // cout << "CIAO4" << endl;
                    if(hProjY->Integral() > 0) {    
                        hProjY->Scale(hSEMothersPairs[hSEMothersPair.first]->GetBinContent(iBin+1)/hProjY->Integral());
                        hSEMothersPairsSmeared[hSEMothersPair.first]->Add(hProjY);
                    }
                }
                for(int iBin=0; iBin<hSEMothersPairs[hSEMothersPair.first]->GetNbinsX(); iBin++) {
                   hSEMothersPairsSmeared[hSEMothersPair.first]->SetBinError(iBin+1, TMath::Sqrt(hSEMothersPairsSmeared[hSEMothersPair.first]->GetBinContent(iBin+1)));
                }

                mergedHistosSmeared[iDir]->Add(hSEMothersPairsSmeared[hSEMothersPair.first]);
                
                oFile.cd(smearedDirectories[iDir].data());
                hSEMothersPairsSmeared[hSEMothersPair.first]->Write();
                
                if(!hSEMothersPair.first.Contains("DIRECT_SIGMA")) {
                    mergedHistosNoDirect[iDir]->Add(hSEMothersPairs[hSEMothersPair.first]);
                    for(int iBin=0; iBin<hSEMothersPairs[hSEMothersPair.first]->GetNbinsX(); iBin++) {
                        TH1D *hProjY = smearMatrices[smearMatrixKey]->ProjectionY("smear_bin", iBin+1, iBin+1);
                        if(hProjY->Integral() > 0) {
                            hProjY->Scale(hSEMothersPairs[hSEMothersPair.first]->GetBinContent(iBin+1)/hProjY->Integral());
                            mergedHistosSmearedNoDirect[iDir]->Add(hProjY);
                        }
                    }
                    for(int iBin=0; iBin<hSEMothersPairs[hSEMothersPair.first]->GetNbinsX(); iBin++) {
                        mergedHistosSmearedNoDirect[iDir]->SetBinError(iBin+1, TMath::Sqrt(mergedHistosSmearedNoDirect[iDir]->GetBinContent(iBin+1)));
                    }                
                }
            }
        }
    }
    
    for(int iDir=0; iDir<directories.size(); iDir++) {
        oFile.cd(directories[iDir].data());
        mergedHistos[iDir]->Write();
        mergedHistosNoDirect[iDir]->Write();
    }

    // cout << "WRITING MERGED" << endl;
    
    for(int iDir=0; iDir<smearedDirectories.size(); iDir++) {
        oFile.cd(smearedDirectories[iDir].data());
        mergedHistosSmeared[iDir]->Write();
        mergedHistosSmearedNoDirect[iDir]->Write();
    }

    oFile.Close();
}
