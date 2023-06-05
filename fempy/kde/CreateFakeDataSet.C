#include "Riostream.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"

void CreateFakeDataSet(TString inFileName = "", TString oFileName = "", Double_t scaleFactor = 1, bool isMC = false) {
    if (isMC) {
        inFileName = "/home/daniel/alice/CharmingAnalyses/DKDpi/oton_mctruth/mcgp/AnalysisResults_all.root";
        oFileName = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/mcgp/AnalysisResults.root";
    } else {
        inFileName = "/home/daniel/alice/CharmingAnalyses/DKDpi/oton_mctruth/data/AnalysisResults_all.root";
        oFileName = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults.root";
    }

    TFile *inFile = new TFile(inFileName);

    std::vector<TString> pairCombs = {"Particle0_Particle2", "Particle1_Particle3"};
    std::map<TString, TString> regions = {{"sgn", ""}, {"sbl", "SBLeft_"}, {"sbr", "SBRight_"}};

    std::map<TString, TH2D *> mapSE;
    std::map<TString, TH2D *> mapME;

    for (auto reg : regions) {
        TString dir = Form("HM_CharmFemto_%sDpion_Results0", reg.second.Data());

        // common for data and MC
        for (auto pair : pairCombs) {
            auto histSE = (TH2D *)inFile->Get(Form("%s/%s", dir.Data(), dir.Data()))
                              ->FindObject(pair)
                              ->FindObject(Form("SEMultDist_%s", pair.Data()));
            mapSE.insert({Form("%s_%s", reg.first.Data(), pair.Data()), histSE});

            auto histME = (TH2D *)inFile->Get(Form("%s/%s", dir.Data(), dir.Data()))
                              ->FindObject(pair)
                              ->FindObject(Form("MEMultDist_%s", pair.Data()));
            mapME.insert({Form("%s_%s", reg.first.Data(), pair.Data()), histME});
        }
        // save MC truth at reco and gen levels
    }
    if (isMC) {
        // Reconstruction level
        TString dirRec = Form("HM_CharmFemto_DpionMCtruthReco_Results0");
        for (auto pair : pairCombs) {
            auto histSE = (TH2D *)inFile->Get(Form("%s/%s", dirRec.Data(), dirRec.Data()))
                              ->FindObject(pair)
                              ->FindObject(Form("SEMultDist_%s", pair.Data()));
            mapSE.insert({Form("MCRec_sgn_%s", pair.Data()), histSE});

            auto histME = (TH2D *)inFile->Get(Form("%s/%s", dirRec.Data(), dirRec.Data()))
                              ->FindObject(pair)
                              ->FindObject(Form("MEMultDist_%s", pair.Data()));
            mapME.insert({Form("MCRec_sgn_%s", pair.Data()), histME});
        }

        // Generation level
        TString dirGen = Form("HM_CharmFemto_DpionMCtruthGen_Results0");
        for (auto pair : pairCombs) {
            auto histSE = (TH2D *)inFile->Get(Form("%s/%s", dirGen.Data(), dirGen.Data()))
                              ->FindObject(pair)
                              ->FindObject(Form("SEMultDist_%s", pair.Data()));
            mapSE.insert({Form("MCGen_sgn_%s", pair.Data()), histSE});

            auto histME = (TH2D *)inFile->Get(Form("%s/%s", dirGen.Data(), dirGen.Data()))
                              ->FindObject(pair)
                              ->FindObject(Form("MEMultDist_%s", pair.Data()));
            mapME.insert({Form("MCGen_sgn_%s", pair.Data()), histME});
        }
    }
    inFile->Close();

    TFile *oFile = new TFile(oFileName, "recreate");

    Double_t kStar = 0;
    Double_t mult = 0;

    std::cout << "Creating Same-Event distributions ..." << std::endl;
    for (const auto &[label, histSE] : mapSE) {

        int nSamplSE = int(scaleFactor * gRandom->Poisson(histSE->GetEntries()));
        TString treeName = Form("treeDpi_SE_%s", label.Data());
        TTree *treeSE = new TTree(Form("treeDpi_SE_%s", label.Data()), Form("treeDpi SE %s", label.Data()));
        treeSE->SetAutoSave(0); // avoid saving backup cycles
        treeSE->Branch("kStar", &kStar, "kStar/D");
        treeSE->Branch("mult", &mult, "mult/D");

        for (int i = 0; i < nSamplSE; i++) {
            histSE->GetRandom2(kStar, mult); // transorm to MeV/c;
            kStar *= 1000;
            treeSE->Fill();
        }
        treeSE->Write();
        std::cout << label << "   " << histSE << "   entries: " << treeSE->GetEntries() << std::endl;
    }

    std::cout << "\nCreating Mixed-Event distributions ..." << std::endl;
    // todo: too many events in the mc truth gen level? cannot sample from 2D distrib
    for (const auto &[label, histME] : mapME) {

        int nSamplME = int(scaleFactor * gRandom->Poisson(histME->GetEntries()));
        TTree *treeME = new TTree(Form("treeDpi_ME_%s", label.Data()), Form("treeDpi ME %s", label.Data()));
        treeME->SetAutoSave(0); // avoid saving backup cycles
        treeME->Branch("kStar", &kStar, "kStar/D");
        treeME->Branch("mult", &mult, "mult/D");

        for (int i = 0; i < nSamplME; i++) {
            histME->GetRandom2(kStar, mult); // transorm to MeV/c;
            kStar *= 1000;
            treeME->Fill();
        }
        treeME->Write();
        std::cout << label << "   " << histME << "   entries: " << treeME->GetEntries()
                  << "   histME entries: " << histME->GetEntries() << std::endl;
        // histME->Write(Form("hhDpi_Mult_kStar_ME_%s", label));
    }

    for (const auto &[label, histSE] : mapSE) {
        histSE->Clone(Form("hhDpi_Mult_kStar_SE_%s", label.Data()))->Write();
    }
    for (const auto &[label, histME] : mapME) {
        histME->Clone(Form("hhDpi_Mult_kStar_ME_%s", label.Data()))->Write();
    }
    oFile->Close();
}
