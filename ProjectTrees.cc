#include <map>
#include <string>

#include "ROOT/RDataFrame.hxx"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TTree.h"
#include "fempy/utils/Analysis.hxx"


void ProjectTrees(std::string inFileName = "/data/DstarPi/tree_pc/mcgp/AnalysisResults_3998.root",
                  std::string oFileName = "/home/daniel/an/DstarPi/proj/Projection_mcgp_test.root",
                  std::string pair = "DstarPi", std::string version = "");

bool MassSelection(const double &mass, const double &pt, const int &hpdg, const std::string &massRegion,
                   const double fNSigmaMass = 2., double fNSigmaOffsetSideband = 5., double fSidebandWidth = 0.2,
                   double fLowerDstarRemoval = 1.992, double fUpperDstarRemoval = 2.028) {
    if (massRegion == "any") {
        return true;
    }

    // simple parametrisation from D+ in 5.02 TeV
    double massMean =
        TDatabasePDG::Instance()->GetParticle(hpdg)->Mass() + 0.0025;  // mass shift observed in all Run2 data samples
                                                                       // for all D-meson species
    double massWidth = 0.;
    switch (hpdg) {
        case 411:
            massWidth = 0.006758 + pt * 0.0005124;
            break;
        case 413:
            Double_t mDstarPDG = TDatabasePDG::Instance()->GetParticle(413)->Mass();
            Double_t mD0PDG = TDatabasePDG::Instance()->GetParticle(421)->Mass();
            massMean = mDstarPDG - mD0PDG;  // no extra mass shift because it is deltamass
            massWidth = 0.00124673 - pt * 0.000340426 + pt * pt * 4.40729e-05;
            if (pt > 4 && pt < 5)
                massWidth = 0.00104329 - 0.000113275 * pt;
            else if (pt >= 5)
                massWidth = 0.000519861 - 8.58874e-06 * pt;
            break;
    }

    // select D mesons mass window
    double fLowerMassSelection, fUpperMassSelection;

    if (massRegion == "sgn") {
        fLowerMassSelection = massMean - fNSigmaMass * massWidth;
        fUpperMassSelection = massMean + fNSigmaMass * massWidth;
    } else if (massRegion == "sbl") {
        fLowerMassSelection = massMean - fNSigmaOffsetSideband * massWidth - fSidebandWidth;
        fUpperMassSelection = massMean - fNSigmaOffsetSideband * massWidth;
    } else if (massRegion == "sbr") {
        fLowerMassSelection = massMean + fNSigmaOffsetSideband * massWidth;
        fUpperMassSelection = massMean + fNSigmaOffsetSideband * massWidth + fSidebandWidth;

        if (hpdg == 411) {
            // additional removal of D*
            if (mass > fLowerDstarRemoval && mass < fUpperDstarRemoval) {
                return false;
            }
        }
    } else {
        printf("region not supported. Exit!\n");
        exit(1);
    }

    if (mass > fLowerMassSelection && mass < fUpperMassSelection) {
        return true;
    }

    return false;
}

std::string GetRegionFromKey(std::string key) {
    if (key.find("_SBRight_DstarPion_") != std::string::npos)
        return std::string("sbr");
    else if (key.find("_DstarPion_") != std::string::npos)
        return std::string("sgn");
    else {
        printf("Error: region not implemented. Exit!");
        exit(1);
    }
}

bool replace(std::string &str, const std::string &from, const std::string &to) {
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos) return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

std::string TranslateSelection(const std::string selection, const std::map<std::string, std::string> aliases) {
    std::string translatedSelection = selection;

    bool done = true;
    for (auto alias : aliases) {
        replace(translatedSelection, alias.first, alias.second);
        done = false;
    }
    if (done) {
        return translatedSelection;
    } else
        return TranslateSelection(translatedSelection, aliases);
}

void ProjectTrees(std::string inFileName, std::string oFileName, std::string pair, std::string version) {
    int hpdg;
    std::string pairLabel;

    Analysis analysis(pair, version);

    std::vector<const char *> regions;
    if (pair == "DstarPi") {
        hpdg = 413;
        regions = {"sgn", "sbr"};
        pairLabel = "DstarPion";
    } else if (pair == "DPi") {
        hpdg = 411;
        regions = {"sgn", "sbl", "sbr"};
        pairLabel = "DplusPion";
    } else if (pair == "DK") {
        hpdg = 411;
        regions = {"sgn", "sbl", "sbr"};
        pairLabel = "Dkaon";
    } else {
        printf("not impl\n");
        return;
    }


    // don't overwrite output files: skip if existing
    if(!gSystem->AccessPathName(oFileName.data())) {
        printf("Warning: the file %s already exists. Skipping...\n", oFileName.data());
        return;
    }

    auto inFile = TFile::Open(inFileName.data());
    auto oFile = TFile::Open(oFileName.data(), "recreate");
    oFile->mkdir("distr");
    for (auto &[checkName, _] : analysis.selections) {
        oFile->mkdir(Form("distr/%s", checkName.data()));
    }

    std::vector<std::string> keys = {};
    for (auto key : *inFile->GetListOfKeys()) {
        if (std::string(key->GetName()).find("_Trees") != std::string::npos) keys.push_back(std::string(key->GetName()));
    }

    for (auto key : keys) {
        auto dir = (TDirectory *)inFile->Get(key.data());

        
        for (const char *event : {"SE", "ME"}) {
            for (const char *comb : {"pp", "mm", "pm", "mp"}) {
                
                std::cout << "Projecting: " << key << " " << event << " " << comb << "...";
                auto tree = (TTree *)dir->Get(Form("t%s_%s", event, comb));
                ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager> df = ROOT::RDataFrame(*tree);

                df = df.Define("kStarMeV", "kStar * 1000");
                for (auto &[aliasName, alias] : analysis.aliases) {
                    if (!tree->FindBranch(aliasName.data()))
                        df = df.Define(aliasName.data(), alias.data());
                }

                for (const char* region : regions) {
                    auto lamMassSelection = [hpdg, region](float mass, float pt) { return MassSelection(mass, pt, hpdg, region); };

                    auto dfMass = df.Filter(lamMassSelection, {"heavy_invmass", "heavy_pt"});

                    // std::string region = GetRegionFromKey(key);


                    for (auto &[varName, varSel] : analysis.selections) {

                        oFile->cd(Form("distr/%s", varName.data()));
                        auto dfVar = dfMass.Filter(varSel);

                        const char *histName = Form("h%s_%s_%s", event, comb, region);
                        const char *histTitle = ";#it{k}* (MeV/#it{c});Mult;Counts";
                        auto hDistr = dfVar.Histo2D<float, int>({histName, histTitle, 3000u, 0., 3000., 180u, 0.5, 180.5},
                                                                "kStarMeV", "mult");
                        hDistr->Write();
                    }

                }
                std::cout << " done!" << std::endl;
            }
        }
    }

    oFile->Close();
    printf("Distributions saved to %s\n", oFileName.data());
}
