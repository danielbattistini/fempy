#include <map>
#include <string>

#include "ROOT/RDataFrame.hxx"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

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

void ProjectTrees(std::string inFileName = "/data/DstarPi/tree_pc/mcgp/AnalysisResults_3998.root",
                  std::string oFileName = "/home/daniel/an/DstarPi/proj/Projection_mcgp_test.root",
                  std::string pair = "DstarPi") {
    ROOT::EnableImplicitMT(16);
    int hpdg;
    std::string pairLabel;

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

    const char *centr =
        "(abs(light_eta) < 0.8) & "
        "(light_pt > 0.14) & "
        "(light_pt <4.) & "
        "(((abs(light_nsigtpc) < 3) & (light_pt < 0.5)) | ((abs(light_nsigcomb) < 3) & (light_pt > 0.5))) &"
        "(light_ncrossed > 70) & "
        "(light_ncls > 80) & "
        "(abs(light_dcaz) < 0.3) & "
        "(abs(light_dcaxy) < 0.3)";


    std::map<const char *, std::string> checks = {
        {"oldpckept", Form("(is_oldpcrm == 0) & %s", centr)},
        {"oldpckept_motherpi_eq413", Form("abs(light_motherpdg) == 413 & is_oldpcrm == 0 & %s", centr)},
        {"newpckept_motherpi_eq413", Form("abs(light_motherpdg) == 413 & is_newpcrm == 0 & %s", centr)},
        {"oldpckept_motherpi_neq413", Form("abs(light_motherpdg) != 413 & is_oldpcrm == 0 & %s", centr)},
        {"newpckept_motherpi_neq413", Form("abs(light_motherpdg) != 413 & is_newpcrm == 0 & %s", centr)},
        {"oldpcrm", Form("(is_oldpcrm == 1) & %s", centr)},
        {"oldpcrm_motherpi_eq413", Form("abs(light_motherpdg) == 413 & is_oldpcrm == 1 & %s", centr)},
        {"newpcrm_motherpi_eq413", Form("abs(light_motherpdg) == 413 & is_newpcrm == 1 & %s", centr)},
        {"oldpcrm_motherpi_neq413", Form("abs(light_motherpdg) != 413 & is_oldpcrm == 1 & %s", centr)},
        {"newpcrm_motherpi_neq413", Form("abs(light_motherpdg) != 413 & is_newpcrm == 1 & %s", centr)},
        {"oldpcrm_multDeq1", Form("heavy_mult == 1 & is_oldpcrm == 0 & %s", centr)},
        {"oldpcrm_multDgt1", Form("heavy_mult > 1 & is_oldpcrm == 0 & %s", centr)},
        {"newpcrm_multDeq1", Form("heavy_mult == 1 & is_newpcrm == 0 & %s", centr)},
        {"newpcrm_multDgt1", Form("heavy_mult > 1 & is_newpcrm == 0 & %s", centr)},
    };

    auto inFile = TFile::Open(inFileName.data());
    auto oFile = TFile::Open(oFileName.data(), "recreate");
    oFile->mkdir("distr");
    oFile->mkdir("checks");
    for (auto &[checkName, _]: checks) {
        oFile->mkdir(Form("checks/%s", checkName));
    }
    oFile->cd("distr");
    auto dir = (TDirectory *)inFile->Get(Form("HM_CharmFemto_%s_Trees0", pairLabel.data()));

    for (const char *event : {"SE", "ME"}) {
        for (const char *comb : {"pp", "mm", "pm", "mp"}) {
            for (const char *region : regions) {
                auto df = ROOT::RDataFrame(*(TTree *)dir->Get(Form("t%s_%s", event, comb)));

                // filter treefor standard analysis and systematics
                oFile->cd("distr");
                auto massSel = [hpdg, region](float mass, float pt) { return MassSelection(mass, pt, hpdg, region); };
                auto dfMass =
                    df.Filter(massSel, {"heavy_invmass", "heavy_pt"})
                        .Define("kStarMeV", "kStar * 1000")
                        .Define("light_pt", "pow(light_px*light_px + light_py*light_py, 0.5)")
                        .Define("light_nsigcomb", "pow(light_nsigtpc*light_nsigtpc + light_nsigtof*light_nsigtof, 0.5)")
                        .Filter(centr);


                const char *histName = Form("h%s_%s_%s", event, comb, region);
                const char *histTitle = ";#it{k}* (MeV/#it{c});Mult;Counts";
                auto hDistr = dfMass.Histo2D<float, int>({histName, histTitle, 3000u, 0., 3000., 180u, 0.5, 180.5},
                                                         "kStarMeV", "mult");
                hDistr->Write();

                // filter tree for additional checks
                oFile->cd("checks");
                for (auto &[checkName, checkSel] : checks) {
                    oFile->cd(Form("checks/%s", checkName));
                    std::cout<<checkName<<std::endl;

                    auto dfCheck = dfMass.Filter(checkSel.data());
                    auto hDistrCheck = dfCheck.Histo2D<float, int>(
                        {Form("%s_%s", histName, checkName), histTitle, 3000u, 0., 3000., 180u, 0.5, 180.5},
                        "kStarMeV", "mult");
                    hDistrCheck->Write(histName);
                }
            }
        }
    }
    oFile->Close();
    printf("Distributions saved to %s\n", oFileName.data());
}