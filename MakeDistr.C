#include <filesystem>
#include <map>
#include <string>

#include "ROOT/RDataFrame.hxx"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TTree.h"
#include "yaml-cpp/yaml.h"

namespace fs = std::filesystem;

template <typename T>
void print(std::vector<T>);

template <typename T>
std::vector<std::vector<T>> Combinations(std::vector<std::vector<T>>);

bool MassSelection(const double &mass, const double &pt, const int &hpdg, const std::string &massRegion,
                   const double fNSigmaMass = 2., double fNSigmaOffsetSideband = 5., double fSidebandWidth = 0.2,
                   double fLowerDstarRemoval = 1.992, double fUpperDstarRemoval = 2.028);

std::map<std::string, std::string> LoadAliases(YAML::Node);
std::vector<std::string> LoadSelections(YAML::Node);

void MakeDistr(std::string inFileName = "/data/DstarPi/tree_pc/mcgp/AnalysisResults_3998.root",
               std::string cfgFileName = "/home/daniel/phsw/fempy/selections.yml",
               fs::path oDir = "/home/daniel/an/DstarPi", std::string pair = "DstarPi",
               std::string suffix = "test",
               std::string treeDirName = "HM_CharmFemto_DstarPion_Trees0");

void MakeDistr(std::string inFileName, std::string cfgFileName, fs::path oDir, std::string pair, std::string suffix, std::string treeDirName) {
    // configure analysis
    float charmMassMin;
    float charmMassMax;
    std::vector<const char *> regions;
    std::string heavy_mass_label;
    int hpdg;

    if (pair == "DstarPi") {
        hpdg = 413;
        regions = {"sgn", "sbr"};
        heavy_mass_label  = "#it{M}(K#pi#pi) #minus #it{M}(K#pi) (GeV/#it{c})";
        charmMassMin = 0.14;
        charmMassMax = 0.24;
    } else if (pair == "DPi" || pair == "DK") {
        hpdg = 411;
        regions = {"sgn", "sbl", "sbr"};
        heavy_mass_label = "#it{M}(K#pi#pi) (GeV/#it{c})";
        charmMassMin = 1.7;
        charmMassMax = 2.0;
    } else {
        printf("\033[31mAnalysis not implemented. Exit!\033[0m\n");
        exit(1);
    }

    // For conversion to FD
    std::map<std::string, std::string> pairsToFD = {
        {"pp", "Particle0_Particle2"},
        {"mm", "Particle1_Particle3"},
        {"pm", "Particle1_Particle2"},
        {"mp", "Particle0_Particle3"},
        {"bs1", "Particle0_Particle0"},
        {"bs2", "Particle0_Particle1"},
        {"bs3", "Particle1_Particle1"},
        {"bs4", "Particle2_Particle2"},
        {"bs5", "Particle2_Particle3"},
        {"bs6", "Particle3_Particle3"},
    };
    std::map<std::string, std::string> regionToFD = {
        {"sgn", ""},
        {"sbl", "SBLeft"},
        {"sbr", "SBRight"},
    };
    // load configuration for systematic variations
    YAML::Node config = YAML::LoadFile(cfgFileName.data());
    auto aliases = LoadAliases(config["aliases"]);
    auto selections = LoadSelections(config["selections"]);

    // open input file
    auto inFile = TFile::Open(inFileName.data());

    // open output file
    fs::path oFileBaseName(suffix == "" ? "Distr.root" : Form("Distr_%s.root", suffix.data()));
    auto oFileName = std::string(oDir / oFileBaseName);
    auto oFile = TFile::Open(oFileName.data(), "recreate");

    std::vector<TList*> listFD; // one for each syst variation

    std::vector<std::map<std::string, TList*>> listPairsFD;
    for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++) {
        listFD.push_back(new TList());
        listPairsFD.push_back({
                {"pp", new TList()},
                {"mm", new TList()},
                {"pm", new TList()},
                {"mp", new TList()},
                {"bs1", new TList()},
                {"bs2", new TList()},
                {"bs3", new TList()},
                {"bs4", new TList()},
                {"bs5", new TList()},
                {"bs6", new TList()},
            });
        for (std::string comb : {"pp", "mm", "pm", "mp", "bs1", "bs2", "bs3", "bs4", "bs5", "bs6"}) {
            listPairsFD[iSelection][comb.data()]->SetName(pairsToFD[comb].data());
        }
    };

    // TList *myList = new TList();

    for (const char *comb : {"pp", "mm", "pm", "mp"}) {
        auto dir = (TDirectory *)inFile->Get(treeDirName.data());
        if (!dir) {
            std::cerr << "\033[31mError: Could not load the key \"" << treeDirName << "\". Available keys:\033[0m" << std::endl;
            inFile->ls();
            exit(1);
        }

        for (const char *event : {"SE", "ME"}) {
            auto tree = (TTree *)dir->Get(Form("t%s_%s", event, comb));
            ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager> df = ROOT::RDataFrame(*tree);

            // define aliases
            df = df.Define("kStarMeV", "kStar * 1000");
            for (auto &[name, expression] : aliases) {
                df = df.Define(name.data(), expression.data());
            }

            // project charm mass histograms
            oFile->mkdir(Form("%s/%s", comb, event));
            oFile->cd(Form("%s/%s", comb, event));
            for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++) {
                auto hCharmMassVsKStar =
                    df.Histo2D<float, float>({Form("hCharmMassVsKStar%lu", iSelection),
                                              Form(";#it{k}* (MeV/#it{c});%s;Counts", heavy_mass_label.data()), 3000u, 0.,
                                              3000., 1000u, charmMassMin, charmMassMax},
                                             "kStarMeV", "heavy_invmass");
                hCharmMassVsKStar->Write();
            }

            for (const char *region : regions) {
                std::cout << "\033[34mApplying selections to " << comb << "  " << event << "  " << region << "\033[0m"
                          << std::endl;
                oFile->mkdir(Form("%s/%s/%s", comb, event, region));
                oFile->cd(Form("%s/%s/%s", comb, event, region));

                auto lamMassSelection = [hpdg, region](float mass, float pt) {
                    return MassSelection(mass, pt, hpdg, region);
                };
                auto dfMass = df.Filter(lamMassSelection, {"heavy_invmass", "heavy_pt"});

                // apply selections
                for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++) {
                    std::cout << "selection " << iSelection << ": " << selections[iSelection] << std::endl;
                    auto dfSel = dfMass.Filter(selections[iSelection].data());

                    auto hCharmMassVsKStar =
                        dfSel.Histo2D<float, float>({Form("hCharmMassVsKStar%lu", iSelection),
                                                     Form(";#it{k}* (MeV/#it{c});%s;Counts", heavy_mass_label.data()), 3000u,
                                                     0., 3000., 1000u, charmMassMin, charmMassMax},
                                                    "kStarMeV", "heavy_invmass");
                    auto hMultVsKStar = dfSel.Histo2D<float, int>(
                        {Form("hMultVsKStar%lu", iSelection), ";#it{k}* (MeV/#it{c});Multiplicity;Counts", 3000u, 0.,
                         3000., 180u, 0.5, 180.5},
                        "kStarMeV", "mult");
                    auto hCharmMassVsPt =
                        dfSel.Histo2D<float, float>({Form("hCharmMassVsPt%lu", iSelection),
                                                     Form(";#it{p}_{T} (GeV/#it{c});%s;Counts", heavy_mass_label.data()),
                                                     100u, 0., 10., 1000u, charmMassMin, charmMassMax},
                                                    "heavy_pt", "heavy_invmass");

                    hCharmMassVsKStar->Write();
                    hMultVsKStar->Write();
                    hCharmMassVsPt->Write();

                    // // QA histograms
                    // dfSel.Histo1D<float>(
                    //     {Form("hHeavyBkgScore%lu", iSelection), ";Background score;Counts",
                    //      500u, 0, 0.1}, "heavy_bkg_score")->Write();
                    // dfSel.Histo2D<float, float>(
                    //     {Form("hLightNSigmaTOFVsNSigmaTPC%lu", iSelection), ";#it{n}_{#sigma}^{TPC};#it{n}_{#sigma}^{TOF};Counts",
                    //      200u, -10, 10, 200u, -10, 10}, "light_nsigtpc", "light_nsigtof")->Write();
                    // dfSel.Histo2D<double, float>(
                    //     {Form("hLightEtaVsPt%lu", iSelection), ";#it{p}_{T} (GeV/#it{c});#eta;Counts",
                    //      200u, -1, 1, 200u, 0, 5}, "light_pt", "light_eta")->Write();
                    // dfSel.Histo1D<int>(
                    //     {Form("hLightNCls%lu", iSelection), ";#it{n}_{clusters};Counts",
                    //      100, 59.5, 159.5}, "light_ncls")->Write();



                    // save as FD output

                    // myInnerList->Add(hMultVsKStar.GetPtr()->Clone(Form("h%s_%lu", event, iSelection)));


                    auto hMultVsKStarGeV = dfSel.Histo2D<float, int>(
                        {Form("hMultVsKStar%lu", iSelection), ";#it{k}* (GeV/#it{c});Multiplicity;Counts", 3000u, 0.,
                         3., 180u, 0.5, 180.5},
                        "kStar", "mult");
                    
                    listPairsFD[iSelection][comb]->Add(hMultVsKStarGeV.GetPtr()->ProjectionX()->Clone(Form("%sDist_%s", event, pairsToFD[comb].data())));
                    listPairsFD[iSelection][comb]->Add(hMultVsKStarGeV.GetPtr()->Clone(Form("%sMultDist_%s", event, pairsToFD[comb].data())));
                    // std::cout << listPairsFD[iSelection][comb]->GetSize() <<std::endl;
                    break;
                } // selections
                break;
            } // regions
            // break;
        } // SE, ME
        // break;
    } // charge comb

    // save output in FD format
    for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++){
        std::string dirNameFD = Form("HM_CharmFemto_DplusPion_Results%lu", iSelection);
        oFile->mkdir(dirNameFD.data());
        oFile->cd(dirNameFD.data());

        // save the data
        for (std::string comb : {"pp", "mm", "pm", "mp"}) {
            listFD[iSelection]->Add(listPairsFD[iSelection][comb]);
        }
        
        // give bullshit to GentleFemto
        for (std::string comb : {"bs1", "bs2", "bs3", "bs4", "bs5", "bs6"}) {
            for (const char *event : {"SE", "ME"}) {
                TH1D *hBS = new TH1D(Form("%sDist_%s", event, pairsToFD[comb].data()), "", 3000, 0, 3000);
                listPairsFD[iSelection][comb]->Add(hBS);
                TH1D *hMultBS = new TH1D(Form("%sMultDist_%s", event, pairsToFD[comb].data()), "", 3000, 0, 3000);
                listPairsFD[iSelection][comb]->Add(hMultBS);
            }
            listFD[iSelection]->Add(listPairsFD[iSelection][comb]);
        }
        
        // write to disk
        listFD[iSelection]->Write(dirNameFD.data(), 1);
    }

    oFile->Close();
    std::cout << "Output saved to: " << oFileName << std::endl;
}

// print vector
template <typename T>
void print(std::vector<T> vec) {
    for (auto elem : vec) std::cout << elem << "  ";
    std::cout << std::endl;
}

// all possible combinations of elements of sub-arrays
// arbitrary number of nested loops
template <typename T>
std::vector<std::vector<T>> Combinations(std::vector<std::vector<T>> items) {
    std::vector<std::vector<T>> combs = {};
    std::vector<long unsigned int> idxs = {};

    // initialize indeces
    for (const auto &item : items) idxs.push_back(0);

    // count how many loops must be done
    int nCombs = 1;
    for (auto list : items) nCombs *= list.size();

    // compute combinations
    int iComb = 0;
    do {
        std::vector<T> comb = {};
        for (long unsigned int iVar = 0; iVar < items.size(); iVar++) comb.push_back(items[iVar][idxs[iVar]]);
        combs.push_back(comb);

        int iVar = 0;
        while (idxs[iVar] == items[iVar].size() - 1) {  // reset idx when max is reached
            idxs[iVar++] = 0;
        }
        idxs[iVar]++;
    } while (++iComb < nCombs);

    return combs;
}

bool MassSelection(const double &mass, const double &pt, const int &hpdg, const std::string &massRegion,
                   const double fNSigmaMass, double fNSigmaOffsetSideband, double fSidebandWidth,
                   double fLowerDstarRemoval, double fUpperDstarRemoval) {
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

// load the aliases in a yaml node
std::map<std::string, std::string> LoadAliases(YAML::Node config) {
    std::map<std::string, std::string> aliases;
    for (YAML::iterator itAlias = config.begin(); itAlias != config.end(); ++itAlias) {
        for (auto alias : *itAlias) {
            aliases.insert({alias.first.as<std::string>(), alias.second.as<std::string>()});
        }
    }
    return aliases;
}

/*
Given a yaml node configuration with the systematic selections for each
variable, returns a vector of all the possible
combinations of the elementary selections.
*/
std::vector<std::string> LoadSelections(YAML::Node config) {
    std::vector<std::vector<std::string>> elemSelections = {};

    for (YAML::iterator itSel = config.begin(); itSel != config.end(); ++itSel) {
        auto var = (*itSel)["var_name"].as<std::string>();
        auto sign = (*itSel)["sign"].as<std::string>();
        auto values = (*itSel)["values"].as<std::vector<float>>();

        elemSelections.push_back({});
        for (auto value : values) {
            elemSelections.back().push_back(Form("%s %s %f", var.data(), sign.data(), float(value)));
        }
    }

    // concatenate the selections
    std::vector<std::string> tot_selections = {};
    for (auto sel : Combinations(elemSelections)) {
        std::string tot_sel = sel[0];
        for (long unsigned int iSel = 1; iSel < sel.size(); iSel++) tot_sel += Form(" && %s", sel[iSel].data());
        tot_selections.push_back(tot_sel);
    }
    return tot_selections;
}