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

    // load configuration for systematic variations
    YAML::Node config = YAML::LoadFile(cfgFileName.data());
    auto aliases = LoadAliases(config["aliases"]);
    auto selections = LoadSelections(config["selections"]);

    // open input file
    auto inFile = TFile::Open(inFileName.data());

    // open output file
    fs::path oFileName(suffix == "" ? "Distr.root" : Form("Distr_%s.root", suffix.data()));
    auto oFile = TFile::Open(std::string(oDir / oFileName).data(), "recreate");


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
                                              Form(";k* (MeV/#it{c});%s;Counts", heavy_mass_label.data()), 3000u, 0.,
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
                                                     Form(";k* (MeV/#it{c});%s;Counts", heavy_mass_label.data()), 3000u,
                                                     0., 3000., 1000u, charmMassMin, charmMassMax},
                                                    "kStarMeV", "heavy_invmass");
                    auto hMultVsKStar = dfSel.Histo2D<float, int>(
                        {Form("hMultVsKStar%lu", iSelection), ";#it{k}* (MeV/#it{c});Multiplicity;Counts", 3000u, 0.,
                         3000., 180u, 0.5, 180.5},
                        "kStarMeV", "mult");
                    auto hCharmMassVsPt =
                        dfSel.Histo2D<float, float>({Form("hCharmMassVsPt%lu", iSelection),
                                                     Form(";#it{p}_{T}(GeV/#it{c});%s;Counts", heavy_mass_label.data()),
                                                     100u, 0., 10., 1000u, charmMassMin, charmMassMax},
                                                    "heavy_pt", "heavy_invmass");

                    hCharmMassVsKStar->Write();
                    hMultVsKStar->Write();
                    hCharmMassVsPt->Write();
                }
            }
        }
    }
    oFile->Close();
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