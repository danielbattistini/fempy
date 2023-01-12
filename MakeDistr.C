// #include <map>
// #include <string>

// #include "ROOT/RDataFrame.hxx"
// #include "Riostream.h"
// #include "TDatabasePDG.h"
// #include "TFile.h"
// #include "TH1.h"
// #include "TH2.h"
// #include "TSystem.h"
// #include "TTree.h"
// // #include "fempy/utils/Analysis.hxx"
// #include "yaml-cpp/yaml.h"

// void MakeDistr(std::string inFileName = "/data/DstarPi/tree_pc/mcgp/AnalysisResults_3998.root",
//                std::string oFileName = "/home/daniel/an/DstarPi/test_proj_syst.root",
//                std::string pair = "DstarPi");

void MakeDistr(std::string inFileName = "/data/DstarPi/tree_pc/mcgp/AnalysisResults_3998.root",
               std::string cfgFileName = "/home/daniel/phsw/fempy/selections.yml",
               std::string oFileName = "/home/daniel/an/DstarPi/test_proj_syst.root", std::string pair = "DstarPi");


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
    std::vector<std::vector<T>> combs = {{}};
    std::vector<int> idxs = {};

    // initialize indeces
    for (int i; i < items.size(); i++) idxs.push_back(0);

    // count how many loops must be done
    int nCombs = 1;
    for (auto list : items) nCombs *= list.size();

    // compute combinations
    int iComb = 0;
    do {
        std::vector<T> comb = {};
        for (int iVar = 0; iVar < items.size(); iVar++) comb.push_back(items[iVar][idxs[iVar]]);
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

// void MakeDistr(std::string inFileName, std::string oFileName, std::string pair) {
std::map<std::string, std::string> LoadAliases(YAML::Node config) {
    std::map<std::string, std::string> aliases;
    for (YAML::iterator itAlias = config.begin(); itAlias != config.end(); ++itAlias) {
        for (auto alias : *itAlias) {
            aliases.insert({alias.first.as<std::string>(), alias.second.as<std::string>()});
        }
    }
    return aliases;
}



std::vector<std::string> LoadSelections(YAML::Node config) {
    std::vector<std::vector<std::string>> elemSelections = {};

    for (YAML::iterator itSel = config.begin(); itSel != config.end(); ++itSel) {
        auto var = (*itSel)["var_name"].as<std::string>();
        auto sign = (*itSel)["sign"].as<std::string>();
        auto values = (*itSel)["values"].as<std::vector<float>>();

        // std::cout<< values[0]<< std::endl;

        elemSelections.push_back({});
        for (auto value : values) {
            elemSelections.back().push_back(Form("%s %s %f", var.data(), sign.data(), float(value)));
        }
        // for (auto alias : *itSel) {
        //     aliases.insert({alias.first.as<std::string>(), alias.second.as<std::string>()});
        // }
    }
    // std::vector<std::string> selections = {};
    // for (auto elemSel : elemSelections) std::string selection;

    // for (auto sel : elemSel) std::cout << sel << std::endl;
    // std::cout << std::endl;

    // selections.push_back(selection);

    // // combine elementary selections

    return {};  // Combinations(elemSelections);
}

void MakeDistr(std::string inFileName, std::string cfgFileName, std::string oFileName, std::string pair) {
    // std::vector<std::vector<int>> elems = {{1, 2, 3}, {4, 5}, {6, 7, 8}};
    // std::vector<std::vector<int>> v = Combinations(elems);
    // for (auto vv : v) print(vv);
    // exit(1);

    YAML::Node config = YAML::LoadFile(cfgFileName.data());

    auto inFile = TFile::Open(inFileName.data());
    auto oFile = TFile::Open(oFileName.data(), "recreate");

    auto aliases = LoadAliases(config["aliases"]);
    auto selections = LoadSelections(config["selections"]);

    std::vector<const char *> regions = {"sgn", "sbl", "sbr"};

    // for (auto al : selections)
    //     std::cout<<al<< std::endl;
    std::string heavy_mass_label = "#it{M}(K#pi#pi) #minus #it{M}(K#pi) (GeV/#it{c})";
    const float charmMassMin = 0.14;
    const float charmMassMax = 0.24;

    for (const char *comb : {"pp", "mm", "pm", "mp"}) {
        // oFile->mkdir(comb);
        // oFile->cd(comb);

        // oFile->ls();
        auto dir = (TDirectory *)inFile->Get("HM_CharmFemto_DstarPion_Trees0");

        for (const char *event : {"SE", "ME"}) {
            auto tree = (TTree *)dir->Get(Form("t%s_%s", event, comb));
            ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager> df = ROOT::RDataFrame(*tree);

            // define aliases
            df = df.Define("kStarMeV", "kStar * 1000");
            for (auto &[name, expression] : aliases) {
                df = df.Define(name.data(), expression.data());
            }

            for (const char *region : regions) {
                oFile->mkdir(Form("%s/%s/%s", comb, region, event));
                oFile->cd(Form("%s/%s/%s", comb, region, event));

                int hpdg = 413;
                auto lamMassSelection = [hpdg, region](float mass, float pt) {
                    return MassSelection(mass, pt, hpdg, region);
                };
                auto dfMass = df.Filter(lamMassSelection, {"heavy_invmass", "heavy_pt"});

                // apply selections
                auto dfSel = dfMass.Filter("kStarMeV>300");

                auto hCharmMassVsKStar = dfSel.Histo2D<float, float>(
                    {"hCharmMassVsKStar", Form(";k* (MeV/#it{c});%s;Counts", heavy_mass_label.data()), 3000u, 0., 3000.,
                     1000u, charmMassMin, charmMassMax},
                    "kStarMeV", "heavy_invmass");
                auto hMultVsKStar = dfSel.Histo2D<float, int>(
                    {"hMultVsKStar", ";#it{k}* (MeV/#it{c});Multiplicity;Counts", 3000u, 0., 3000., 180u, 0.5, 180.5},
                    "kStarMeV", "mult");
                auto hCharmMassVsPt = dfSel.Histo2D<float, float>(
                    {"hCharmMassVsPt", Form(";#it{p}_{T}(GeV/#it{c});%s;Counts", heavy_mass_label.data()), 100u, 0.,
                     10., 1000u, charmMassMin, charmMassMax},
                    "heavy_pt", "heavy_invmass");

                hCharmMassVsKStar->Write();
                hMultVsKStar->Write();
                hCharmMassVsPt->Write();
            }
            break;
        }
    }
    oFile->Close();
}