#include <filesystem>
#include <map>
#include <string>
#include <algorithm>

#include "ROOT/RDataFrame.hxx"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TTree.h"
#include "yaml-cpp/yaml.h"

namespace fs = std::filesystem;

template <typename T>
void print(std::vector<T>, std::string="\n");

template <typename T>
std::vector<std::vector<T>> Combinations(std::vector<std::vector<T>>);

bool IsColInDF(ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void> df, std::string colName) {
    for (auto &&dfColName : df.GetColumnNames()) if (std::string(dfColName)==colName) return true;
    return false;
}

bool MassSelection(const double &mass, const double &pt, const int &hpdg, const std::string &massRegion,
                   const double fNSigmaMass = 2., double fNSigmaOffsetSideband = 5., double fSidebandWidth = 0.2,
                   double fLowerDstarRemoval = 1.992, double fUpperDstarRemoval = 2.028);

std::map<std::string, std::string> LoadAliases(YAML::Node);
std::vector<std::string> LoadSelections(YAML::Node, int n=0);

void MakeDistr(
    std::string inFileName = "/data/DstarPi/tree_pc/mcgp/AnalysisResults_3998.root",
    std::string cfgFileName = "/home/daniel/an/DPi/cfg_selection_nosel.yml",
    bool isMC = false,
    fs::path oDir = "/home/daniel/an/DstarPi",
    std::string pair = "DstarPi",
    std::string suffix = "test",
    std::string treeDirName = "auto",
    double nSigmaMass = 2.,
    double nSigmaOffsetSideband = 5., 
    double sidebandWidth = 0.2,
    bool uniq = false,
    bool doKStarSlices = false,
    double kStarBW = 50, // MeV/c
    bool doQA=false,
    bool doSyst = false,
    unsigned int nSelToKeep = 0,
    int nJobs=1);

void MakeDistr(
    std::string inFileName,
    std::string cfgFileName,
    bool isMC,
    fs::path oDir,
    std::string pair,
    std::string suffix,
    std::string treeDirName,
    double nSigmaMass,
    double nSigmaOffsetSideband, 
    double sidebandWidth,
    bool uniq,
    bool doKStarSlices,
    double kStarBW,
    bool doQA,
    bool doSyst,
    unsigned int nSelToKeep,
    int nJobs) {
    srand(42);

    if (nJobs>1 && ! uniq)
        ROOT::EnableImplicitMT(nJobs); 

    // configure analysis
    float charmMassMin;
    float charmMassMax;
    std::vector<const char *> regions;
    std::string heavy_mass_label;
    int hpdg;

    if (pair == "DstarPi" || pair == "DstarK") {
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

    //todo: switch to regions as parameter and remove isMC
    if (isMC) regions = {"sgn"};
    int lpdg;
    if (pair == "DstarPi" || pair == "DPi") {
        lpdg = 211;
    } else if (pair == "DK" || pair == "DstarK") {
        lpdg = 321;
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
        {"sbl", "SBLeft_"},
        {"sbr", "SBRight_"},
    };
    // load configuration for systematic variations
    YAML::Node config = YAML::LoadFile(cfgFileName.data());
    auto aliases = LoadAliases(config["aliases"]);
    auto selections = LoadSelections(config["selections"], nSelToKeep);
    if (!doSyst) selections = {selections[0]};
    // open input file
    auto inFile = TFile::Open(inFileName.data());

    // open output file
    fs::path oFileBaseName(suffix == "" ? "Distr.root" : Form("Distr_%s.root", suffix.data()));
    auto oFileName = std::string(oDir / oFileBaseName);
    auto oFile = TFile::Open(oFileName.data(), "recreate");

    std::string logFileName = Form("%s.log", oFileName.substr(0, oFileName.find_last_of(".")).data()); 
    ofstream log;
    log.open(logFileName.data());

    std::map<std::string, std::vector<TList*>> listFD; // one for each syst variation

    std::vector<std::map<std::string, std::map<std::string, TList*>>> listPairsFD;
    for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++) {
        for (const auto region : regions)
            listFD[region].push_back(new TList());
        
        std::map<std::string, std::map<std::string, TList*>> mapLists;
        for (const auto region : regions) {
            for (std::string comb : {"pp", "mm", "pm", "mp", "bs1", "bs2", "bs3", "bs4", "bs5", "bs6"}) {
                mapLists[region][comb] = new TList();
            }
        }
        
        listPairsFD.push_back(mapLists);
        for (const auto region : regions) {
            for (std::string comb : {"pp", "mm", "pm", "mp", "bs1", "bs2", "bs3", "bs4", "bs5", "bs6"}) {
                listPairsFD[iSelection][region][comb.data()]->SetName(pairsToFD[comb].data());
            }
        }
    };

    // autodetect treeDirname
    if (treeDirName == "auto") {
        std::vector<std::string> treeDirNames;
        for (auto key : *inFile->GetListOfKeys()) {
            if (std::string(key->GetName()).find("_Trees") != std::string::npos)
                treeDirNames.push_back(std::string(key->GetName()));
        }
        if (treeDirNames.size() > 1) {
            printf("\033[31mError: the input file contains more than 1 directory with trees. "
                   "Split the file with rootmv. Content of the file:\033[0m]");
            inFile->ls();
            return;
        }
        if (treeDirNames.size() == 0) {
            printf("\033[31mError: the input doesn't contain any directory with trees. Content of the file:\033[0m]");
            inFile->ls();
            return;
        }
        treeDirName = treeDirNames[0];
    }
    // remove the HM_CharmFemto_ (14 characters) and the suffix _TreesX
    static std::string pairName = treeDirName.substr(14, treeDirName.find("_Trees")-14);
    
    static std::string name, title;
    using namespace ROOT::RDF;

    for (const char *comb : {"pp", "mm", "pm", "mp"}) {
        auto dir = (TDirectory *)inFile->Get(treeDirName.data());
        if (!dir) {
            std::cerr << "\033[31mError: Could not load the key \"" << treeDirName << "\". Available keys:\033[0m" << std::endl;
            inFile->ls();
            exit(1);
        }

        for (std::string event : {"SE", "ME"}) {
            auto tree = (TTree *)dir->Get(Form("t%s_%s", event.data(), comb));
            if (!tree) {
                std::cerr << "Error: tree not loaded properly. Exit!" << std::endl;
                exit(1);
            }
            if (tree->GetEntries() == 0) {
                std::cerr << "Error: the tree has 0 entries. Exit!" << std::endl;
                exit(1);
            }
            ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager> df = ROOT::RDataFrame(*tree);

            // define aliases
            df = df.Define("kStarMeV", "kStar * 1000");
            for (auto &[name, expression] : aliases) {
                df = df.Define(name.data(), expression.data());
            }

            // project charm mass histograms
            oFile->mkdir(Form("%s/%s", comb, event.data()));
            oFile->cd(Form("%s/%s", comb, event.data()));
            
            for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++) {
                if (!doKStarSlices) break;
                auto dfSel = df.Filter(selections[iSelection].data());

                name = Form("hCharmMassVsKStar%lu", iSelection);
                title = Form(";#it{k}* (MeV/#it{c});%s;Counts", heavy_mass_label.data());
                TH2DModel model = {name.data(), title.data(), 3000u, 0., 3000., 1000u, charmMassMin, charmMassMax};
                auto hCharmMassVsKStar = dfSel.Histo2D<float, float>(model, "kStarMeV", "heavy_invmass");
                hCharmMassVsKStar->Write();
                int nKStarBins = round(3000/kStarBW);
                for (int iKStarBin = 0; iKStarBin < nKStarBins; iKStarBin++){
                    double kStarMin = kStarBW * iKStarBin;
                    double kStarMax = kStarBW * (iKStarBin + 1);

                    int firstBin = hCharmMassVsKStar->GetXaxis()->FindBin(kStarMin * 1.0001);
                    int lastBin = hCharmMassVsKStar->GetXaxis()->FindBin(kStarMax * 0.9999);

                    name = Form("hCharmMass%lu_kStar%.0f_%.0f", iSelection, kStarMin, kStarMax);
                    auto hCharmMassKStarBin = hCharmMassVsKStar->ProjectionY(name.data(), firstBin, lastBin);
                    hCharmMassKStarBin->SetTitle(Form(";%s;Counts", heavy_mass_label.data()));
                    hCharmMassKStarBin->Write();
                }
            }

            for (const char *region : regions) {
                std::cout << "\033[34mApplying selections to "
                          << comb << "  "
                          << event << "  "
                          << region << "\033[0m"
                          << std::endl;
                log << "Applying selections to "
                    << comb << "  "
                    << event << "  "
                    << region
                    << std::endl;

                oFile->mkdir(Form("%s/%s/%s", comb, event.data(), region));
                oFile->cd(Form("%s/%s/%s", comb, event.data(), region));

                auto lamMassSelection = [hpdg, region, nSigmaMass, nSigmaOffsetSideband, sidebandWidth](float mass, float pt) {
                    return MassSelection(mass, pt, hpdg, region, nSigmaMass, nSigmaOffsetSideband, sidebandWidth);
                };
                auto dfMass = df.Filter(lamMassSelection, {"heavy_invmass", "heavy_pt"});

                // apply selections
                for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++) {
                    std::cout << "selection " << iSelection << ": " << selections[iSelection] << std::endl;
                    log << "selection " << iSelection << ": " << selections[iSelection] << std::endl;

                    auto dfSel = dfMass.Filter(selections[iSelection].data());

                    auto hCharmPtVsKStar = dfSel.Histo2D<float, float>({"hCharmPtVsKStar0",
                                                        ";#it{k}*;#it{p}_{T} (GeV/#it{c});Counts",
                                                        3000u, 0., 3000., 100u, 0, 10},
                                                        "kStarMeV", "heavy_pt");
                    hCharmPtVsKStar->Write();

                    if (doKStarSlices) {
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

                        if (uniq && event == "SE") {
                            static std::vector<double> seenDouble = {};
                            auto EraseDoubleDuplicates = [](double var) {
                                for (const auto &seenMass : seenDouble) {
                                    if (std::abs(seenMass - var)/var < 0.000000001) {
                                        return false;
                                    }
                                }
                                seenDouble.insert(seenDouble.begin(), var); // insert in the front saves time
                                return true;
                            };
                            
                            std::string name, title;

                            // uniq light pt
                            title = ";#it{p}_{T} (GeV/#it{c});Counts";
                            name = Form("hLightPtUniq%lu_kStar0_3000", iSelection);
                            dfSel.Filter(EraseDoubleDuplicates, {"light_pt"})
                                .Histo1D<double>({name.data(), title.data(), 200u, 0, 5}, "light_pt")
                                ->Write();
                            seenDouble.clear();

                            if (IsColInDF(dfSel, "light_pdg")) {
                                // uniq true light pt
                                name = Form("hTrueLightPtUniq%lu_kStar0_3000", iSelection);
                                title = ";#it{p}_{T} (GeV/#it{c});Counts";
                                dfSel.Filter(EraseDoubleDuplicates, {"light_pt"})
                                    .Filter(Form("light_pdg == %d", comb[1] == 'p' ? +lpdg : -lpdg))
                                    .Histo1D<double>({name.data(), title.data(), 200u, 0, 5}, "light_pt")
                                    ->Write();
                                seenDouble.clear();
                            }
                        }
                    }

                    // save as FD output
                    auto hMultVsKStarGeV = dfSel.Histo2D<float, int>(
                        {Form("hMultVsKStar%lu", iSelection), ";#it{k}* (GeV/#it{c});Multiplicity;Counts", 600u, 0.,
                         3., 180u, 0.5, 180.5},
                        "kStar", "mult");
                    
                    listPairsFD[iSelection][region][comb]->Add(hMultVsKStarGeV.GetPtr()->ProjectionX()->Clone(Form("%sDist_%s", event.data(), pairsToFD[comb].data())));
                    listPairsFD[iSelection][region][comb]->Add(hMultVsKStarGeV.GetPtr()->Clone(Form("%sMultDist_%s", event.data(), pairsToFD[comb].data())));

                    // QA histograms
                    if (doQA && IsColInDF(dfSel, "heavy_bkg_score")) {
                        dfSel.Histo2D<float, float>(
                            {Form("hHeavyBkgScoreVsPt%lu", iSelection), ";#it{p}_{T} (GeV/#it{c});Background score;Counts",
                            100u, 0., 10., 500u, 0, 0.1}, "heavy_pt", "heavy_bkg_score")->Write();
                    }
                            
                    if (doQA && IsColInDF(dfSel, "light_nsigtpc") && IsColInDF(dfSel, "light_nsigtof")) {
                        dfSel.Histo2D<float, float>(
                            {Form("hLightNSigmaTOFVsNSigmaTPC%lu", iSelection), ";#it{n}_{#sigma}^{TPC};#it{n}_{#sigma}^{TOF};Counts",
                            200u, -10, 10, 200u, -10, 10}, "light_nsigtpc", "light_nsigtof")->Write();
                    }

                    if (doQA && IsColInDF(dfSel, "light_pt") && IsColInDF(dfSel, "light_eta")) {
                        auto hLightEtaVsPt = dfSel.Histo2D<double, float>(
                            {Form("hLightEtaVsPt%lu", iSelection), ";#it{p}_{T} (GeV/#it{c});#eta;Counts",
                            200u, 0, 5, 200u, -1, 1}, "light_pt", "light_eta");
                        hLightEtaVsPt->SetTitle(";#it{p}_{T} (GeV/#it{c});Counts");
                        hLightEtaVsPt->Write();
                        hLightEtaVsPt->ProjectionX(Form("hLightPt%lu", iSelection))->Write();;
                    }

                    if (doQA && IsColInDF(dfSel, "light_ncls")) {
                        dfSel.Histo1D<int>(
                            {Form("hLightNCls%lu", iSelection), ";#it{n}_{clusters};Counts",
                            100, 59.5, 159.5}, "light_ncls")->Write();
                    }

                    if (doQA && IsColInDF(dfSel, "light_pdg")) {
                        dfSel.Filter(Form("light_pdg == %d", comb[1] == 'p' ? +lpdg : -lpdg)).Histo1D<double>(
                            {Form("hTrueLightPt%lu", iSelection), ";#it{p}_{T} (GeV/#it{c});Counts",
                            200u, 0, 5}, "light_pt")->Write();
                    }
                } // selections
            } // regions
        } // SE, ME
    } // charge comb

    // save output in FD format
    for (long unsigned int iSelection = 0; iSelection < selections.size(); iSelection++){
        for (std::string region : regions) {
            std::string dirNameFD = Form("HM_CharmFemto_%s_%sResults%lu", pairName.data(), regionToFD[region].data(), iSelection);
            oFile->mkdir(dirNameFD.data());
            oFile->cd(dirNameFD.data());

            // save the data
            for (std::string comb : {"pp", "mm", "pm", "mp"}) {
                listFD[region][iSelection]->Add(listPairsFD[iSelection][region][comb]);
            }
            
            // give bullshit to GentleFemto
            for (std::string comb : {"bs1", "bs2", "bs3", "bs4", "bs5", "bs6"}) {
                for (const char *event : {"SE", "ME"}) {
                    TH1D *hBS = new TH1D(Form("%sDist_%s", event, pairsToFD[comb].data()), "", 3000, 0, 3000);
                    listPairsFD[iSelection][region][comb]->Add(hBS);
                    TH1D *hMultBS = new TH1D(Form("%sMultDist_%s", event, pairsToFD[comb].data()), "", 3000, 0, 3000);
                    listPairsFD[iSelection][region][comb]->Add(hMultBS);
                }
                listFD[region][iSelection]->Add(listPairsFD[iSelection][region][comb]);
            }
            listFD[region][iSelection]->Write(dirNameFD.data(), 1);
        }
    }

    oFile->Close();
    log.close();

    std::cout << "Output saved to: " << oFileName << std::endl;
}

// print vector
template <typename T>
void print(std::vector<T> vec, std::string separator) {
    for (auto elem : vec) std::cout << elem << separator;
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
std::vector<std::string> LoadSelections(YAML::Node config, int n) {
    std::vector<std::vector<std::string>> elemSelections = {};

    for (YAML::iterator itSel = config.begin(); itSel != config.end(); ++itSel) {
        auto var = (*itSel)["var_name"].as<std::string>();
        auto sign = (*itSel)["sign"].as<std::string>();
        auto values = (*itSel)["values"].as<std::vector<std::string>>();

        elemSelections.push_back({});
        for (auto value : values) {
            elemSelections.back().push_back(Form("%s %s %s", var.data(), sign.data(), value.data()));
        }
    }

    // concatenate the selections
    std::vector<std::string> tot_selections = {};
    for (auto sel : Combinations(elemSelections)) {
        std::string tot_sel = sel[0];
        for (long unsigned int iSel = 1; iSel < sel.size(); iSel++) tot_sel += Form(" && %s", sel[iSel].data());
        tot_selections.push_back(tot_sel);
    }
    if (n<0) return tot_selections;

    // reduce the number of selections to n
    int idx;
    std::vector<int> indeces = {};
    std::vector<std::string> red_selections = {tot_selections[0]};
    for (unsigned long int iSel = 1; iSel < n && iSel < tot_selections.size(); iSel++) {
        do {
            idx = rand()%(tot_selections.size()-1)+1; // always keep the central selection (idx=0)
        } while (std::find(indeces.begin(), indeces.end(), idx) != indeces.end());
        indeces.push_back(idx);
        red_selections.push_back(tot_selections[idx]);
    }
    return red_selections;
}