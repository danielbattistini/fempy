#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"

void ComputeRawCF(TString inFileName = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/densities/Densities_test_kern-Gaussian_iter-Adaptive_mirror-NoMirror_binning-RelaxedBinning.root",
                  TString oFileName = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/RawCF.root") {
    TFile *inFile = new TFile(inFileName);

    double kStarMin = 0;
    double kStarMax = 3000;

    bool doUseReweighted = kTRUE;
    std::vector<TString> pairs = {"sc"};
    std::vector<TString> regions = {"sgn", "sbl", "sbr"};
    std::vector<TString> types = {"centr"};
    // std::vector<TString> types = {"centr", "upper", "lower"};

    TFile *oFile = new TFile(oFileName, "recreate");

    std::map<TString, TF1 *> mapCorrelationsUnc;
    std::map<TString, TF1 *> mapDensities;
    for (auto pairId : pairs) {
        for (auto regionId : regions) {
            double normFactor = 1;
            double normFactorRew = 1;
            TString baseName = Form("%s_%s", pairId.Data(), regionId.Data());
            std::cout << "analyzing " << pairId << "   " << regionId << "   " << baseName.Data() << std::endl;
            inFile->ls();
            // load se/me: central upper and lower
            for (auto type : types) {
                TString mapName = Form("%s_%s", baseName.Data(), type.Data());

                TF1 *fSE = (TF1 *)inFile->Get(Form("fSE_%s_%s", baseName.Data(), type.Data()));
                TF1 *fME = (TF1 *)inFile->Get(Form("fME_%s_%s", baseName.Data(), type.Data()));
                printf("get %s %s\n\n", type.Data(), Form("fMERew_%s_%s", baseName.Data(), type.Data()));
                TF1 *fMERew = (TF1 *)inFile->Get(Form("fMERew_%s_%s", baseName.Data(), type.Data()));
                std::cout << fSE->GetNpx() <<std::endl;
                std::cout << fME->GetNpx() <<std::endl;
                std::cout << fMERew->GetNpx() <<std::endl;

                mapDensities.insert({Form("SE_%s_%s", baseName.Data(), type.Data()), fSE});
                mapDensities.insert({Form("ME_%s_%s", baseName.Data(), type.Data()), fME});
                mapDensities.insert({Form("MERew_%s_%s", baseName.Data(), type.Data()), fMERew});

                if (type == "centr") {
                    normFactor = fME->Integral(1000, 1500, 1e-4) / fSE->Integral(1000, 1500, 1e-4);
                    normFactorRew = fMERew->Integral(1000, 1500, 1e-4) / fSE->Integral(1000, 1500, 1e-4);
                }

                TString cfName = Form("fCF_%s_%s", baseName.Data(), type.Data());
                auto cfFcn = [&](double *x, double *) {
                    return normFactor * fSE->EvalPar(x, nullptr) / fME->EvalPar(x, nullptr);
                };
                mapCorrelationsUnc.insert({mapName, new TF1(cfName, cfFcn, kStarMin, kStarMax, 0)});
                mapCorrelationsUnc[mapName]->SetTitle(";#it{k}* (MeV/#it{c});#it{C}");
                mapCorrelationsUnc[mapName]->SetNpx(1000);
                mapCorrelationsUnc[mapName]->Write();
                
                // rew cf
                TString cfRewName = Form("fCFRew_%s_%s", baseName.Data(), type.Data());
                auto cfRewFcn = [&](double *x, double *) {
                    return normFactorRew * fSE->EvalPar(x, nullptr) / fMERew->EvalPar(x, nullptr);
                };                
                mapCorrelationsUnc.insert({Form("Rew_%s", mapName.Data()), new TF1(cfRewName, cfRewFcn, kStarMin, kStarMax, 0)});
                mapCorrelationsUnc[Form("Rew_%s", mapName.Data())]->SetTitle(";#it{k}* (MeV/#it{c});#it{C}");
                mapCorrelationsUnc[Form("Rew_%s", mapName.Data())]->SetNpx(1000);
                mapCorrelationsUnc[Form("Rew_%s", mapName.Data())]->Write();

                // break;
            }
            std::cout << std::endl;

            // mapCorrelationsUnc[Form("%s_upper", baseName.Data())]->Draw();

            // auto cfAbsUncUpperFcn = [&](double *x, double *) {
            //     std::cout << mapCorrelationsUnc[Form("%s_upper", baseName.Data())]->EvalPar(x, nullptr) << "  "
            //               << mapCorrelationsUnc[Form("%s_centr", baseName.Data())]->EvalPar(x, nullptr) << std::endl;

            //     double upper = mapCorrelationsUnc[Form("%s_upper", baseName.Data())]->EvalPar(x, nullptr);
            //     double central = mapCorrelationsUnc[Form("%s_centr", baseName.Data())]->EvalPar(x, nullptr);
            //     // // return 3;
            //     // // return normFactor * fSE->EvalPar(x, nullptr) / fME->EvalPar(x, nullptr);

            //     return upper / central;
            //     // return 0.1;
            // };

            // auto name1 = "sc_sgn_upper";
            // auto name2 = "sc_sgn_centr";
            // // TCanvas *c = new TCanvas("c", "c", 600, 600);

            // double x[] = {20};
            // printf("%f%", mapCorrelationsUnc[name1]->Eval(20));
            // printf("%f%", mapCorrelationsUnc[name2]->Eval(20));
            // // inFile->Close();

            // for (auto p : mapCorrelationsUnc) {
            //     std::cout << p.first << "   " << p.second << std::endl;
            // }
            // oFile->Close();

            // std::cout << "ddd: " << baseName.Data() << "   " << mapCorrelationsUnc["sc_sbl_upper"] << std::endl;
            cout << "Corr func" << endl;

            for (const auto &[a, b] : mapCorrelationsUnc)
                std::cout << "map key: " << a << "   name: " << b << "   " << b->GetName() << std::endl;
            cout << " densities" << endl << endl;
            for (const auto &[a, b] : mapDensities)
                std::cout << "map key: " << a << "   name: " << b << "   " << b->GetName() << std::endl;

            // const char * aaa = "sc_sgn_lower";
            // std::cout << aaa <<  " dfvdfvdfvd  " << mapCorrelationsUnc[aaa]<<std::endl;;
            // continue;
            // upper unc

            // std::cout << "basenameeee: " << baseName.Data() << "   " << mapCorrelationsUnc["sc_sgn_lower"]->GetName() <<
            // std::endl;

            // double x[] = {20};
            // std::cout << "------->>> " << cfAbsUncUpperFcn(x, nullptr) << std::endl;
            // printf("%.10f", cfAbsUncUpperFcn(x, nullptr));
            // continue;
            // mapCorrelationsUnc[Form("%s_absunc_upper", baseName.Data())] =
            //     new TF1(Form("%s_unc_upper", baseName.Data()), cfAbsUncUpperFcn, kStarMin, kStarMax, 0);
            // mapCorrelationsUnc[Form("%s_absunc_upper", baseName.Data())]->Write();
            // std::cout << Form("%s_centr", baseName.Data()) << "    "
            //           << mapDensities[TString(Form("SE_sc_sgn_centr", baseName.Data()))] << std::endl;

            // auto cfRelUncUpperFcn = [&](double *x, double *) {
            //     // std::cout << Form("%s_centr", baseName.Data()) << "    "
            //     //           << mapDensities[Form("SE_sc_sgn_centr", baseName.Data())]->EvalPar(x, nullptr) << std::endl;
            //     if (!mapDensities[TString(Form("SE_%s_upper", baseName.Data()))])
            //         printf("---> %s is NULL", Form("SE_%s_upper", baseName.Data()));
            //     if (!mapDensities[TString(Form("SE_%s_centr", baseName.Data()))])
            //         printf("---> %s is NULL", Form("SE_%s_centr", baseName.Data()));

            //     double upper = mapDensities[TString(Form("SE_%s_upper", baseName.Data()))]->EvalPar(x, nullptr);
            //     double centr = mapDensities[TString(Form("SE_%s_centr", baseName.Data()))]->EvalPar(x, nullptr);
            //     return (upper - centr) / centr;
            // };
            // auto cfRelUncLowerFcn = [&](double *x, double *) {
            //     double lower = mapDensities[TString(Form("SE_%s_lower", baseName.Data()))]->EvalPar(x, nullptr);
            //     double centr = mapDensities[TString(Form("SE_%s_centr", baseName.Data()))]->EvalPar(x, nullptr);
            //     return (centr - lower) / centr;
            // };
            // auto fff = new TF1(Form("%s_relunc_upper", baseName.Data()), cfRelUncUpperFcn, kStarMin, kStarMax, 0);
            // auto fffL = new TF1(Form("%s_relunc_lower", baseName.Data()), cfRelUncLowerFcn, kStarMin, kStarMax, 0);
            // // // fff->Draw();
            // // mapCorrelationsUnc[mapName]->SetTitle(";#it{k}* (MeV/#it{c});#it{C}");
            // fff->SetNpx(1000);
            // fff->Write();
            // // mapCorrelationsUnc[mapName]->SetTitle(";#it{k}* (MeV/#it{c});#it{C}");
            // fffL->SetNpx(1000);
            // fffL->Write();
            // mapCorrelationsUnc.insert({Form("%s_relunc_upper", baseName.Data()), fff });

            // // mapCorrelationsUnc[Form("%s_relunc_upper", baseName.Data())] =
            // //     new TF1(Form("%s_unc_upper", baseName.Data()), cfRelUncUpperFcn, kStarMin, kStarMax, 0);
            // mapCorrelationsUnc[Form("%s_relunc_upper", baseName.Data())]->Write();

            // // lower unc
            // auto cfAbsUncLowerFcn = [&](double *x, double *) {
            //     double central = mapCorrelationsUnc[Form("%s_centr", baseName.Data())]->EvalPar(x, nullptr);
            //     double lower = mapCorrelationsUnc[Form("%s_lower", baseName.Data())]->EvalPar(x, nullptr);
            //     return central - lower;
            // };
            // mapCorrelationsUnc[Form("%s_absunc_lower", baseName.Data())] =
            //     new TF1(Form("%s_unc_lower", baseName.Data()), cfAbsUncLowerFcn, kStarMin, kStarMax, 0);
            // mapCorrelationsUnc[Form("%s_absunc_lower", baseName.Data())]->Write();

            //     auto cfRelUncLowerFcn = [&](double *x, double *) {
            //         double lower = mapCorrelationsUnc[Form("%s_lower", baseName.Data())]->EvalPar(x, nullptr);
            //         double central = mapCorrelationsUnc[Form("%s_centr", baseName.Data())]->EvalPar(x, nullptr);
            //         return (central - lower) / central;
            //     };
            //     mapCorrelationsUnc[Form("%s_relunc_lower", baseName.Data())] =
            //         new TF1(Form("%s_unc_lower", baseName.Data()), cfRelUncLowerFcn, kStarMin, kStarMax, 0);
            //     mapCorrelationsUnc[Form("%s_relunc_lower", baseName.Data())]->Write();
            break;
        }
    }
    oFile->Close();
    std::cout << "Raw CF were saved in " << oFileName << std::endl;
}
