#include <math.h> /* sqrt */

#include "yaml-cpp/yaml.h"

#include "Riostream.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TNamed.h"
#include "TSpline.h"
#include "TString.h"

double kStarRange[2] = {0, 3000};

YAML::Node lambdaParDB;

TString heavyTypes[] = {"true", "nonprompt", "Dstar", "combbkg"};
TString lightTypes[] = {"true", "secondaries", "strong", "misidentified"};
TString treatements[] = {"gen", "flat", "sideband", "coulomb"};
std::vector<TString> pairs = {"sc"};
std::vector<TString> types = {"sgn", "sbl", "sbr"};
// std::vector<TString> types = {"sgn", "sbl", "sbr", "mj", "fd"};
std::vector<TString> uncLabels = {"centr", "upper", "lower"};

std::map<TString, TNamed *> mapCF;
std::map<TString, double> summedLambdaPars;
TString currentPairId;

double wsbl = 0.51;
double wsbr = 1. - wsbl;

double GenCorrelationFcn(double *x, double *par) {
    double kStar = x[0];

    double csgn = ((TF1 *)mapCF[Form("%s_sgn_centr", currentPairId.Data())])->Eval(kStar);
    double csbl = ((TF1 *)mapCF[Form("%s_sbl_centr", currentPairId.Data())])->Eval(kStar);
    double csbr = ((TF1 *)mapCF[Form("%s_sbr_centr", currentPairId.Data())])->Eval(kStar);
    // double cfd = kStar < 900 ? ((TSpline3 *)mapCF[Form("%s_fd_centr", currentPairId.Data())])->Eval(kStar) : 1;

    double csb = wsbl * csbl + wsbr * csbr;

    double lsb = summedLambdaPars["sideband"];
    double ldstar = summedLambdaPars["coulomb"];
    double lgen = summedLambdaPars["gen"];
    double lflat = summedLambdaPars["flat"];

    // printf("~~~~~~~~~~~~~%f", lsb + ldstar + lgen + lflat);
    // double cf = (csgn - 0.3 * csb) / 0.7;
    double cf = (csgn - lsb * csb - lflat - ldstar)  / lgen;
    // double cf = ((csgn - lsb * csb) / 1 - ldstar * cfd - lflat) / lgen;
    return cf;
}

double GenCorrelationAbsUncFcn(double *x, double *par) {
    double kStar = x[0];
    auto name = Form("%s_sgn_absunc_upper", currentPairId.Data());
    printf("%s", name);

    double usgn = ((TF1 *)mapCF[Form("%s_sgn_absunc_upper", currentPairId.Data())])->Eval(kStar);
    double usbl = ((TF1 *)mapCF[Form("%s_sbl_absunc_upper", currentPairId.Data())])->Eval(kStar);
    double usbr = ((TF1 *)mapCF[Form("%s_sbr_absunc_upper", currentPairId.Data())])->Eval(kStar);
    // double cfd = kStar < 900 ? ((TSpline3 *)mapCF[Form("%s_fd_centr", currentPairId.Data())])->Eval(kStar) : 1;

    printf("%f\n", usgn);
    printf("%f\n", usbl);
    printf("%f\n\n", usbr);
    double usb = sqrt(wsbl * wsbl * usbl * usbl + wsbr * wsbr * usbr * usbr);

    double lsb = summedLambdaPars["sideband"];
    // double ldstar = summedLambdaPars["coulomb"];
    double lgen = summedLambdaPars["gen"];
    // double lflat = summedLambdaPars["flat"];

    // printf("~~~~~~~~~~~~~%f", lsb + ldstar + lgen + lflat);
    // double cf = (csgn - 0.3 * csb) / 0.7;
    // double cf = (csgn - lsb * csb - lflat - ldstar * cf)  / lgen;
    double unc = sqrt(usgn * usgn + lsb * lsb * usb * usb);
    // printf("%f", unc);
    return unc;
}

double GenCorrelationUpperUncFcn(double *x, double *par) {
    return GenCorrelationFcn(x, par) + GenCorrelationAbsUncFcn(x, par);
}

double GenCorrelationLowerUncFcn(double *x, double *par) {
    return GenCorrelationFcn(x, par) - GenCorrelationAbsUncFcn(x, par);
}
// double GenCorrelationUpperFcn(double *x, double *par) {
//     double kStar = x[0];
//     TStringcfId = Form("%s_%s_centr", pairId.Data(), type.Data());

//     for (auto pairId : pairs) {

//         // double csgn = ((TF1 *)mapCF[Form("%s_sgn", pairId.Data())])->Eval(kStar);
//         // double csbl = ((TF1 *)mapCF[Form("%s_sbl", pairId.Data())])->Eval(kStar);
//         // double csbr = ((TF1 *)mapCF[Form("%s_sbr", pairId.Data())])->Eval(kStar);
//         // double cfd = kStar < 900 ? ((TSpline3 *)mapCF[Form("%s_fd", pairId.Data())])->Eval(kStar) : 1;

//         double csgnUpper = ((TF1 *)mapCF[Form("%s_sgn_upper", pairId.Data())])->Eval(kStar);
//         double csblUpper = ((TF1 *)mapCF[Form("%s_sbl_upper", pairId.Data())])->Eval(kStar);
//         double csbrUpper = ((TF1 *)mapCF[Form("%s_sbr_upper", pairId.Data())])->Eval(kStar);

//         // double csgnUpperUnc = csgnUpper - csgn;
//         // double csblUpperUnc = csblUpper - csbl;
//         // double csbrUpperUnc = csbrUpper - csbr;

//         std::map<TString, double> summedLambdaPars;
//         for (auto t : treatements)
//             summedLambdaPars.insert({t, 0});

//         for (auto ht : heavyTypes)
//             for (auto lt : lightTypes) {
//                 auto treatment = lambdaParDB["Dpi_HMpp13TeV"][ht.Data()][lt.Data()].as<std::pair<double, TString>>();
//                 summedLambdaPars[treatment.second.Data()] += treatment.first;
//             }

//         double csb = 0.51 * csbl + 0.49 * csbr;

//         double lsb = summedLambdaPars["sideband"];
//         double ldstar = summedLambdaPars["coulomb"];
//         double lgen = summedLambdaPars["gen"];
//         double lflat = summedLambdaPars["flat"];

//         double cf = ((csgn - lsb * csb) / 1 - ldstar * cfd - lflat) / lgen;
//     }
//     return cf;
// }

void ComputeGenCF(TString inFileNameSgnSB = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/RawCF.root",
                  TString inFileNameFD = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/fd/Result_PipDp.root",
                  TString oFileName = "/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/GenCF.root",
                  TString lambdaParamDBFileName = "/home/daniel/alice/femtokde/utils/lambda_param_db.yml") {

    TFile *inFileSgnSB = new TFile(inFileNameSgnSB);
    TFile *inFileFD = new TFile(inFileNameFD);
    for (auto pairId : pairs) {
        for (auto type : types) {
            for (auto uncLabel : uncLabels) {
                TNamed *correlationFunction;
                TString cfId = Form("%s_%s_%s", pairId.Data(), type.Data(), uncLabel.Data());

                if (type == "sgn" || type == "sbl" || type == "sbr") {
                    correlationFunction = (TF1 *)inFileSgnSB->Get(Form("fCF_%s", cfId.Data()));
                    // todo: add minijets!
                } else if (type == "fd") {
                    auto gFD = (TGraphErrors *)inFileFD->Get("Dstar_Coulomb_smearedFULL1");
                    correlationFunction =
                        new TSpline3(gFD->GetName(), gFD, 0, gFD->GetPointY(0), gFD->GetPointY(gFD->GetN()));
                    correlationFunction->SetName(Form("splCF_%s", cfId.Data()));
                }
                mapCF.insert({cfId, correlationFunction});
            }
        }
    }

    inFileSgnSB->Close();
    inFileFD->Close();

    for (auto pairId : pairs) {
        for (auto type : types) {
            TString cfId = Form("%s_%s", pairId.Data(), type.Data());

            auto cfAbsUncUpperFcn = [&](double *x, double *) {
                double upper = 0;
                double centr = 0;

                if (type == "sgn" || type == "sbl" || type == "sbr") {
                    upper = ((TF1 *)mapCF[Form("%s_upper", cfId.Data())])->Eval(x[0]);
                    centr = ((TF1 *)mapCF[Form("%s_centr", cfId.Data())])->Eval(x[0]);
                    // } else if (type == "fd") {
                    //     centr = ((TSpline3 *)mapCF[Form("%s_centr", cfId.Data())])->Eval(x[0]);
                    //     upper = ((TSpline3 *)mapCF[Form("%s_upper", cfId.Data())])->Eval(x[0]);
                    // } else {
                    //     printf("Not implemented\n");
                    // printf("hahaha");
                }

                return upper - centr;
            };

            // double x[] = {20};
            // std::cout << cfId << "    " << cfAbsUncUpperFcn(x, nullptr) << std::endl;
            auto fCFAbsUncUpper =
                new TF1(Form("fCF_%s_absunc_upper", cfId.Data()), cfAbsUncUpperFcn, kStarRange[0], kStarRange[1], 0);
            mapCF.insert({Form("%s_absunc_upper", cfId.Data()), fCFAbsUncUpper});

            auto cfRelUncUpperFcn = [&](double *x, double *) {
                double delta = 0;
                double centr = 1;

                if (type == "sgn" || type == "sbl" || type == "sbr") {
                    delta = ((TF1 *)mapCF[Form("%s_absunc_upper", cfId.Data())])->Eval(x[0]);
                    centr = ((TF1 *)mapCF[Form("%s_centr", cfId.Data())])->Eval(x[0]);
                }

                return delta / centr;
            };
            auto fCFRelUncUpper =
                new TF1(Form("fCF_%s_relunc_upper", cfId.Data()), cfRelUncUpperFcn, kStarRange[0], kStarRange[1], 0);
            mapCF.insert({Form("%s_relunc_upper", cfId.Data()), fCFRelUncUpper});
        }
    }

    // first string is particle combinations, second is type of CF (signal, sb, feed-down etc.)
    // std::map<TString, std::map<TString, TF1 *>> correlations;

    // load correlation functions

    // for (auto pairId : pairs) {
    //     // central
    //     for (auto type : types) {
    //         TStringcfId = Form("%s_%s_centr", pairId.Data(), type.Data());

    //         if (type == "sgn" || type == "sbl" || type == "sbr") {
    //             TString name = Form("fCF_%s", cfId);
    //             mapCF.insert({cfId, (TF1 *)inFileSgnSB->Get(name.Data())});
    //             // todo: add minijets!
    //         } else if (type == "fd") {
    //             TStringname = "Dstar_Coulomb_smearedFULL1";
    //             auto gFD = (TGraphErrors *)inFileFD->Get(name);
    //             TSpline3 *splFD = new TSpline3(name, gFD, 0, gFD->GetPointY(0), gFD->GetPointY(gFD->GetN()));

    //             mapCF.insert({cfId, splFD});
    //         }
    //     }
    //     mapCF.insert({"centr", new TF1()});
    //     // mapCF.insert({"centr", new TF1(Form("fCF_%s_centr", pairId.Data()), GenCorrelationFcn, 0, 1000)});
    //     mapCF["centr"]->Write();

    //     // uncertainty
    //     for (auto type : types) {
    //         TStringcfId = Form("%s_%s_upper", pairId.Data(), type.Data());

    //         if (type == "sgn" || type == "sbl" || type == "sbr") {
    //             TString name = Form("fCF_%s_upper", cfId);
    //             mapCF.insert({cfId, (TF1 *)inFileSgnSB->Get(name.Data())});
    //             // todo: add minijets!
    //             // todo: add feed down uncertainty???
    //             // } else if (type == "fd") {
    //             //     TStringname = "Dstar_Coulomb_smearedFULL1";
    //             //     auto gFD = (TGraphErrors *)inFileFD->Get(name);
    //             //     TSpline3 *splFD = new TSpline3(name, gFD, 0, gFD->GetPointY(0), gFD->GetPointY(gFD->GetN()));

    //             //     mapCF.insert({cfId, splFD});
    //         }
    //     }
    //     // mapCF.insert({"upper", new TF1(Form("fCF_%s_upper", pairId.Data()), GenCorrelationUpperFcn, 0, 1000)});
    //     // mapCF.insert({"lower", new TF1()});
    //     // mapCF["lower"]->Write();

    //     for (auto type : types) {
    //         TStringcfId = Form("%s_%s_lower", pairId.Data(), type.Data());

    //         if (type == "sgn" || type == "sbl" || type == "sbr") {
    //             TString name = Form("fCF_%s_lower", cfId);
    //             mapCF.insert({cfId, (TF1 *)inFileSgnSB->Get(name.Data())});
    //             // todo: add minijets!
    //             // todo: add feed down uncertainty???
    //             // } else if (type == "fd") {
    //             //     TStringname = "Dstar_Coulomb_smearedFULL1";
    //             //     auto gFD = (TGraphErrors *)inFileFD->Get(name);
    //             //     TSpline3 *splFD = new TSpline3(name, gFD, 0, gFD->GetPointY(0), gFD->GetPointY(gFD->GetN()));

    //             //     mapCF.insert({cfId, splFD});
    //         }
    //     }
    //     // mapCF.insert({"lower", new TF1(Form("fCF_%s_lower", pairId.Data()), GenCorrelationUpperFcn, 0, 1000)});
    //     mapCF.insert({"lower", new TF1()});
    //     // mapCF["lower"]->Write();
    // }

    for (auto t : treatements)
        summedLambdaPars.insert({t, 0});

    lambdaParDB = YAML::LoadFile(lambdaParamDBFileName.Data());
    for (auto ht : heavyTypes) {
        for (auto lt : lightTypes) {
            auto treatment = lambdaParDB["Dpi_HMpp13TeV"][ht.Data()][lt.Data()].as<std::pair<double, TString>>();
            summedLambdaPars[treatment.second.Data()] += treatment.first;
        }
    }

    TFile *oFile = new TFile(oFileName, "recreate");
    oFile->cd();
    // comppute the gen cf
    TF1 *cf;
    TF1 *cfUnc;
    TF1 *cfUncUpper;
    TF1 *cfUncLower;

    for (auto pairId : pairs) {
        currentPairId = pairId;

        cf = new TF1("cf", GenCorrelationFcn, 0, 3000, 0);
        cfUnc = new TF1("cfUnc", GenCorrelationAbsUncFcn, 0, 3000, 0);
        cfUncUpper = new TF1("cfUncUpper", GenCorrelationUpperUncFcn, 0, 3000, 0);
        cfUncLower = new TF1("cfUncLower", GenCorrelationLowerUncFcn, 0, 3000, 0);
        cf->Write();
        cfUnc->Write();
        cfUncUpper->Write();
        cfUncLower->Write();
    }

    for (const auto &[a, b] : mapCF) {
        std::cout << a << "   " << b << std::endl;
        b->Write();
    }

    // for (const auto &[a, b] : summedLambdaPars) {

    //     std::cout << a << "   " << b << std::endl;
    // }

    cfUnc->Draw();
    new TBrowser();
    // oFile->Close();
}