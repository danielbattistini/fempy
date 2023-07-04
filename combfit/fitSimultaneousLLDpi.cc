#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Fit/Fitter.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>

#include "./functions.h"

// definition of shared parameter
// same charge function
const int nParsSC = 8;
int iparSC[nParsSC] = {
    0,  // first source radius
    1,  // second source radius
    2,  // normalisation of two contributions
    3,  // scattering length singlet
    4,  // effective range singlet
    7,  // QS
    8,  // charge product
    10  // reduced mass
};

// opposite charge function
const int nParsOC = 10;
int iparOC[nParsOC] = {
    0,  // first source radius
    1,  // second source radius
    2,  // normalisation of two contributions
    3,  // scattering length singlet
    4,  // effective range singlet
    5,  // scattering length triplet
    6,  // effective range triplet
    7,  // QS
    9,  // charge product
    10  // reduced mass
};

// Create the GlobalCHi2 structure
struct GlobalChi2 {
    GlobalChi2(ROOT::Math::IMultiGenFunction &f1, ROOT::Math::IMultiGenFunction &f2) : fChi2_1(&f1), fChi2_2(&f2) {}

    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator()(const double *par) const {
        double p1[nParsSC];  // NOLINT [runtime/arrays]
        for (int i = 0; i < nParsSC; ++i) {
            p1[i] = par[iparSC[i]];
        }
        double p2[nParsOC];  // NOLINT [runtime/arrays]
        for (int i = 0; i < nParsOC; ++i) {
            p2[i] = par[iparOC[i]];
        }

        return (*fChi2_1)(p1) + (*fChi2_2)(p2);
    }

    const ROOT::Math::IMultiGenFunction *fChi2_1;
    const ROOT::Math::IMultiGenFunction *fChi2_2;
};

std::array<double, 2> ComputeChi2(TGraphAsymmErrors *graphFirst, TGraphAsymmErrors *graphSecond, TF1 *funcFirst,
                                  TF1 *funcSecond, double kStarMax = 200., int npars = 2) {
    double chi2{0.}, ndf{0};
    for (int iPt{0}; iPt < graphFirst->GetN(); ++iPt) {
        double kStar{-1.}, cfFirst{-1.}, cfSecond{-1.};
        graphFirst->GetPoint(iPt, kStar, cfFirst);
        if (kStar <= kStarMax) {
            graphSecond->GetPoint(iPt, kStar, cfSecond);
            double deltaFirst = cfFirst - funcFirst->Eval(kStar);
            double deltaSecond = cfSecond - funcSecond->Eval(kStar);
            double uncFirst = graphFirst->GetErrorYhigh(iPt);
            double uncSecond = graphSecond->GetErrorYhigh(iPt);
            chi2 += deltaFirst * deltaFirst / (uncFirst * uncFirst);
            chi2 += deltaSecond * deltaSecond / (uncSecond * uncSecond);
            ndf += 2;
        } else {
            break;
        }
    }
    return std::array<double, 2>{chi2, ndf - npars};
}

double GetMinHisto(TH1 *histo) {
    for (int iBin{1}; iBin <= histo->GetNbinsX(); ++iBin) {
        if (histo->GetBinContent(iBin) > 0) {
            return histo->GetBinCenter(iBin);
        }
    }
    return -1;
}

ROOT::Fit::FitResult fitSimultaneousLL(TGraphAsymmErrors *gSC, TGraphAsymmErrors *gOC, double radiusFirst = 0.97,
                                       double radiusSecond = 2.52, double weightFirst = 0.66, double maxRange = 300.,
                                       bool deleteFile = false, std::string outName = "LLSimFit_Dpi",
                                       std::string suffix = "", bool savePDF = false, bool fixRadii = true) {
    gStyle->SetPadRightMargin(0.035);
    gStyle->SetPadTopMargin(0.035);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetTitleSize(0.045, "xy");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    auto kAzureMy = TColor::GetFreeColorIndex();
    TColor *cAzureMy = new TColor(kAzureMy, 159. / 255, 191. / 255, 223. / 255, "kAzureMy", 1.0);

    double massD = TDatabasePDG::Instance()->GetParticle(411)->Mass() * 1000;
    double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
    double redMass = massD * massPi / (massD + massPi);

    // now perform global fit
    auto fSCSim = new TF1("fSCSim", GeneralCoulombLednickyTwoRadii, 0., maxRange, nParsSC);
    auto fOCSim = new TF1("fOCSim", GeneralCoulombLednickySecondTwoRadii, 0., maxRange, nParsOC);

    ROOT::Math::WrappedMultiTF1 wfSC(*fSCSim, 1);
    ROOT::Math::WrappedMultiTF1 wfOC(*fOCSim, 1);

    ROOT::Fit::DataOptions opt;
    ROOT::Fit::DataRange rangeSC;
    rangeSC.SetRange(0., maxRange);
    ROOT::Fit::BinData dataSC(opt, rangeSC);
    ROOT::Fit::FillData(dataSC, gSC);

    ROOT::Fit::DataRange rangeOC;
    rangeOC.SetRange(0., maxRange);
    ROOT::Fit::BinData dataOC(opt, rangeOC);
    ROOT::Fit::FillData(dataOC, gOC);

    ROOT::Fit::Chi2Function chi2_SC(dataSC, wfSC);
    ROOT::Fit::Chi2Function chi2_OC(dataOC, wfOC);
    GlobalChi2 globalChi2(chi2_SC, chi2_OC);

    const int nPars = 11;
    double par0[nPars] = {radiusFirst, radiusSecond, weightFirst, -0.08, 0., -0.08, 0., 0., 1., -1., redMass};

    ROOT::Fit::Fitter fitter;
    fitter.Config().SetParamsSettings(nPars, par0);
    fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).Fix();
    fitter.Config().ParSettings(2).Fix();
    fitter.Config().ParSettings(3).SetLowerLimit(-1);
    fitter.Config().ParSettings(3).SetUpperLimit(+1);
    fitter.Config().ParSettings(4).Fix();
    fitter.Config().ParSettings(5).SetLowerLimit(-1);
    fitter.Config().ParSettings(5).SetUpperLimit(+1);
    fitter.Config().ParSettings(6).Fix();
    fitter.Config().ParSettings(7).Fix();
    fitter.Config().ParSettings(8).Fix();
    fitter.Config().ParSettings(9).Fix();
    fitter.Config().ParSettings(10).Fix();

    fitter.Config().MinimizerOptions().SetPrintLevel(0);
    fitter.Config().SetMinimizer("Minuit2", "Migrad");

    // // fit FCN function directly
    // // (specify optionally data size and flag to indicate that is a chi2 fit)
    // // repeat with only stat unc
    fitter.FitFCN(nPars, globalChi2, 0, dataSC.Size() + dataOC.Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();

    return result;

    // const double *params = result.GetParams();
    // const double *parErrors = result.GetErrors();

    // auto lineAtOne = new TLine(0., 1., 300., 1.);
    // lineAtOne->SetLineWidth(1);
    // lineAtOne->SetLineColor(kGray + 1);
    // lineAtOne->SetLineStyle(2);

    // auto lat = new TLatex();
    // lat->SetTextColor(kBlack);
    // lat->SetTextFont(42);
    // lat->SetTextSize(0.045);
    // lat->SetNDC();

    // auto cSimFit = new TCanvas("cSimFit", "", 1000, 500);
    // cSimFit->Divide(2, 1);
    // cSimFit->cd(1)->DrawFrame(0., 0.5, 300., 1.5,
    //                           ";#it{k}* (MeV/#it{c});#it{C}(#it{k}*)");
    // lineAtOne->Draw();
    // fSCSim->SetFitResult(result, iparSC);
    // fSCSim->SetRange(rangeSC().first, rangeSC().second);
    // fSCSim->SetLineColor(kAzure + 4);
    // fSCSimStat->SetFitResult(resultStat, iparSC);
    // fSCSimStat->SetRange(rangeSC().first, rangeSC().second);
    // fSCSimStat->SetLineColor(kAzure + 4);

    // MyFitResult myResultSCStat(fitter.Result(), fSCSim, iparSC);
    // TH1F *hCISCStat = new TH1F("hCISCStat", "", 1000, 0., maxRange);
    // hCISCStat->SetFillColor(kAzureMy);
    // hCISCStat->SetLineColor(kAzureMy);
    // hCISCStat->SetFillColor(kAzureMy);
    // for (auto i{1}; i <= 1000; i++) {
    //   hCISCStat->SetBinContent(i, fSCSim->Eval(hCISCStat->GetBinCenter(i)));
    // }
    // ROOT::Fit::BinData dataSC4CIStat(opt, rangeSC);
    // ROOT::Fit::FillData(dataSC4CIStat, hCISCStat);
    // double ciSCStat[1000];
    // myResultSCStat.GetConfidenceIntervals(dataSC4CIStat, ciSCStat, 0.68);
    // for (auto i{1}; i <= 1000; i++) {
    //   hCISCStat->SetBinError(i, ciSCStat[i - 1]);
    // }

    // auto leg = new TLegend(0.5, 0.75, 0.8, 0.9);
    // leg->SetTextSize(0.045);
    // leg->SetBorderSize(0);
    // leg->SetFillStyle(0);
    // leg->AddEntry(hCISCStat, "stat", "f");

    // hCISC->Draw("e2same");
    // hCISCStat->Draw("e2same");
    // fSCSim->Draw("same");
    // gSC->Draw("p");
    // lat->DrawLatex(
    //     0.2, 0.7,
    //     Form("a_{I=3/2} = %0.2f #pm %0.2f fm", params[3], parErrors[3]));
    // lat->DrawLatex(
    //     0.2, 0.8,
    //     Form("a_{I=1/2} = %0.2f #pm %0.2f fm", params[5], parErrors[5]));

    // cSimFit->cd(2)->DrawFrame(0., 0.5, 300., 1.5,
    //                           ";#it{k}* (MeV/#it{c});#it{C}(#it{k}*)");
    // lineAtOne->Draw();
    // fOCSim->SetFitResult(result, iparOC);
    // fOCSim->SetRange(rangeOC().first, rangeOC().second);
    // fOCSim->SetLineColor(kAzure + 4);
    // fOCSimStat->SetFitResult(resultStat, iparOC);
    // fOCSimStat->SetRange(rangeOC().first, rangeOC().second);
    // fOCSimStat->SetLineColor(kAzure + 4);

    // MyFitResult myResultOCStat(fitter.Result(), fOCSim, iparOC);
    // TH1F *hCIOCStat = new TH1F("hCIOCStat", "", 1000, 0., maxRange);
    // hCIOCStat->SetFillColor(kAzureMy);
    // hCIOCStat->SetLineColor(kAzureMy);
    // hCIOCStat->SetFillColor(kAzureMy);
    // for (auto i{1}; i <= 1000; i++) {
    //   hCIOCStat->SetBinContent(i, fOCSim->Eval(hCIOCStat->GetBinCenter(i)));
    // }
    // ROOT::Fit::BinData dataOC4CIStat(opt, rangeOC);
    // ROOT::Fit::FillData(dataOC4CIStat, hCIOC);
    // double ciOCStat[1000];
    // myResultOCStat.GetConfidenceIntervals(dataOC4CIStat, ciOCStat, 0.68);
    // for (auto i{1}; i <= 1000; i++) {
    //   hCIOCStat->SetBinError(i, ciOCStat[i - 1]);
    // }

    // hCIOCStat->Draw("e2same");
    // fOCSim->Draw("same");
    // gOC->Draw("p");
    // leg->Draw();

    // auto lineVert = new TLine(0., -0.29, 0., 0.29);
    // lineVert->SetLineWidth(1);
    // lineVert->SetLineColor(kGray + 1);
    // lineVert->SetLineStyle(2);

    // auto lineHor = new TLine(-0.29, 0., 0.29, 0.);
    // lineHor->SetLineWidth(1);
    // lineHor->SetLineColor(kGray + 1);
    // lineHor->SetLineStyle(2);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // auto cContour = new TCanvas("cContour", "", 500, 500);

    // TFitResult fitResForContour(result);
    // TGraph *gContour1Sigma = new TGraph(200);
    // gContour1Sigma->SetName("gContour1Sigma");
    // fitResForContour.Contour(3, 5, gContour1Sigma, 0.34);
    // TGraph *gContour2Sigma = new TGraph(200);
    // gContour2Sigma->SetName("gContour2Sigma");
    // fitResForContour.Contour(3, 5, gContour2Sigma, 0.86);
    // gContour1Sigma->SetTitle("");
    // gContour2Sigma->SetTitle("");
    // gContour1Sigma->SetFillColor(kAzureMy);
    // gContour2Sigma->SetFillColor(kAzure + 4);
    // gContour1Sigma->SetLineColor(kAzureMy);
    // gContour2Sigma->SetLineColor(kAzure + 4);

    // TFitResult fitResForContourStat(resultStat);
    // TGraph *gContour1SigmaStat = new TGraph(200);
    // gContour1SigmaStat->SetName("gContour1SigmaStat");
    // fitResForContourStat.Contour(3, 5, gContour1SigmaStat, 0.34);
    // TGraph *gContour2SigmaStat = new TGraph(200);
    // gContour2SigmaStat->SetName("gContour2SigmaStat");
    // fitResForContourStat.Contour(3, 5, gContour2SigmaStat, 0.86);
    // gContour1SigmaStat->SetTitle("");
    // gContour2SigmaStat->SetTitle("");
    // gContour1SigmaStat->SetFillColor(kAzureMy);
    // gContour2SigmaStat->SetFillColor(kAzure + 4);
    // gContour1SigmaStat->SetLineColor(kAzureMy);
    // gContour2SigmaStat->SetLineColor(kAzure + 4);

    // auto legContour = new TLegend(0.7, 0.75, 0.9, 0.9);
    // legContour->SetTextSize(0.045);
    // legContour->SetBorderSize(0);
    // legContour->SetFillStyle(0);
    // legContour->AddEntry(gContour1Sigma, "68% CI", "f");
    // legContour->AddEntry(gContour2Sigma, "95% CI", "f");
    // auto hFrame = cContour->DrawFrame(-0.29, -0.29, 0.29, 0.29,
    //                                   ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
    // lineHor->Draw();
    // lineVert->Draw();
    // hFrame->GetXaxis()->SetDecimals();
    // hFrame->GetYaxis()->SetDecimals();
    // hFrame->GetXaxis()->SetNdivisions(508);
    // hFrame->GetYaxis()->SetNdivisions(508);
    // hFrame->GetYaxis()->SetTitleOffset(1.4);
    // gContour2Sigma->Draw("cf");
    // gContour1Sigma->Draw("cf");
    // legContour->Draw();

    // if (saveResult != kNot) {
    //   std::string outNameTot = "";
    //   std::string outNamePDF = "";
    //   std::string option = "";
    //   if (saveResult == kRecreate) {
    //     outNameTot =
    //         Form("%s_radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s.root",
    //              outName.data(), radiusFirst, radiusSecond, weightFirst,
    //              suffix.data());
    //     outNamePDF =
    //         Form("%s_radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s.pdf",
    //              outName.data(), radiusFirst, radiusSecond, weightFirst,
    //              suffix.data());
    //     option = "create";
    //   } else if (saveResult == kUpdate) {
    //     outNameTot = Form("%s.root", outName.data());
    //     outNamePDF =
    //         Form("%s_radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s.pdf",
    //              outName.data(), radiusFirst, radiusSecond, weightFirst,
    //              suffix.data());
    //     option = "update";
    //     if (deleteFile && !gSystem->Exec(Form("ls %s", outNameTot.data()))) {
    //       gSystem->Exec(Form("rm %s", outNameTot.data()));
    //     }
    //   }

    //   TFile outFile(outNameTot.data(), option.data());
    //   TDirectoryFile *dir = nullptr;
    //   if (saveResult == kUpdate) {
    //     dir = new TDirectoryFile(
    //         Form("radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s",
    //              radiusFirst, radiusSecond, weightFirst, suffix.data()),
    //         Form("radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s",
    //              radiusFirst, radiusSecond, weightFirst, suffix.data()));
    //     outFile.cd();
    //     dir->Write();
    //     dir->cd();
    //   }
    //   cContour->Write();
    //   cSimFit->Write();
    //   gContour2Sigma->Write();
    //   gContour1Sigma->Write();
    //   gContour2SigmaStat->Write();
    //   gContour1SigmaStat->Write();
    //   gSC->Write();
    //   gOC->Write();
    //   gSCStat->Write();
    //   gOCStat->Write();
    //   hCISC->Write();
    //   hCIOC->Write();
    //   hCISCStat->Write();
    //   hCIOCStat->Write();
    //   fSCSim->Write();
    //   fOCSim->Write();
    //   fSCSimStat->Write();
    //   fOCSimStat->Write();
    //   outFile.Close();

    //   if (savePDF) {
    //     cSimFit->SaveAs(outNamePDF.data());
    //   }
    // }

    // delete hCISC;
    // delete hCIOC;
    // delete hCISCStat;
    // delete hCIOCStat;
    // delete cContour;
    // delete cSimFit;

    // return std::map<std::string, double>{
    //   {"par0", params[3]},
    //   {"par1", params[5]},
    //   {"unc0", parErrors[3]},
    //   {"unc1", parErrors[5]},
    //   {"chi2", result.Chi2()},
    //   {"ndf", result.Ndf()}
    // };
}

// void mainSingleFit(TGraphErrors *gSC, TGraphErrors *gOC) {
//   std::array<double, 3> radiusFirst = {0.97, 0.97 + 0.09, 0.97 - 0.08};
//   std::array<double, 3> radiusSecond = {2.52, 2.52 + 0.36, 2.52 - 0.20};
//   std::array<double, 3> weightsFirst = {0.66, 0.66 + 0.03, 0.66 - 0.02};
//   bool deleteFile = true;
//   int counter{0};

//   auto inFileSCStat = TFile::Open("Result_PipDp.root");
//   auto gSCStat = (TGraphAsymmErrors *)inFileSCStat->Get("CF_corr_N_truth");
//   gSCStat->SetName("gSC_corr_N_truth_stat");
//   for (int iPt{0}; iPt < gSCStat->GetN(); ++iPt) {
//     double kStar, cf;
//     gSC->GetPoint(iPt, kStar, cf);
//     gSCStat->SetPoint(iPt, kStar, cf);
//   }
//   auto inFileOCStat = TFile::Open("Result_PipDm.root");
//   auto gOCStat = (TGraphAsymmErrors *)inFileOCStat->Get("CF_corr_N_truth");
//   gOCStat->SetName("gOC_corr_N_truth_stat");
//   for (int iPt{0}; iPt < gOCStat->GetN(); ++iPt) {
//     double kStar, cf;
//     gOC->GetPoint(iPt, kStar, cf);
//     gOCStat->SetPoint(iPt, kStar, cf);
//   }

//   auto hScatterLenFirst =
//       new TH1F("hScatterLenFirst", ";;scattering length (fm)", 28, -0.5, 27.5);
//   auto hScatterLenSecond =
//       new TH1F("hScatterLenSecond", ";;scattering length (fm)", 28, -0.5, 27.5);
//   auto hScatterLenStatFirst = new TH1F(
//       "hScatterStatLenFirst", ";;scattering length (fm)", 28, -0.5, 27.5);
//   auto hScatterLenStatSecond = new TH1F(
//       "hScatterStatLenSecond", ";;scattering length (fm)", 28, -0.5, 27.5);
//   auto hChi2 = new TH1F("hChi2", ";;#chi^{2}/ndf", 28, -0.5, 27.5);
//   hScatterLenFirst->SetLineColor(kAzure + 4);
//   hScatterLenFirst->SetLineWidth(2);
//   hScatterLenFirst->SetMarkerColor(kAzure + 4);
//   hScatterLenFirst->SetMarkerStyle(kFullCircle);
//   hScatterLenSecond->SetLineColor(kRed + 1);
//   hScatterLenSecond->SetLineWidth(2);
//   hScatterLenSecond->SetMarkerColor(kRed + 1);
//   hScatterLenSecond->SetMarkerStyle(kFullCircle);
//   hScatterLenStatFirst->SetLineColor(kAzure + 4);
//   hScatterLenStatFirst->SetLineWidth(2);
//   hScatterLenStatFirst->SetMarkerColor(kAzure + 4);
//   hScatterLenStatFirst->SetMarkerStyle(kFullCircle);
//   hScatterLenStatSecond->SetLineColor(kRed + 1);
//   hScatterLenStatSecond->SetLineWidth(2);
//   hScatterLenStatSecond->SetMarkerColor(kRed + 1);
//   hScatterLenStatSecond->SetMarkerStyle(kFullCircle);
//   hChi2->SetLineColor(kBlack);
//   hChi2->SetLineWidth(2);
//   hChi2->SetMarkerColor(kBlack);
//   hChi2->SetMarkerStyle(kFullCircle);

//   std::string suffix = "";

//   gROOT->SetBatch(true);
//   for (auto &rad1 : radiusFirst) {
//     for (auto &rad2 : radiusSecond) {
//       for (auto &wei : weightsFirst) {
//         if (counter > 0)
//           deleteFile = false;
//         std::map<std::string, double> res = fitSimultaneousLL(
//             gSCStat, gOCStat, gSC, gOC, rad1, rad2, wei, 200., kUpdate,
//             deleteFile, Form("LLSimFit_Dpi%s", suffix.data()), "", true, true);
//         counter++;
//         hScatterLenFirst->SetBinContent(counter, res["par0"]);
//         hScatterLenFirst->SetBinError(counter, res["unc0"]);
//         hScatterLenSecond->SetBinContent(counter, res["par1"]);
//         hScatterLenSecond->SetBinError(counter, res["unc1"]);
//         hScatterLenStatFirst->SetBinContent(counter, res["par0stat"]);
//         hScatterLenStatFirst->SetBinError(counter, res["unc0stat"]);
//         hScatterLenStatSecond->SetBinContent(counter, res["par1stat"]);
//         hScatterLenStatSecond->SetBinError(counter, res["unc1stat"]);
//         hChi2->SetBinContent(counter, res["chi2"] / res["ndf"]);
//         hChi2->SetBinError(counter, 1.e-20);
//         hScatterLenFirst->GetXaxis()->SetBinLabel(
//             counter,
//             Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
//                  rad1, rad2, wei));
//         hScatterLenSecond->GetXaxis()->SetBinLabel(
//             counter,
//             Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
//                  rad1, rad2, wei));
//         hScatterLenStatFirst->GetXaxis()->SetBinLabel(
//             counter,
//             Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
//                  rad1, rad2, wei));
//         hScatterLenStatSecond->GetXaxis()->SetBinLabel(
//             counter,
//             Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
//                  rad1, rad2, wei));
//         hChi2->GetXaxis()->SetBinLabel(
//             counter,
//             Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
//                  rad1, rad2, wei));
//       }
//     }
//   }
//   double rad1 = 1.15;
//   double rad2 = 2.69;
//   double wei = 0.68;
//   std::map<std::string, double> res = fitSimultaneousLL(
//       gSCStat, gOCStat, gSC, gOC, rad1, rad2, wei, 200., kUpdate, deleteFile,
//       Form("LLSimFit_Dpi%s", suffix.data()), "", true, true);
//   hScatterLenFirst->SetBinContent(28, res["par0"]);
//   hScatterLenFirst->SetBinError(28, res["unc0"]);
//   hScatterLenSecond->SetBinContent(28, res["par1"]);
//   hScatterLenSecond->SetBinError(28, res["unc1"]);
//   hScatterLenStatFirst->SetBinContent(28, res["par0stat"]);
//   hScatterLenStatFirst->SetBinError(28, res["unc0stat"]);
//   hScatterLenStatSecond->SetBinContent(28, res["par1stat"]);
//   hScatterLenStatSecond->SetBinError(28, res["unc1stat"]);
//   hChi2->SetBinContent(28, res["chi2"] / res["ndf"]);
//   hChi2->SetBinError(28, 1.e-20);
//   hScatterLenFirst->GetXaxis()->SetBinLabel(
//       28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
//                rad2, wei));
//   hScatterLenSecond->GetXaxis()->SetBinLabel(
//       28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
//                rad2, wei));
//   hScatterLenStatFirst->GetXaxis()->SetBinLabel(
//       28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
//                rad2, wei));
//   hScatterLenStatSecond->GetXaxis()->SetBinLabel(
//       28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
//                rad2, wei));
//   hChi2->GetXaxis()->SetBinLabel(
//       28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
//                rad2, wei));
//   gROOT->SetBatch(false);

//   auto cScattLen = new TCanvas("cScattLen", "", 1500, 500);
//   auto leg = new TLegend(0.2, 0.2, 0.4, 0.4);
//   leg->SetTextSize(0.045);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//   leg->AddEntry(hScatterLenFirst, "a_{I=3/2}");
//   leg->AddEntry(hScatterLenSecond, "a_{I=1/2}");
//   cScattLen->SetLeftMargin(0.08);
//   hScatterLenFirst->GetYaxis()->SetRangeUser(-0.3, 0.1);
//   hScatterLenFirst->Draw();
//   hScatterLenSecond->Draw("same");
//   leg->Draw();
//   cScattLen->SaveAs("ScatteringLength_Dpi_fit_allradii.pdf");

//   auto cChi2 = new TCanvas("cChi2", "", 1500, 500);
//   cChi2->SetLeftMargin(0.08);
//   hChi2->Draw();
//   cChi2->SaveAs("ScatteringLength_Dpi_fitchi2_allradii.pdf");

//   auto lineVert = new TLine(0., -0.29, 0., 0.29);
//   lineVert->SetLineWidth(1);
//   lineVert->SetLineColor(kGray + 1);
//   lineVert->SetLineStyle(2);

//   auto lineHor = new TLine(-0.29, 0., 0.29, 0.);
//   lineHor->SetLineWidth(1);
//   lineHor->SetLineColor(kGray + 1);
//   lineHor->SetLineStyle(2);

//   std::vector<TGraph *> gContour1Sigma{}, gContour2Sigma{};
//   auto inFile = TFile::Open(Form("LLSimFit_Dpi%s.root", suffix.data()));
//   TDirectoryFile *dir = nullptr;
//   TKey *key;
//   TIter next(inFile->GetListOfKeys());
//   while ((key = (TKey *)next())) {
//     dir = dynamic_cast<TDirectoryFile *>(inFile->Get(key->GetName()));
//     gContour1Sigma.push_back(
//         dynamic_cast<TGraph *>(dir->Get("gContour1Sigma")));
//     gContour2Sigma.push_back(
//         dynamic_cast<TGraph *>(dir->Get("gContour2Sigma")));
//   }

//   auto cContourAll = new TCanvas("cContourAll", "", 900, 900);
//   cContourAll->DrawFrame(-0.29, -0.29, 0.29, 0.29,
//                          ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
//   lineHor->Draw();
//   lineVert->Draw();

//   for (int i{0}; i < 28; ++i) {
//     gContour2Sigma[i]->Draw("cf");
//   }
//   for (int i{0}; i < 28; ++i) {
//     gContour1Sigma[i]->Draw("cf");
//   }

//   auto legContour = new TLegend(0.7, 0.75, 0.9, 0.9);
//   legContour->SetTextSize(0.045);
//   legContour->SetBorderSize(0);
//   legContour->SetFillStyle(0);
//   legContour->AddEntry(gContour1Sigma[0], "68% CI", "f");
//   legContour->AddEntry(gContour2Sigma[0], "95% CI", "f");
//   legContour->Draw();

//   auto lat = new TLatex();
//   lat->SetTextColor(kBlack);
//   lat->SetTextFont(42);
//   lat->SetTextSize(0.07);
//   lat->SetNDC();
//   cContourAll->SaveAs("ContourPlot_Dpi_fit_allradii_convolution.pdf");

//   auto cContour = new TCanvas("cContour", "", 1920, 1080);
//   cContour->Divide(10, 3);
//   int iPad{0};
//   for (auto &rad1 : radiusFirst) {
//     for (auto &rad2 : radiusSecond) {
//       for (auto &wei : weightsFirst) {
//         ++iPad;
//         auto hFrame = cContour->cd(iPad)->DrawFrame(
//             -0.29, -0.29, 0.29, 0.29, ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
//         hFrame->GetXaxis()->SetTitleSize(0.07);
//         hFrame->GetYaxis()->SetTitleSize(0.07);
//         hFrame->GetXaxis()->SetLabelSize(0.07);
//         hFrame->GetYaxis()->SetLabelSize(0.07);
//         hFrame->GetXaxis()->SetTitleOffset(0.9);
//         hFrame->GetYaxis()->SetTitleOffset(0.9);
//         lineHor->Draw();
//         lineVert->Draw();
//         gContour2Sigma[iPad - 1]->Draw("cf");
//         gContour1Sigma[iPad - 1]->Draw("cf");
//         lat->DrawLatex(0.2, 0.9, Form("#it{r}_{1} = %0.2f fm", rad1));
//         lat->DrawLatex(0.2, 0.8, Form("#it{r}_{2} = %0.2f fm", rad2));
//         lat->DrawLatex(0.2, 0.7, Form("#it{w}_{1} = %0.2f fm", wei));
//       }
//     }
//   }
//   auto hFrame = cContour->cd(iPad + 1)->DrawFrame(
//       -0.29, -0.29, 0.29, 0.29, ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
//   hFrame->GetXaxis()->SetTitleSize(0.07);
//   hFrame->GetYaxis()->SetTitleSize(0.07);
//   hFrame->GetXaxis()->SetLabelSize(0.07);
//   hFrame->GetYaxis()->SetLabelSize(0.07);
//   hFrame->GetXaxis()->SetTitleOffset(0.9);
//   hFrame->GetYaxis()->SetTitleOffset(0.9);
//   lineHor->Draw();
//   lineVert->Draw();
//   gContour2Sigma[iPad]->Draw("cf");
//   gContour1Sigma[iPad]->Draw("cf");
//   lat->DrawLatex(0.2, 0.9, Form("#it{r}_{1} = %0.2f fm", rad1));
//   lat->DrawLatex(0.2, 0.8, Form("#it{r}_{2} = %0.2f fm", rad2));
//   lat->DrawLatex(0.2, 0.7, Form("#it{w}_{1} = %0.2f fm", wei));
//   cContour->SaveAs("ContourPlot_Dpi_fit_allradii.pdf");

//   // add some info to the output files
//   auto fileToUpdate =
//       TFile::Open(Form("LLSimFit_Dpi%s.root", suffix.data()), "update");
//   fileToUpdate->cd();
//   TDirectoryFile *dirScatt =
//       new TDirectoryFile("scattering_length", "scattering_length");
//   dirScatt->Write();
//   dirScatt->cd();
//   hScatterLenFirst->Write();
//   hScatterLenSecond->Write();
//   hScatterLenStatFirst->Write();
//   hScatterLenStatSecond->Write();
//   hChi2->Write();
//   cScattLen->Write();
//   cChi2->Write();
//   dirScatt->Close();

//   fileToUpdate->cd();
//   TDirectoryFile *dirContour = new TDirectoryFile("contours", "contours");
//   dirContour->Write();
//   dirContour->cd();
//   cContourAll->Write();
//   cContour->Write();
//   dirContour->Close();

//   fileToUpdate->cd();
//   TDirectoryFile *dirCoulomb = new TDirectoryFile("Coulomb", "Coulomb");
//   dirCoulomb->Write();
//   dirCoulomb->cd();
//   gCoulombSCmin->Write();
//   gCoulombSCmax->Write();
//   gCoulombOCmin->Write();
//   gCoulombOCmax->Write();
//   dirCoulomb->Close();
//   fileToUpdate->Close();
// }

void fitSimultaneousLLDpi(const char *inFileName, const char *oFileName, int nIter) {
    TFile *inFile = new TFile(inFileName);
    TFile *oFile = new TFile(oFileName, "create");

    for (auto uncType : {"tot", "stat"}) {
        oFile->mkdir(Form("%s/iters", uncType));
        oFile->cd(Form("%s/iters", uncType));
        TNtuple *tResults = new TNtuple("tResults", "", "a0sin:a0tri:chi2ndf:r1:r2:w1");

        TCanvas *cCombFit = new TCanvas("cCombFit", "", 1200, 600);
        cCombFit->Divide(3, 1);
        cCombFit->SaveAs(
            Form("/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs%d%s.pdf[", nIter, uncType));

        std::vector<double> weights1 = {0.66, 0.69, 0.64};
        std::vector<double> radii1 = {0.97, 1.06, 0.89};
        std::vector<double> radii2 = {2.52, 2.88, 2.32};
        std::vector<std::vector<double>> fitRanges = {{10, 450}, {10, 400}, {10, 500}};

        for (int iIter = 0; iIter < nIter; iIter++) {
            TGraphAsymmErrors *gSC =
                reinterpret_cast<TGraphAsymmErrors *>(inFile->Get(Form("sc/%s/gCFGen%d", uncType, iIter)));
            gSC->SetName(Form("gSC%s%d", uncType, iIter));
            gSC->Write();
            TGraphAsymmErrors *gOC =
                reinterpret_cast<TGraphAsymmErrors *>(inFile->Get(Form("oc/%s/gCFGen%d", uncType, iIter)));
            gOC->SetName(Form("gOC%s%d", uncType, iIter));
            gOC->Write();

            int iRadius = gRandom->Integer(3);
            int iFitRange = gRandom->Integer(3);

            ROOT::Fit::FitResult result = fitSimultaneousLL(gSC, gOC, radii1[iRadius], radii2[iRadius],
                                                            weights1[iRadius], fitRanges[iFitRange][1]);
            const double *params = result.GetParams();
            const double *parErrors = result.GetErrors();

            double a0singlet = params[3];
            double a0triplet = params[5];
            double a0singletUnc = parErrors[3];
            double a0tripletUnc = parErrors[5];

            double chi2ndf = result.Chi2() / result.Ndf();
            if (std::abs(a0singlet) < 1 && std::abs(a0triplet) < 1 && chi2ndf < 10)
                tResults->Fill(a0singlet, a0triplet, chi2ndf, radii1[iRadius], radii2[iRadius], weights1[iRadius]);

            // Draw
            TVirtualPad *pad1 = cCombFit->cd(1);
            pad1->DrawFrame(0, 0.8, 300, 1.4);
            auto fSCSim = new TF1("fSCSim", GeneralCoulombLednickyTwoRadii, 0., 300, nParsSC);

            fSCSim->FixParameter(0, params[0]);
            fSCSim->FixParameter(1, params[1]);
            fSCSim->FixParameter(2, params[2]);
            fSCSim->FixParameter(3, params[3]);
            fSCSim->FixParameter(4, params[4]);
            fSCSim->FixParameter(5, params[7]);
            fSCSim->FixParameter(6, params[8]);
            fSCSim->FixParameter(7, params[10]);
            fSCSim->Draw("same");
            gSC->Draw("same");

            TVirtualPad *pad2 = cCombFit->cd(2);
            pad2->DrawFrame(0, 0.8, 300, 1.4);
            auto fOCSim = new TF1("fSCSim", GeneralCoulombLednickySecondTwoRadii, 0., 300, nParsOC);
            fOCSim->FixParameter(0, params[0]);
            fOCSim->FixParameter(1, params[1]);
            fOCSim->FixParameter(2, params[2]);
            fOCSim->FixParameter(3, params[3]);
            fOCSim->FixParameter(4, params[4]);
            fOCSim->FixParameter(5, params[5]);
            fOCSim->FixParameter(6, params[6]);
            fOCSim->FixParameter(7, params[7]);
            fOCSim->FixParameter(8, params[9]);
            fOCSim->FixParameter(9, params[10]);
            fOCSim->Draw("same");
            gOC->Draw("same");

            TVirtualPad *pad3 = cCombFit->cd(3);
            pad3->DrawFrame(-0.3, -0.3, 0.3, 0.3);
            TGraphErrors *gScattPar = new TGraphErrors(1);
            gScattPar->SetName("gScattPar");
            gScattPar->SetPoint(0, a0singlet, a0triplet);
            gScattPar->SetPointError(0, a0singletUnc, a0tripletUnc);
            gScattPar->SetMarkerStyle(20);
            gScattPar->SetMarkerSize(1);
            gScattPar->Draw("same");
            cCombFit->SaveAs(
                Form("/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs%d%s.pdf", nIter, uncType));
            // cCombFit->Write();
            delete gSC;
            delete gOC;
        }
        cCombFit->SaveAs(
            Form("/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs%d%s.pdf]", nIter, uncType));
        oFile->cd(uncType);
        tResults->Write();
        delete tResults;
        delete cCombFit;
    }
    cCombFit->SaveAs(Form("/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs%d%s.pdf]", nIter, uncType));
    oFile->cd(uncType);
    tResults->Write();
    delete tResults;
    delete cCombFit;
}
