#include "functions.h"

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

enum { kRecreate = 0, kUpdate, kNot };

// definition of shared parameter
// same charge function
const int nParsSC = 8;
int iparSC[nParsSC] = {
    0, // first source radius
    1, // second source radius
    2, // normalisation of two contributions
    3, // scattering length singlet
    4, // effective range singlet
    7, // QS
    8, // charge product
    10 // reduced mass
};

// opposite charge function
const int nParsOC = 10;
int iparOC[nParsOC] = {
    0, // first source radius
    1, // second source radius
    2, // normalisation of two contributions
    3, // scattering length singlet
    4, // effective range singlet
    5, // scattering length triplet
    6, // effective range triplet
    7, // QS
    9, // charge product
    10 // reduced mass
};

// Create the GlobalCHi2 structure
struct GlobalChi2 {
  GlobalChi2(ROOT::Math::IMultiGenFunction &f1,
             ROOT::Math::IMultiGenFunction &f2)
      : fChi2_1(&f1), fChi2_2(&f2) {}

  // parameter vector is first background (in common 1 and 2)
  // and then is signal (only in 2)
  double operator()(const double *par) const {

    double p1[nParsSC];
    for (int i = 0; i < nParsSC; ++i) {
      p1[i] = par[iparSC[i]];
    }
    double p2[nParsOC];
    for (int i = 0; i < nParsOC; ++i) {
      p2[i] = par[iparOC[i]];
    }

    return (*fChi2_1)(p1) + (*fChi2_2)(p2);
  }

  const ROOT::Math::IMultiGenFunction *fChi2_1;
  const ROOT::Math::IMultiGenFunction *fChi2_2;
};

std::array<double, 2> ComputeChi2(TGraphAsymmErrors *graphFirst,
                                  TGraphAsymmErrors *graphSecond,
                                  TF1 *funcFirst, TF1 *funcSecond,
                                  double kStarMax = 200., int npars = 2) {
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

std::map<std::string, double>
fitSimultaneousLL(TGraphAsymmErrors *gSCStat, TGraphAsymmErrors *gOCStat,
                  TGraphAsymmErrors *gSC, TGraphAsymmErrors *gOC,
                  double radiusFirst = 0.97, double radiusSecond = 2.52,
                  double weightFirst = 0.66, double maxRange = 200.,
                  int saveResult = kRecreate, bool deleteFile = false,
                  std::string outName = "LLSimFit_Dpi", std::string suffix = "",
                  bool savePDF = false, bool fixRadii = true) {
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
  TColor *cAzureMy =
      new TColor(kAzureMy, 159. / 255, 191. / 255, 223. / 255, "kAzureMy", 1.0);

  double massD = TDatabasePDG::Instance()->GetParticle(411)->Mass() * 1000;
  double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
  double redMass = massD * massPi / (massD + massPi);

  // now perform global fit
  auto fSCSim =
      new TF1("fSCSim", GeneralCoulombLednickyTwoRadii, 0., maxRange, nParsSC);
  auto fOCSim = new TF1("fOCSim", GeneralCoulombLednickySecondTwoRadii, 0.,
                        maxRange, nParsOC);
  auto fSCSimStat = new TF1("fSCSimStat", GeneralCoulombLednickyTwoRadii, 0.,
                            maxRange, nParsSC);
  auto fOCSimStat = new TF1("fOCSimStat", GeneralCoulombLednickySecondTwoRadii,
                            0., maxRange, nParsOC);
  ROOT::Math::WrappedMultiTF1 wfSC(*fSCSim, 1);
  ROOT::Math::WrappedMultiTF1 wfOC(*fOCSim, 1);
  ROOT::Math::WrappedMultiTF1 wfSCStat(*fSCSimStat, 1);
  ROOT::Math::WrappedMultiTF1 wfOCStat(*fOCSimStat, 1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeSC;
  rangeSC.SetRange(0., maxRange);
  ROOT::Fit::BinData dataSC(opt, rangeSC);
  ROOT::Fit::FillData(dataSC, gSC);
  ROOT::Fit::BinData dataSCStat(opt, rangeSC);
  ROOT::Fit::FillData(dataSCStat, gSCStat);

  ROOT::Fit::DataRange rangeOC;
  rangeOC.SetRange(0., maxRange);
  ROOT::Fit::BinData dataOC(opt, rangeOC);
  ROOT::Fit::FillData(dataOC, gOC);
  ROOT::Fit::BinData dataOCStat(opt, rangeOC);
  ROOT::Fit::FillData(dataOCStat, gOCStat);

  ROOT::Fit::Chi2Function chi2_SC(dataSC, wfSC);
  ROOT::Fit::Chi2Function chi2_OC(dataOC, wfOC);
  GlobalChi2 globalChi2(chi2_SC, chi2_OC);

  ROOT::Fit::Chi2Function chi2_SCStat(dataSCStat, wfSC);
  ROOT::Fit::Chi2Function chi2_OCStat(dataOCStat, wfOC);
  GlobalChi2 globalChi2Stat(chi2_SCStat, chi2_OCStat);

  const int nPars = 11;
  double par0[nPars] = {radiusFirst, radiusSecond, weightFirst, -0.08,
                        0.,          -0.08,        0.,          0.,
                        1.,          -1.,          redMass};

  ROOT::Fit::Fitter fitter;
  fitter.Config().SetParamsSettings(nPars, par0);
  if (!fixRadii) {
    fitter.Config().ParSettings(0).SetLimits(radiusFirst - 0.08,
                                             radiusFirst + 0.09);
    fitter.Config().ParSettings(1).SetLimits(radiusSecond - 0.20,
                                             radiusSecond + 0.36);
    fitter.Config().ParSettings(2).SetLimits(weightFirst - 0.02,
                                             weightFirst + 0.03);
  } else {
    fitter.Config().ParSettings(0).Fix();
    fitter.Config().ParSettings(1).Fix();
    fitter.Config().ParSettings(2).Fix();
  }
  fitter.Config().ParSettings(4).Fix();
  fitter.Config().ParSettings(6).Fix();
  fitter.Config().ParSettings(7).Fix();
  fitter.Config().ParSettings(8).Fix();
  fitter.Config().ParSettings(9).Fix();
  fitter.Config().ParSettings(10).Fix();

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit2", "Migrad");

  ROOT::Fit::Fitter fitterStat;
  fitterStat.Config().SetParamsSettings(nPars, par0);
  if (!fixRadii) {
    fitterStat.Config().ParSettings(0).SetLimits(radiusFirst - 0.08,
                                                 radiusFirst + 0.09);
    fitterStat.Config().ParSettings(1).SetLimits(radiusSecond - 0.20,
                                                 radiusSecond + 0.32);
    fitterStat.Config().ParSettings(2).SetLimits(weightFirst - 0.02,
                                                 weightFirst + 0.03);
  } else {
    fitterStat.Config().ParSettings(0).Fix();
    fitterStat.Config().ParSettings(1).Fix();
    fitterStat.Config().ParSettings(2).Fix();
  }
  fitterStat.Config().ParSettings(4).Fix();
  fitterStat.Config().ParSettings(6).Fix();
  fitterStat.Config().ParSettings(7).Fix();
  fitterStat.Config().ParSettings(8).Fix();
  fitterStat.Config().ParSettings(9).Fix();
  fitterStat.Config().ParSettings(10).Fix();

  fitterStat.Config().MinimizerOptions().SetPrintLevel(0);
  fitterStat.Config().SetMinimizer("Minuit2", "Migrad");

  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(nPars, globalChi2, 0, dataSC.Size() + dataOC.Size(), true);
  ROOT::Fit::FitResult result = fitter.Result();
  // repeat with only stat unc
  fitterStat.FitFCN(nPars, globalChi2Stat, 0,
                    dataSCStat.Size() + dataOCStat.Size(), true);
  ROOT::Fit::FitResult resultStat = fitterStat.Result();

  const double *params = result.GetParams();
  const double *parErrors = result.GetErrors();

  const double *paramsStat = resultStat.GetParams();
  const double *parErrorsStat = resultStat.GetErrors();

  auto lineAtOne = new TLine(0., 1., 300., 1.);
  lineAtOne->SetLineWidth(1);
  lineAtOne->SetLineColor(kGray + 1);
  lineAtOne->SetLineStyle(2);

  auto lat = new TLatex();
  lat->SetTextColor(kBlack);
  lat->SetTextFont(42);
  lat->SetTextSize(0.045);
  lat->SetNDC();

  auto cSimFit = new TCanvas("cSimFit", "", 1000, 500);
  cSimFit->Divide(2, 1);
  cSimFit->cd(1)->DrawFrame(0., 0.5, 300., 1.5,
                            ";#it{k}* (MeV/#it{c});#it{C}(#it{k}*)");
  lineAtOne->Draw();
  fSCSim->SetFitResult(result, iparSC);
  fSCSim->SetRange(rangeSC().first, rangeSC().second);
  fSCSim->SetLineColor(kAzure + 4);
  fSCSimStat->SetFitResult(resultStat, iparSC);
  fSCSimStat->SetRange(rangeSC().first, rangeSC().second);
  fSCSimStat->SetLineColor(kAzure + 4);

  MyFitResult myResultSC(fitter.Result(), fSCSim, iparSC);
  TH1F *hCISC = new TH1F("hCISC", "", 1000, 0., maxRange);
  hCISC->SetFillColor(kAzure + 4);
  hCISC->SetLineColor(kAzure + 4);
  hCISC->SetFillColor(kAzure + 4);
  for (auto i{1}; i <= 1000; i++) {
    hCISC->SetBinContent(i, fSCSim->Eval(hCISC->GetBinCenter(i)));
  }
  ROOT::Fit::BinData dataSC4CI(opt, rangeSC);
  ROOT::Fit::FillData(dataSC4CI, hCISC);
  double ciSC[1000];
  myResultSC.GetConfidenceIntervals(dataSC4CI, ciSC, 0.68);
  for (auto i{1}; i <= 1000; i++) {
    hCISC->SetBinError(i, ciSC[i - 1]);
  }

  MyFitResult myResultSCStat(fitterStat.Result(), fSCSim, iparSC);
  TH1F *hCISCStat = new TH1F("hCISCStat", "", 1000, 0., maxRange);
  hCISCStat->SetFillColor(kAzureMy);
  hCISCStat->SetLineColor(kAzureMy);
  hCISCStat->SetFillColor(kAzureMy);
  for (auto i{1}; i <= 1000; i++) {
    hCISCStat->SetBinContent(i, fSCSim->Eval(hCISCStat->GetBinCenter(i)));
  }
  ROOT::Fit::BinData dataSC4CIStat(opt, rangeSC);
  ROOT::Fit::FillData(dataSC4CIStat, hCISCStat);
  double ciSCStat[1000];
  myResultSCStat.GetConfidenceIntervals(dataSC4CIStat, ciSCStat, 0.68);
  for (auto i{1}; i <= 1000; i++) {
    hCISCStat->SetBinError(i, ciSCStat[i - 1]);
  }

  auto leg = new TLegend(0.5, 0.75, 0.8, 0.9);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hCISCStat, "stat", "f");
  leg->AddEntry(hCISC, "#sqrt{stat^{2}+syst^{2}}", "f");

  hCISC->Draw("e2same");
  hCISCStat->Draw("e2same");
  fSCSim->Draw("same");
  gSC->Draw("p");
  lat->DrawLatex(
      0.2, 0.7,
      Form("a_{I=3/2} = %0.2f #pm %0.2f fm", params[3], parErrors[3]));
  lat->DrawLatex(
      0.2, 0.8,
      Form("a_{I=1/2} = %0.2f #pm %0.2f fm", params[5], parErrors[5]));

  cSimFit->cd(2)->DrawFrame(0., 0.5, 300., 1.5,
                            ";#it{k}* (MeV/#it{c});#it{C}(#it{k}*)");
  lineAtOne->Draw();
  fOCSim->SetFitResult(result, iparOC);
  fOCSim->SetRange(rangeOC().first, rangeOC().second);
  fOCSim->SetLineColor(kAzure + 4);
  fOCSimStat->SetFitResult(resultStat, iparOC);
  fOCSimStat->SetRange(rangeOC().first, rangeOC().second);
  fOCSimStat->SetLineColor(kAzure + 4);

  MyFitResult myResultOC(fitter.Result(), fOCSim, iparOC);
  TH1F *hCIOC = new TH1F("hCIOC", "", 1000, 0., maxRange);
  hCIOC->SetFillColor(kAzure + 4);
  hCIOC->SetLineColor(kAzure + 4);
  hCIOC->SetFillColor(kAzure + 4);
  for (auto i{1}; i <= 1000; i++) {
    hCIOC->SetBinContent(i, fOCSim->Eval(hCIOC->GetBinCenter(i)));
  }
  ROOT::Fit::BinData dataOC4CI(opt, rangeOC);
  ROOT::Fit::FillData(dataOC4CI, hCIOC);
  double ciOC[1000];
  myResultOC.GetConfidenceIntervals(dataOC4CI, ciOC, 0.68);
  for (auto i{1}; i <= 1000; i++) {
    hCIOC->SetBinError(i, ciOC[i - 1]);
  }

  MyFitResult myResultOCStat(fitterStat.Result(), fOCSim, iparOC);
  TH1F *hCIOCStat = new TH1F("hCIOCStat", "", 1000, 0., maxRange);
  hCIOCStat->SetFillColor(kAzureMy);
  hCIOCStat->SetLineColor(kAzureMy);
  hCIOCStat->SetFillColor(kAzureMy);
  for (auto i{1}; i <= 1000; i++) {
    hCIOCStat->SetBinContent(i, fOCSim->Eval(hCIOCStat->GetBinCenter(i)));
  }
  ROOT::Fit::BinData dataOC4CIStat(opt, rangeOC);
  ROOT::Fit::FillData(dataOC4CIStat, hCIOC);
  double ciOCStat[1000];
  myResultOCStat.GetConfidenceIntervals(dataOC4CIStat, ciOCStat, 0.68);
  for (auto i{1}; i <= 1000; i++) {
    hCIOCStat->SetBinError(i, ciOCStat[i - 1]);
  }

  hCIOC->Draw("e2same");
  hCIOCStat->Draw("e2same");
  fOCSim->Draw("same");
  gOC->Draw("p");
  leg->Draw();

  auto lineVert = new TLine(0., -0.29, 0., 0.29);
  lineVert->SetLineWidth(1);
  lineVert->SetLineColor(kGray + 1);
  lineVert->SetLineStyle(2);

  auto lineHor = new TLine(-0.29, 0., 0.29, 0.);
  lineHor->SetLineWidth(1);
  lineHor->SetLineColor(kGray + 1);
  lineHor->SetLineStyle(2);

  auto cContour = new TCanvas("cContour", "", 500, 500);

  TFitResult fitResForContour(result);
  TGraph *gContour1Sigma = new TGraph(200);
  gContour1Sigma->SetName("gContour1Sigma");
  fitResForContour.Contour(3, 5, gContour1Sigma, 0.34);
  TGraph *gContour2Sigma = new TGraph(200);
  gContour2Sigma->SetName("gContour2Sigma");
  fitResForContour.Contour(3, 5, gContour2Sigma, 0.86);
  gContour1Sigma->SetTitle("");
  gContour2Sigma->SetTitle("");
  gContour1Sigma->SetFillColor(kAzureMy);
  gContour2Sigma->SetFillColor(kAzure + 4);
  gContour1Sigma->SetLineColor(kAzureMy);
  gContour2Sigma->SetLineColor(kAzure + 4);

  TFitResult fitResForContourStat(resultStat);
  TGraph *gContour1SigmaStat = new TGraph(200);
  gContour1SigmaStat->SetName("gContour1SigmaStat");
  fitResForContourStat.Contour(3, 5, gContour1SigmaStat, 0.34);
  TGraph *gContour2SigmaStat = new TGraph(200);
  gContour2SigmaStat->SetName("gContour2SigmaStat");
  fitResForContourStat.Contour(3, 5, gContour2SigmaStat, 0.86);
  gContour1SigmaStat->SetTitle("");
  gContour2SigmaStat->SetTitle("");
  gContour1SigmaStat->SetFillColor(kAzureMy);
  gContour2SigmaStat->SetFillColor(kAzure + 4);
  gContour1SigmaStat->SetLineColor(kAzureMy);
  gContour2SigmaStat->SetLineColor(kAzure + 4);

  auto legContour = new TLegend(0.7, 0.75, 0.9, 0.9);
  legContour->SetTextSize(0.045);
  legContour->SetBorderSize(0);
  legContour->SetFillStyle(0);
  legContour->AddEntry(gContour1Sigma, "68% CI", "f");
  legContour->AddEntry(gContour2Sigma, "95% CI", "f");
  auto hFrame = cContour->DrawFrame(-0.29, -0.29, 0.29, 0.29,
                                    ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
  lineHor->Draw();
  lineVert->Draw();
  hFrame->GetXaxis()->SetDecimals();
  hFrame->GetYaxis()->SetDecimals();
  hFrame->GetXaxis()->SetNdivisions(508);
  hFrame->GetYaxis()->SetNdivisions(508);
  hFrame->GetYaxis()->SetTitleOffset(1.4);
  gContour2Sigma->Draw("cf");
  gContour1Sigma->Draw("cf");
  legContour->Draw();

  if (saveResult != kNot) {
    std::string outNameTot = "";
    std::string outNamePDF = "";
    std::string option = "";
    if (saveResult == kRecreate) {
      outNameTot =
          Form("%s_radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s.root",
               outName.data(), radiusFirst, radiusSecond, weightFirst,
               suffix.data());
      outNamePDF =
          Form("%s_radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s.pdf",
               outName.data(), radiusFirst, radiusSecond, weightFirst,
               suffix.data());
      option = "recreate";
    } else if (saveResult == kUpdate) {
      outNameTot = Form("%s.root", outName.data());
      outNamePDF =
          Form("%s_radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s.pdf",
               outName.data(), radiusFirst, radiusSecond, weightFirst,
               suffix.data());
      option = "update";
      if (deleteFile && !gSystem->Exec(Form("ls %s", outNameTot.data()))) {
        gSystem->Exec(Form("rm %s", outNameTot.data()));
      }
    }

    TFile outFile(outNameTot.data(), option.data());
    TDirectoryFile *dir = nullptr;
    if (saveResult == kUpdate) {
      dir = new TDirectoryFile(
          Form("radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s",
               radiusFirst, radiusSecond, weightFirst, suffix.data()),
          Form("radiusFirst%0.2f_radiusSecond%0.2f_weightFirst%0.2f_%s",
               radiusFirst, radiusSecond, weightFirst, suffix.data()));
      outFile.cd();
      dir->Write();
      dir->cd();
    }
    cContour->Write();
    cSimFit->Write();
    gContour2Sigma->Write();
    gContour1Sigma->Write();
    gContour2SigmaStat->Write();
    gContour1SigmaStat->Write();
    gSC->Write();
    gOC->Write();
    gSCStat->Write();
    gOCStat->Write();
    hCISC->Write();
    hCIOC->Write();
    hCISCStat->Write();
    hCIOCStat->Write();
    fSCSim->Write();
    fOCSim->Write();
    fSCSimStat->Write();
    fOCSimStat->Write();
    outFile.Close();

    if (savePDF) {
      cSimFit->SaveAs(outNamePDF.data());
    }
  }

  delete hCISC;
  delete hCIOC;
  delete hCISCStat;
  delete hCIOCStat;
  delete cContour;
  delete cSimFit;

  return std::map<std::string, double>{{"par0", params[3]},
                                       {"par1", params[5]},
                                       {"par0stat", paramsStat[3]},
                                       {"par1stat", paramsStat[5]},
                                       {"unc0", parErrors[3]},
                                       {"unc1", parErrors[5]},
                                       {"unc0stat", parErrorsStat[3]},
                                       {"unc1stat", parErrorsStat[5]},
                                       {"chi2", result.Chi2()},
                                       {"ndf", result.Ndf()}};
}

void mainSingleFit(bool fix = true) {
  std::array<double, 3> radiusFirst = {0.97, 0.97 + 0.09, 0.97 - 0.08};
  std::array<double, 3> radiusSecond = {2.52, 2.52 + 0.36, 2.52 - 0.20};
  std::array<double, 3> weightsFirst = {0.66, 0.66 + 0.03, 0.66 - 0.02};
  bool deleteFile = true;
  int counter{0};

  auto inFileSC = TFile::Open("Results_PipDp_corrected.root");
  auto gSC = (TGraphAsymmErrors *)inFileSC->Get("CF_corr_N_truth");
  auto gCoulombSCmin = (TGraphAsymmErrors *)inFileSC->Get("D_Coulomb1");
  auto gCoulombSCmax = (TGraphAsymmErrors *)inFileSC->Get("D_Coulomb2");
  gSC->SetName("gSC_corr_N_truth_stat_plus_sys");
  gCoulombSCmin->SetName("gSC_coulomb_min");
  gCoulombSCmax->SetName("gSC_coulomb_max");
  auto inFileOC = TFile::Open("Results_PipDm_corrected.root");
  auto gOC = (TGraphAsymmErrors *)inFileOC->Get("CF_corr_N_truth");
  auto gCoulombOCmin = (TGraphAsymmErrors *)inFileOC->Get("D_Coulomb1");
  auto gCoulombOCmax = (TGraphAsymmErrors *)inFileOC->Get("D_Coulomb2");
  gOC->SetName("gOC_corr_N_truth_stat_plus_sys");
  gCoulombOCmin->SetName("gOC_coulomb_min");
  gCoulombOCmax->SetName("gOC_coulomb_max");

  auto inFileSCStat = TFile::Open("Result_PipDp.root");
  auto gSCStat = (TGraphAsymmErrors *)inFileSCStat->Get("CF_corr_N_truth");
  gSCStat->SetName("gSC_corr_N_truth_stat");
  for (int iPt{0}; iPt < gSCStat->GetN(); ++iPt) {
    double kStar, cf;
    gSC->GetPoint(iPt, kStar, cf);
    gSCStat->SetPoint(iPt, kStar, cf);
  }
  auto inFileOCStat = TFile::Open("Result_PipDm.root");
  auto gOCStat = (TGraphAsymmErrors *)inFileOCStat->Get("CF_corr_N_truth");
  gOCStat->SetName("gOC_corr_N_truth_stat");
  for (int iPt{0}; iPt < gOCStat->GetN(); ++iPt) {
    double kStar, cf;
    gOC->GetPoint(iPt, kStar, cf);
    gOCStat->SetPoint(iPt, kStar, cf);
  }

  auto hScatterLenFirst =
      new TH1F("hScatterLenFirst", ";;scattering length (fm)", 28, -0.5, 27.5);
  auto hScatterLenSecond =
      new TH1F("hScatterLenSecond", ";;scattering length (fm)", 28, -0.5, 27.5);
  auto hScatterLenStatFirst = new TH1F(
      "hScatterStatLenFirst", ";;scattering length (fm)", 28, -0.5, 27.5);
  auto hScatterLenStatSecond = new TH1F(
      "hScatterStatLenSecond", ";;scattering length (fm)", 28, -0.5, 27.5);
  auto hChi2 = new TH1F("hChi2", ";;#chi^{2}/ndf", 28, -0.5, 27.5);
  hScatterLenFirst->SetLineColor(kAzure + 4);
  hScatterLenFirst->SetLineWidth(2);
  hScatterLenFirst->SetMarkerColor(kAzure + 4);
  hScatterLenFirst->SetMarkerStyle(kFullCircle);
  hScatterLenSecond->SetLineColor(kRed + 1);
  hScatterLenSecond->SetLineWidth(2);
  hScatterLenSecond->SetMarkerColor(kRed + 1);
  hScatterLenSecond->SetMarkerStyle(kFullCircle);
  hScatterLenStatFirst->SetLineColor(kAzure + 4);
  hScatterLenStatFirst->SetLineWidth(2);
  hScatterLenStatFirst->SetMarkerColor(kAzure + 4);
  hScatterLenStatFirst->SetMarkerStyle(kFullCircle);
  hScatterLenStatSecond->SetLineColor(kRed + 1);
  hScatterLenStatSecond->SetLineWidth(2);
  hScatterLenStatSecond->SetMarkerColor(kRed + 1);
  hScatterLenStatSecond->SetMarkerStyle(kFullCircle);
  hChi2->SetLineColor(kBlack);
  hChi2->SetLineWidth(2);
  hChi2->SetMarkerColor(kBlack);
  hChi2->SetMarkerStyle(kFullCircle);

  std::string suffix = "";
  if (!fix) {
    suffix = "_radii_limited";
  }
  gROOT->SetBatch(true);
  for (auto &rad1 : radiusFirst) {
    for (auto &rad2 : radiusSecond) {
      for (auto &wei : weightsFirst) {
        if (counter > 0)
          deleteFile = false;
        std::map<std::string, double> res = fitSimultaneousLL(
            gSCStat, gOCStat, gSC, gOC, rad1, rad2, wei, 200., kUpdate,
            deleteFile, Form("LLSimFit_Dpi%s", suffix.data()), "", true, fix);
        counter++;
        hScatterLenFirst->SetBinContent(counter, res["par0"]);
        hScatterLenFirst->SetBinError(counter, res["unc0"]);
        hScatterLenSecond->SetBinContent(counter, res["par1"]);
        hScatterLenSecond->SetBinError(counter, res["unc1"]);
        hScatterLenStatFirst->SetBinContent(counter, res["par0stat"]);
        hScatterLenStatFirst->SetBinError(counter, res["unc0stat"]);
        hScatterLenStatSecond->SetBinContent(counter, res["par1stat"]);
        hScatterLenStatSecond->SetBinError(counter, res["unc1stat"]);
        hChi2->SetBinContent(counter, res["chi2"] / res["ndf"]);
        hChi2->SetBinError(counter, 1.e-20);
        hScatterLenFirst->GetXaxis()->SetBinLabel(
            counter,
            Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
                 rad1, rad2, wei));
        hScatterLenSecond->GetXaxis()->SetBinLabel(
            counter,
            Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
                 rad1, rad2, wei));
        hScatterLenStatFirst->GetXaxis()->SetBinLabel(
            counter,
            Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
                 rad1, rad2, wei));
        hScatterLenStatSecond->GetXaxis()->SetBinLabel(
            counter,
            Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
                 rad1, rad2, wei));
        hChi2->GetXaxis()->SetBinLabel(
            counter,
            Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f",
                 rad1, rad2, wei));
      }
    }
  }
  double rad1 = 1.15;
  double rad2 = 2.69;
  double wei = 0.68;
  std::map<std::string, double> res = fitSimultaneousLL(
      gSCStat, gOCStat, gSC, gOC, rad1, rad2, wei, 200., kUpdate, deleteFile,
      Form("LLSimFit_Dpi%s", suffix.data()), "", true, fix);
  hScatterLenFirst->SetBinContent(28, res["par0"]);
  hScatterLenFirst->SetBinError(28, res["unc0"]);
  hScatterLenSecond->SetBinContent(28, res["par1"]);
  hScatterLenSecond->SetBinError(28, res["unc1"]);
  hScatterLenStatFirst->SetBinContent(28, res["par0stat"]);
  hScatterLenStatFirst->SetBinError(28, res["unc0stat"]);
  hScatterLenStatSecond->SetBinContent(28, res["par1stat"]);
  hScatterLenStatSecond->SetBinError(28, res["unc1stat"]);
  hChi2->SetBinContent(28, res["chi2"] / res["ndf"]);
  hChi2->SetBinError(28, 1.e-20);
  hScatterLenFirst->GetXaxis()->SetBinLabel(
      28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
               rad2, wei));
  hScatterLenSecond->GetXaxis()->SetBinLabel(
      28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
               rad2, wei));
  hScatterLenStatFirst->GetXaxis()->SetBinLabel(
      28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
               rad2, wei));
  hScatterLenStatSecond->GetXaxis()->SetBinLabel(
      28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
               rad2, wei));
  hChi2->GetXaxis()->SetBinLabel(
      28, Form("#it{r}_{1} = %0.2f #it{r}_{2} = %0.2f #it{w}_{1} = %0.2f", rad1,
               rad2, wei));
  gROOT->SetBatch(false);

  auto cScattLen = new TCanvas("cScattLen", "", 1500, 500);
  auto leg = new TLegend(0.2, 0.2, 0.4, 0.4);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hScatterLenFirst, "a_{I=3/2}");
  leg->AddEntry(hScatterLenSecond, "a_{I=1/2}");
  cScattLen->SetLeftMargin(0.08);
  hScatterLenFirst->GetYaxis()->SetRangeUser(-0.3, 0.1);
  hScatterLenFirst->Draw();
  hScatterLenSecond->Draw("same");
  leg->Draw();
  cScattLen->SaveAs("ScatteringLength_Dpi_fit_allradii.pdf");

  auto cChi2 = new TCanvas("cChi2", "", 1500, 500);
  cChi2->SetLeftMargin(0.08);
  hChi2->Draw();
  cChi2->SaveAs("ScatteringLength_Dpi_fitchi2_allradii.pdf");

  auto lineVert = new TLine(0., -0.29, 0., 0.29);
  lineVert->SetLineWidth(1);
  lineVert->SetLineColor(kGray + 1);
  lineVert->SetLineStyle(2);

  auto lineHor = new TLine(-0.29, 0., 0.29, 0.);
  lineHor->SetLineWidth(1);
  lineHor->SetLineColor(kGray + 1);
  lineHor->SetLineStyle(2);

  std::vector<TGraph *> gContour1Sigma{}, gContour2Sigma{};
  auto inFile = TFile::Open(Form("LLSimFit_Dpi%s.root", suffix.data()));
  TDirectoryFile *dir = nullptr;
  TKey *key;
  TIter next(inFile->GetListOfKeys());
  while ((key = (TKey *)next())) {
    dir = dynamic_cast<TDirectoryFile *>(inFile->Get(key->GetName()));
    gContour1Sigma.push_back(
        dynamic_cast<TGraph *>(dir->Get("gContour1Sigma")));
    gContour2Sigma.push_back(
        dynamic_cast<TGraph *>(dir->Get("gContour2Sigma")));
  }

  auto cContourAll = new TCanvas("cContourAll", "", 900, 900);
  cContourAll->DrawFrame(-0.29, -0.29, 0.29, 0.29,
                         ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
  lineHor->Draw();
  lineVert->Draw();

  for (int i{0}; i < 28; ++i) {
    gContour2Sigma[i]->Draw("cf");
  }
  for (int i{0}; i < 28; ++i) {
    gContour1Sigma[i]->Draw("cf");
  }

  auto legContour = new TLegend(0.7, 0.75, 0.9, 0.9);
  legContour->SetTextSize(0.045);
  legContour->SetBorderSize(0);
  legContour->SetFillStyle(0);
  legContour->AddEntry(gContour1Sigma[0], "68% CI", "f");
  legContour->AddEntry(gContour2Sigma[0], "95% CI", "f");
  legContour->Draw();

  auto lat = new TLatex();
  lat->SetTextColor(kBlack);
  lat->SetTextFont(42);
  lat->SetTextSize(0.07);
  lat->SetNDC();
  cContourAll->SaveAs("ContourPlot_Dpi_fit_allradii_convolution.pdf");

  auto cContour = new TCanvas("cContour", "", 1920, 1080);
  cContour->Divide(10, 3);
  int iPad{0};
  for (auto &rad1 : radiusFirst) {
    for (auto &rad2 : radiusSecond) {
      for (auto &wei : weightsFirst) {
        ++iPad;
        auto hFrame = cContour->cd(iPad)->DrawFrame(
            -0.29, -0.29, 0.29, 0.29, ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
        hFrame->GetXaxis()->SetTitleSize(0.07);
        hFrame->GetYaxis()->SetTitleSize(0.07);
        hFrame->GetXaxis()->SetLabelSize(0.07);
        hFrame->GetYaxis()->SetLabelSize(0.07);
        hFrame->GetXaxis()->SetTitleOffset(0.9);
        hFrame->GetYaxis()->SetTitleOffset(0.9);
        lineHor->Draw();
        lineVert->Draw();
        gContour2Sigma[iPad - 1]->Draw("cf");
        gContour1Sigma[iPad - 1]->Draw("cf");
        lat->DrawLatex(0.2, 0.9, Form("#it{r}_{1} = %0.2f fm", rad1));
        lat->DrawLatex(0.2, 0.8, Form("#it{r}_{2} = %0.2f fm", rad2));
        lat->DrawLatex(0.2, 0.7, Form("#it{w}_{1} = %0.2f fm", wei));
      }
    }
  }
  auto hFrame = cContour->cd(iPad + 1)->DrawFrame(
      -0.29, -0.29, 0.29, 0.29, ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
  hFrame->GetXaxis()->SetTitleSize(0.07);
  hFrame->GetYaxis()->SetTitleSize(0.07);
  hFrame->GetXaxis()->SetLabelSize(0.07);
  hFrame->GetYaxis()->SetLabelSize(0.07);
  hFrame->GetXaxis()->SetTitleOffset(0.9);
  hFrame->GetYaxis()->SetTitleOffset(0.9);
  lineHor->Draw();
  lineVert->Draw();
  gContour2Sigma[iPad]->Draw("cf");
  gContour1Sigma[iPad]->Draw("cf");
  lat->DrawLatex(0.2, 0.9, Form("#it{r}_{1} = %0.2f fm", rad1));
  lat->DrawLatex(0.2, 0.8, Form("#it{r}_{2} = %0.2f fm", rad2));
  lat->DrawLatex(0.2, 0.7, Form("#it{w}_{1} = %0.2f fm", wei));
  cContour->SaveAs("ContourPlot_Dpi_fit_allradii.pdf");

  // add some info to the output files
  auto fileToUpdate =
      TFile::Open(Form("LLSimFit_Dpi%s.root", suffix.data()), "update");
  fileToUpdate->cd();
  TDirectoryFile *dirScatt =
      new TDirectoryFile("scattering_length", "scattering_length");
  dirScatt->Write();
  dirScatt->cd();
  hScatterLenFirst->Write();
  hScatterLenSecond->Write();
  hScatterLenStatFirst->Write();
  hScatterLenStatSecond->Write();
  hChi2->Write();
  cScattLen->Write();
  cChi2->Write();
  dirScatt->Close();

  fileToUpdate->cd();
  TDirectoryFile *dirContour = new TDirectoryFile("contours", "contours");
  dirContour->Write();
  dirContour->cd();
  cContourAll->Write();
  cContour->Write();
  dirContour->Close();

  fileToUpdate->cd();
  TDirectoryFile *dirCoulomb = new TDirectoryFile("Coulomb", "Coulomb");
  dirCoulomb->Write();
  dirCoulomb->cd();
  gCoulombSCmin->Write();
  gCoulombSCmax->Write();
  gCoulombOCmin->Write();
  gCoulombOCmax->Write();
  dirCoulomb->Close();
  fileToUpdate->Close();
}

void mainBootstrapFit(int nIter = 2000) {
  gStyle->SetPalette(kDeepSea);

  auto inFileSC = TFile::Open("DPiPlusOutput_nolambda.root");
  auto inFileOC = TFile::Open("DPiMinusOutput_nolambda.root");

  auto hDistrScattLen =
      new TH2F("hDistrScattLen", ";a_{I=3/2} (fm);a_{I=1/2} (fm)", 50, -0.25,
               0.25, 50, -0.25, 0.25);
  bool deleteFile = true;
  double tot_time_elapsed = 0.;
  gROOT->SetBatch(true);
  std::cout << "\nPerform simultaneous fits" << std::endl;
  for (int iter{0}; iter < nIter; ++iter) {
    std::clock_t c_start = std::clock();
    auto gSC = (TGraphAsymmErrors *)inFileSC->FindObjectAny(
        Form("CF_corr_MC_truth_N_%d", iter));
    auto gOC = (TGraphAsymmErrors *)inFileOC->FindObjectAny(
        Form("CF_corr_MC_truth_N_%d", iter));
    if (iter > 0)
      deleteFile = false;
    std::map<std::string, double> resBootstrap = fitSimultaneousLL(
        gSC, gOC, gSC, gOC, 0.97, 2.52, 0.66, 200., kUpdate, deleteFile,
        "LLSimFit_Dpi_bootstrap", Form("_%d", iter));
    hDistrScattLen->Fill(resBootstrap["par0"], resBootstrap["par1"]);
    std::clock_t c_end = std::clock();
    double time_elapsed = 1.e0 * (c_end - c_start) / CLOCKS_PER_SEC;
    tot_time_elapsed += time_elapsed;
    std::cout << Form("\033[33m\rIteration %04d/%04d, time elapsed per "
                      "iteration = %0.2f (s), total time elapsed = %0.2f (s)",
                      iter + 1, nIter, time_elapsed, tot_time_elapsed)
              << std::flush;
  }
  std::cout << "\033[32m   -> Done\033[0m\n\n" << std::flush;
  gROOT->SetBatch(false);

  auto fGaus = new TF1("fGaus", "gaus", -0.25, 0.25);
  auto hDistrScattLenProjX = hDistrScattLen->ProjectionX("hDistrScattLenProjX");
  hDistrScattLenProjX->Fit("fGaus", "Q0");
  double meanFirst = fGaus->GetParameter(1);
  double sigmaFirst = fGaus->GetParameter(2);
  auto hDistrScattLenProjY = hDistrScattLen->ProjectionY("hDistrScattLenProjY");
  hDistrScattLenProjY->Fit("fGaus", "Q0");
  double meanSecond = fGaus->GetParameter(1);
  double sigmaSecond = fGaus->GetParameter(2);

  auto inFileSCCentral = TFile::Open("Results_PipDp_corrected.root");
  auto gSCCentral =
      (TGraphAsymmErrors *)inFileSCCentral->Get("CF_corr_N_truth");
  gSCCentral->SetName("gSC_corr_N_truth_stat_plus_sys");
  auto inFileOCCentral = TFile::Open("Results_PipDm_corrected.root");
  auto gOCCentral =
      (TGraphAsymmErrors *)inFileOCCentral->Get("CF_corr_N_truth");
  gOCCentral->SetName("gOC_corr_N_truth_stat_plus_sys");

  auto inFileSCCentralStat = TFile::Open("Result_PipDp.root");
  auto gSCCentralStat =
      (TGraphAsymmErrors *)inFileSCCentralStat->Get("CF_corr_N_truth");
  gSCCentralStat->SetName("gSC_corr_N_truth_stat");
  for (int iPt{0}; iPt < gSCCentralStat->GetN(); ++iPt) {
    double kStar, cf;
    gSCCentral->GetPoint(iPt, kStar, cf);
    gSCCentralStat->SetPoint(iPt, kStar, cf);
  }
  auto inFileOCCentralStat = TFile::Open("Result_PipDm.root");
  auto gOCCentralStat =
      (TGraphAsymmErrors *)inFileOCCentralStat->Get("CF_corr_N_truth");
  gOCCentralStat->SetName("gOC_corr_N_truth_stat");
  for (int iPt{0}; iPt < gOCCentralStat->GetN(); ++iPt) {
    double kStar, cf;
    gOCCentral->GetPoint(iPt, kStar, cf);
    gOCCentralStat->SetPoint(iPt, kStar, cf);
  }

  std::map<std::string, double> res = fitSimultaneousLL(
      gSCCentralStat, gOCCentralStat, gSCCentral, gOCCentral, 0.97, 2.52, 0.66,
      200., kRecreate, deleteFile, "LLSimFit_Dpi_ref4bootstrap", "", true);

  auto gResBootstrap = new TGraphAsymmErrors(1);
  gResBootstrap->SetPoint(0, meanFirst, meanSecond);
  gResBootstrap->SetPointError(0, sigmaFirst, sigmaFirst, sigmaSecond,
                               sigmaSecond);
  gResBootstrap->SetLineWidth(2);
  gResBootstrap->SetLineColor(kRed);
  gResBootstrap->SetMarkerColor(kRed);
  gResBootstrap->SetMarkerStyle(kFullCircle);

  auto gResFit = new TGraphAsymmErrors(1);
  gResFit->SetPoint(0, res["par0"], res["par1"]);
  gResFit->SetPointError(0, res["unc0"], res["unc0"], res["unc1"], res["unc1"]);
  gResFit->SetLineWidth(2);
  gResFit->SetLineColor(kRed + 1);
  gResFit->SetMarkerColor(kRed + 1);
  gResFit->SetMarkerStyle(kOpenCircle);

  gResBootstrap->SetName("gResBootstrap");
  gResFit->SetName("gResFit");

  auto cScattLen = new TCanvas("cScattLen", "", 800, 800);

  auto leg = new TLegend(0.5, 0.75, 0.8, 0.85);
  leg->SetTextSize(0.045);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(gResBootstrap, "bootstrap", "p");
  leg->AddEntry(gResFit, "single fit", "p");

  auto lineVert = new TLine(0., -0.2, 0., 0.35);
  lineVert->SetLineWidth(1);
  lineVert->SetLineColor(kGray + 1);
  lineVert->SetLineStyle(2);

  auto lineHor = new TLine(-0.2, 0., 0.35, 0.);
  lineHor->SetLineWidth(1);
  lineHor->SetLineColor(kGray + 1);
  lineHor->SetLineStyle(2);

  auto lineVertShort = new TLine(0., -0.25, 0., 0.15);
  lineVertShort->SetLineWidth(1);
  lineVertShort->SetLineColor(kGray + 1);
  lineVertShort->SetLineStyle(2);

  auto lineHorShort = new TLine(-0.25, 0., 0.15, 0.);
  lineHorShort->SetLineWidth(1);
  lineHorShort->SetLineColor(kGray + 1);
  lineHorShort->SetLineStyle(2);

  cScattLen->SetRightMargin(0.15);
  hDistrScattLen->GetXaxis()->SetDecimals();
  hDistrScattLen->GetXaxis()->SetNdivisions(505);
  hDistrScattLen->GetYaxis()->SetDecimals();
  hDistrScattLen->GetYaxis()->SetNdivisions(505);
  hDistrScattLen->Draw("colz");
  lineHor->Draw();
  lineVert->Draw();
  gResBootstrap->Draw("pz");
  gResFit->Draw("pz");
  leg->Draw();

  cScattLen->SaveAs("ScattLen_Dpi_fit_vs_bootstrap.pdf");

  // compute <Chi2> for different bootstrap iterations varying the two
  // parameters
  const int nQuantiles = 8;
  double prob4Quant[nQuantiles] = {0.01, 0.02, 0.03, 0.04,
                                   0.05, 0.1,  0.3,  0.5}; // 1sigma and 2sigma
  TH2F *hDistrMeanChi2 = nullptr;
  TH2F *hDistrMinChi2 = nullptr;
  TH2F *hDistrMostProbChi2 = nullptr;
  TH2F *hDistrQuantChi2[nQuantiles];

  if (gSystem->Exec("ls LLSimFit_Dpi_contour_bootstrap_vs_fit.root")) {
    double massD = TDatabasePDG::Instance()->GetParticle(411)->Mass() * 1000;
    double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
    double redMass = massD * massPi / (massD + massPi);
    auto fSCSim =
        new TF1("fSCSim", GeneralCoulombLednickyTwoRadii, 0., 200., nParsSC);
    auto fOCSim = new TF1("fOCSim", GeneralCoulombLednickySecondTwoRadii, 0.,
                          200., nParsOC);
    double parsSC[nParsSC] = {0.97, 2.52, 0.66, -0.08, 0., 0., 1., redMass};
    double parsOC[nParsOC] = {0.97,  2.52, 0.66, -0.08, 0.,
                              -0.08, 0.,   0.,   -1.,   redMass};
    fSCSim->SetParameters(parsSC);
    fOCSim->SetParameters(parsOC);
    const int nBinsX = 40;
    const int nBinsY = 40;
    hDistrMeanChi2 = new TH2F("hDistrMeanChi2",
                              ";a_{I=3/2} (fm);a_{I=1/2} (fm);#LT#chi^{2}#GT",
                              nBinsX, -0.25, 0.15, nBinsY, -0.25, 0.15);
    hDistrMinChi2 =
        new TH2F("hDistrMinChi2", ";a_{I=3/2} (fm);a_{I=1/2} (fm);min #chi^{2}",
                 nBinsX, -0.25, 0.15, nBinsY, -0.25, 0.15);
    hDistrMostProbChi2 =
        new TH2F("hDistrMostProbChi2",
                 ";a_{I=3/2} (fm);a_{I=1/2} (fm);most prob #chi^{2}", nBinsX,
                 -0.25, 0.15, nBinsY, -0.25, 0.15);
    for (int iQ{0}; iQ < nQuantiles; ++iQ) {
      hDistrQuantChi2[iQ] = new TH2F(
          Form("hDistrQuant%0.fPercChi2", prob4Quant[iQ] * 100),
          Form(";a_{I=3/2} (fm);a_{I=1/2} (fm);%0.f perc quantile #chi^{2}",
               prob4Quant[iQ] * 100),
          nBinsX, -0.25, 0.15, nBinsY, -0.25, 0.15);
    }
    std::cout << "\nCompute contour plot" << std::endl;
    tot_time_elapsed = 0.;
    auto hChi2 = new TH1F("hChi2", "", 1000, 0., 500.);
    TF1 *fLandau = new TF1("fLandau", "landau", 0., 500);
    fLandau->SetParLimits(0, 0., nIter);
    fLandau->SetParLimits(1, 0., 1000.);
    fLandau->SetParLimits(2, 0., 1000.);
    int iIter = 0;
    for (int iPar0{0}; iPar0 < nBinsX; ++iPar0) {
      double par0 = hDistrMeanChi2->GetXaxis()->GetBinCenter(iPar0 + 1);
      fSCSim->SetParameter(3, par0);
      fOCSim->SetParameter(3, par0);
      for (int iPar1{0}; iPar1 < nBinsY; ++iPar1) {
        iIter++;
        double par1 = hDistrMeanChi2->GetYaxis()->GetBinCenter(iPar1 + 1);
        fOCSim->SetParameter(5, par1);
        std::clock_t c_start = std::clock();
        for (int iter{0}; iter < nIter; ++iter) {
          auto gSC = (TGraphAsymmErrors *)inFileSC->FindObjectAny(
              Form("CF_corr_MC_truth_N_%d", iter));
          auto gOC = (TGraphAsymmErrors *)inFileOC->FindObjectAny(
              Form("CF_corr_MC_truth_N_%d", iter));
          std::array<double, 2> chi2 = ComputeChi2(gSC, gOC, fSCSim, fOCSim);
          hChi2->Fill(chi2[0]);
        }
        std::clock_t c_end = std::clock();
        double time_elapsed = 1.e0 * (c_end - c_start) / CLOCKS_PER_SEC;
        tot_time_elapsed += time_elapsed;
        std::cout
            << Form(
                   "\033[33m\ra3/2 = %0.4f, a1/2 = %0.4f (%0.f%%) time elapsed "
                   "per iteration = %0.2f (s), total time elapsed = %0.2f (s)",
                   par0, par1, double(iIter) * 100 / (nBinsX * nBinsY),
                   time_elapsed, tot_time_elapsed)
            << std::flush;
        // double minVal = GetMinHisto(hChi2);
        double meanVal = hChi2->GetMean();
        hChi2->Fit("fLandau", "Q0L");
        TF1 fLandauSet("fLandauSet", "landau", 0., 500);
        fLandauSet.SetParameters(fLandau->GetParameters());
        double mostProbVal = fLandauSet.GetMaximumX(0., 500.);
        double quantiles[nQuantiles];
        hChi2->GetQuantiles(nQuantiles, quantiles, prob4Quant);
        hDistrMeanChi2->SetBinContent(iPar0 + 1, iPar1 + 1, meanVal);
        hDistrMinChi2->SetBinContent(iPar0 + 1, iPar1 + 1, quantiles[4]);
        hDistrMostProbChi2->SetBinContent(iPar0 + 1, iPar1 + 1, mostProbVal);
        for (int iQ{0}; iQ < nQuantiles; ++iQ) {
          hDistrQuantChi2[iQ]->SetBinContent(iPar0 + 1, iPar1 + 1,
                                             quantiles[iQ]);
        }
        hChi2->Reset();
      }
    }
    std::cout << "\033[32m   -> Done\033[0m\n\n" << std::flush;
  } else {
    auto inFileDistr =
        TFile::Open("LLSimFit_Dpi_contour_bootstrap_vs_fit.root");
    hDistrMeanChi2 = (TH2F *)inFileDistr->Get("hDistrMeanChi2");
    hDistrMinChi2 = (TH2F *)inFileDistr->Get("hDistrMinChi2");
    hDistrMostProbChi2 = (TH2F *)inFileDistr->Get("hDistrMostProbChi2");
    for (int iQ{0}; iQ < nQuantiles; ++iQ) {
      hDistrQuantChi2[iQ] = (TH2F *)inFileDistr->Get(
          Form("hDistrQuant%0.fPercChi2", prob4Quant[iQ] * 100));
    }
  }

  // build contout plot
  auto gContour1SigmaBootstrapMeanChi2 = new TGraph(1);
  auto gContour2SigmaBootstrapMeanChi2 = new TGraph(1);
  gContour1SigmaBootstrapMeanChi2->SetName("gContour1SigmaBootstrapMeanChi2");
  gContour2SigmaBootstrapMeanChi2->SetName("gContour2SigmaBootstrapMeanChi2");
  gContour1SigmaBootstrapMeanChi2->SetLineColor(kRed);
  gContour1SigmaBootstrapMeanChi2->SetMarkerColor(kRed);
  gContour1SigmaBootstrapMeanChi2->SetLineWidth(2);
  gContour2SigmaBootstrapMeanChi2->SetLineColor(kRed + 2);
  gContour2SigmaBootstrapMeanChi2->SetMarkerColor(kRed + 2);
  gContour2SigmaBootstrapMeanChi2->SetLineWidth(2);
  auto gContour1SigmaBootstrapMinChi2 = new TGraph(1);
  auto gContour2SigmaBootstrapMinChi2 = new TGraph(1);
  gContour1SigmaBootstrapMinChi2->SetName("gContour1SigmaBootstrapMinChi2");
  gContour2SigmaBootstrapMinChi2->SetName("gContour2SigmaBootstrapMinChi2");
  gContour1SigmaBootstrapMinChi2->SetLineColor(kGreen + 2);
  gContour1SigmaBootstrapMinChi2->SetMarkerColor(kGreen + 2);
  gContour1SigmaBootstrapMinChi2->SetLineWidth(2);
  gContour2SigmaBootstrapMinChi2->SetLineColor(kGreen + 3);
  gContour2SigmaBootstrapMinChi2->SetMarkerColor(kGreen + 3);
  gContour2SigmaBootstrapMinChi2->SetLineWidth(2);

  double minMeanChi2 = 1.e20;
  double minMinChi2 = 1.e20;
  auto fPol2 = new TF1("fPol2", "[0] + [1]*x + [2] * x * x", -0.25, 0.25);
  for (int iBin{1}; iBin <= hDistrMeanChi2->GetYaxis()->GetNbins(); ++iBin) {
    auto hProjMean = hDistrMeanChi2->ProjectionX("hProjMean", iBin, iBin);
    hProjMean->Fit(fPol2, "Q0L", "", -0.2, 0.1);
    double c = fPol2->GetParameter(0);
    double b = fPol2->GetParameter(1);
    double a = fPol2->GetParameter(2);
    auto minPol2 = -b / (2 * a);
    if (fPol2->Eval(minPol2) < minMeanChi2)
      minMeanChi2 = fPol2->Eval(minPol2);
    delete hProjMean;

    auto hProjMin = hDistrMinChi2->ProjectionX("hProjMin", iBin, iBin);
    hProjMin->Fit(fPol2, "Q0L", "", -0.2, 0.1);
    c = fPol2->GetParameter(0);
    b = fPol2->GetParameter(1);
    a = fPol2->GetParameter(2);
    minPol2 = -b / (2 * a);
    if (fPol2->Eval(minPol2) < minMinChi2)
      minMinChi2 = fPol2->Eval(minPol2);
    delete hProjMin;
  }

  int iPt1Sigma{0}, iPt2Sigma{0};
  std::vector<double> x1, x2, x3, x4, y1, y2;
  for (int iBin{1}; iBin <= hDistrMeanChi2->GetYaxis()->GetNbins(); ++iBin) {
    auto hProj = hDistrMeanChi2->ProjectionX("hProj", iBin, iBin);
    hProj->Fit(fPol2, "Q0L", "", -0.2, 0.1);
    double c = fPol2->GetParameter(0);
    double b = fPol2->GetParameter(1);
    double a = fPol2->GetParameter(2);
    auto minPol2 = -b / (2 * a);
    if (fPol2->Eval(minPol2) > minMeanChi2 + 4)
      continue;
    if (fPol2->Eval(minPol2) < minMeanChi2 + 1) {
      y1.push_back(hDistrMeanChi2->GetYaxis()->GetBinCenter(iBin));
      x1.push_back((-b + std::sqrt(b * b - 4 * a * (c - (minMeanChi2 + 1)))) /
                   (2 * a));
      x2.push_back((-b - std::sqrt(b * b - 4 * a * (c - (minMeanChi2 + 1)))) /
                   (2 * a));
    }
    y2.push_back(hDistrMeanChi2->GetYaxis()->GetBinCenter(iBin));
    x3.push_back((-b + std::sqrt(b * b - 4 * a * (c - (minMeanChi2 + 4)))) /
                 (2 * a));
    x4.push_back((-b - std::sqrt(b * b - 4 * a * (c - (minMeanChi2 + 4)))) /
                 (2 * a));
    delete hProj;
  }

  std::vector<double> x1min, x2min, x3min, x4min, y1min, y2min;
  for (int iBin{1}; iBin <= hDistrMinChi2->GetYaxis()->GetNbins(); ++iBin) {
    auto hProj = hDistrMinChi2->ProjectionX("hProj", iBin, iBin);
    hProj->Fit(fPol2, "Q0L", "", -0.2, 0.1);
    double c = fPol2->GetParameter(0);
    double b = fPol2->GetParameter(1);
    double a = fPol2->GetParameter(2);
    auto minPol2 = -b / (2 * a);
    if (fPol2->Eval(minPol2) > minMinChi2 + 4)
      continue;
    if (fPol2->Eval(minPol2) < minMinChi2 + 1) {
      y1min.push_back(hDistrMinChi2->GetYaxis()->GetBinCenter(iBin));
      x1min.push_back((-b + std::sqrt(b * b - 4 * a * (c - (minMinChi2 + 1)))) /
                      (2 * a));
      x2min.push_back((-b - std::sqrt(b * b - 4 * a * (c - (minMinChi2 + 1)))) /
                      (2 * a));
    }
    y2min.push_back(hDistrMinChi2->GetYaxis()->GetBinCenter(iBin));
    x3min.push_back((-b + std::sqrt(b * b - 4 * a * (c - (minMinChi2 + 4)))) /
                    (2 * a));
    x4min.push_back((-b - std::sqrt(b * b - 4 * a * (c - (minMinChi2 + 4)))) /
                    (2 * a));
    delete hProj;
  }

  for (size_t iPt{0}; iPt < x1.size(); ++iPt) {
    gContour1SigmaBootstrapMeanChi2->SetPoint(iPt, x1[iPt], y1[iPt]);
  }
  for (size_t iPt{0}; iPt < x2.size(); ++iPt) {
    gContour1SigmaBootstrapMeanChi2->SetPoint(
        iPt + x1.size(), x2[x2.size() - 1 - iPt], y1[x2.size() - 1 - iPt]);
  }
  gContour1SigmaBootstrapMeanChi2->SetPoint(
      gContour1SigmaBootstrapMeanChi2->GetN(), x1[0], y1[0]);
  for (size_t iPt{0}; iPt < x3.size(); ++iPt) {
    gContour2SigmaBootstrapMeanChi2->SetPoint(iPt, x3[iPt], y2[iPt]);
  }
  for (size_t iPt{0}; iPt < x4.size(); ++iPt) {
    gContour2SigmaBootstrapMeanChi2->SetPoint(
        iPt + x3.size(), x4[x4.size() - 1 - iPt], y2[x4.size() - 1 - iPt]);
  }
  gContour2SigmaBootstrapMeanChi2->SetPoint(
      gContour2SigmaBootstrapMeanChi2->GetN(), x3[0], y2[0]);

  for (size_t iPt{0}; iPt < x1min.size(); ++iPt) {
    gContour1SigmaBootstrapMinChi2->SetPoint(iPt, x1min[iPt], y1min[iPt]);
  }
  for (size_t iPt{0}; iPt < x2min.size(); ++iPt) {
    gContour1SigmaBootstrapMinChi2->SetPoint(iPt + x1min.size(),
                                             x2min[x2min.size() - 1 - iPt],
                                             y1min[x2min.size() - 1 - iPt]);
  }
  gContour1SigmaBootstrapMinChi2->SetPoint(
      gContour1SigmaBootstrapMinChi2->GetN(), x1min[0], y1min[0]);
  for (size_t iPt{0}; iPt < x3min.size(); ++iPt) {
    gContour2SigmaBootstrapMinChi2->SetPoint(iPt, x3min[iPt], y2min[iPt]);
  }
  for (size_t iPt{0}; iPt < x4min.size(); ++iPt) {
    gContour2SigmaBootstrapMinChi2->SetPoint(iPt + x3min.size(),
                                             x4min[x4min.size() - 1 - iPt],
                                             y2min[x4min.size() - 1 - iPt]);
  }
  gContour2SigmaBootstrapMinChi2->SetPoint(
      gContour2SigmaBootstrapMinChi2->GetN(), x3min[0], y2min[0]);

  auto cContour = new TCanvas("cContour", "", 1000, 1000);
  cContour->Divide(2, 2);
  cContour->cd(1)->SetLeftMargin(0.15);
  cContour->cd(1)->SetRightMargin(0.15);
  cContour->cd(1)->SetLogz();
  hDistrMeanChi2->GetXaxis()->SetDecimals();
  hDistrMeanChi2->GetXaxis()->SetNdivisions(505);
  hDistrMeanChi2->GetYaxis()->SetDecimals();
  hDistrMeanChi2->GetYaxis()->SetNdivisions(505);
  hDistrMeanChi2->GetZaxis()->SetMoreLogLabels();
  hDistrMeanChi2->Draw("colz");
  gContour2SigmaBootstrapMeanChi2->Draw("same");
  gContour1SigmaBootstrapMeanChi2->Draw("same");
  lineHorShort->Draw();
  lineVertShort->Draw();
  cContour->cd(2)->SetLeftMargin(0.15);
  cContour->cd(2)->SetRightMargin(0.15);
  cContour->cd(2)->SetLogz();
  hDistrMinChi2->GetXaxis()->SetDecimals();
  hDistrMinChi2->GetXaxis()->SetNdivisions(505);
  hDistrMinChi2->GetYaxis()->SetDecimals();
  hDistrMinChi2->GetYaxis()->SetNdivisions(505);
  hDistrMinChi2->GetZaxis()->SetMoreLogLabels();
  hDistrMinChi2->Draw("colz");
  gContour2SigmaBootstrapMinChi2->Draw("same");
  gContour1SigmaBootstrapMinChi2->Draw("same");
  lineHorShort->Draw();
  lineVertShort->Draw();

  auto inFileCentral = TFile::Open("LLSimFit_Dpi_ref4bootstrap_radiusFirst0.97_"
                                   "radiusSecond2.52_weightFirst0.66_.root");
  auto gContour1Sigma =
      (TGraphAsymmErrors *)inFileCentral->Get("gContour1Sigma");
  auto gContour2Sigma =
      (TGraphAsymmErrors *)inFileCentral->Get("gContour2Sigma");
  auto gContour1SigmaStat =
      (TGraphAsymmErrors *)inFileCentral->Get("gContour1SigmaStat");
  auto gContour2SigmaStat =
      (TGraphAsymmErrors *)inFileCentral->Get("gContour2SigmaStat");
  auto legContour = new TLegend(0.2, 0.6, 0.4, 0.9);
  legContour->SetTextSize(0.04);
  legContour->SetBorderSize(0);
  legContour->SetFillStyle(0);
  legContour->AddEntry(gContour1Sigma, "68% CI fit (#sqrt{stat^{2}+syst^{2}})",
                       "f");
  legContour->AddEntry(gContour2Sigma, "95% CI fit (#sqrt{stat^{2}+syst^{2}})",
                       "f");
  legContour->AddEntry(gContour1SigmaBootstrapMeanChi2,
                       "68% CI bootstrap (#LT #chi^{2} #GT)", "f");
  legContour->AddEntry(gContour2SigmaBootstrapMeanChi2,
                       "95% CI bootstrap (#LT #chi^{2} #GT)", "f");
  legContour->AddEntry(gContour1SigmaBootstrapMinChi2,
                       "68% CI bootstrap (min #chi^{2})", "f");
  legContour->AddEntry(gContour2SigmaBootstrapMinChi2,
                       "95% CI bootstrap (min #chi^{2})", "f");
  auto legContourStat = new TLegend(0.2, 0.6, 0.4, 0.9);
  legContourStat->SetTextSize(0.04);
  legContourStat->SetBorderSize(0);
  legContourStat->SetFillStyle(0);
  legContourStat->AddEntry(gContour1SigmaStat, "68% CI fit (stat)", "f");
  legContourStat->AddEntry(gContour2SigmaStat, "95% CI fit (stat)", "f");
  legContourStat->AddEntry(gContour1SigmaBootstrapMeanChi2,
                           "68% CI bootstrap (#LT #chi^{2} #GT)", "f");
  legContourStat->AddEntry(gContour2SigmaBootstrapMeanChi2,
                           "95% CI bootstrap (#LT #chi^{2} #GT)", "f");
  legContourStat->AddEntry(gContour1SigmaBootstrapMinChi2,
                           "68% CI bootstrap (min #chi^{2})", "f");
  legContourStat->AddEntry(gContour2SigmaBootstrapMinChi2,
                           "95% CI bootstrap (min #chi^{2})", "f");

  cContour->cd(3)->SetLeftMargin(0.15);
  auto hFrame = cContour->cd(3)->DrawFrame(-0.2, -0.2, 0.35, 0.35,
                                           ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
  hFrame->GetXaxis()->SetDecimals();
  hFrame->GetXaxis()->SetNdivisions(505);
  hFrame->GetYaxis()->SetDecimals();
  hFrame->GetYaxis()->SetNdivisions(505);
  gContour2Sigma->Draw("cf");
  gContour1Sigma->Draw("cf");
  gContour2SigmaBootstrapMeanChi2->Draw("same");
  gContour1SigmaBootstrapMeanChi2->Draw("same");
  gContour2SigmaBootstrapMinChi2->Draw("same");
  gContour1SigmaBootstrapMinChi2->Draw("same");
  legContour->Draw();
  lineHor->Draw();
  lineVert->Draw();

  cContour->cd(4)->SetLeftMargin(0.15);
  auto hFrameStat = cContour->cd(4)->DrawFrame(
      -0.2, -0.2, 0.35, 0.35, ";a_{I=3/2} (fm);a_{I=1/2} (fm)");
  hFrameStat->GetXaxis()->SetDecimals();
  hFrameStat->GetXaxis()->SetNdivisions(505);
  hFrameStat->GetYaxis()->SetDecimals();
  hFrameStat->GetYaxis()->SetNdivisions(505);
  gContour2SigmaStat->Draw("cf");
  gContour1SigmaStat->Draw("cf");
  gContour2SigmaBootstrapMeanChi2->Draw("same");
  gContour1SigmaBootstrapMeanChi2->Draw("same");
  gContour2SigmaBootstrapMinChi2->Draw("same");
  gContour1SigmaBootstrapMinChi2->Draw("same");
  legContourStat->Draw();
  lineHor->Draw();
  lineVert->Draw();

  cContour->SaveAs("ContourPlot_Dpi_fit_vs_bootstrap.pdf");

  TFile outFileContour("LLSimFit_Dpi_contour_bootstrap_vs_fit.root",
                       "recreate");
  hDistrScattLen->Write();
  gResBootstrap->Write();
  gResFit->Write();
  hDistrMeanChi2->Write();
  hDistrMinChi2->Write();
  hDistrMostProbChi2->Write();
  for (int iQ{0}; iQ < nQuantiles; ++iQ) {
    hDistrQuantChi2[iQ]->Write();
  }
  gContour2Sigma->Write();
  gContour1Sigma->Write();
  gContour2SigmaStat->Write();
  gContour1SigmaStat->Write();
  gContour2SigmaBootstrapMeanChi2->Write();
  gContour1SigmaBootstrapMeanChi2->Write();
  gContour2SigmaBootstrapMinChi2->Write();
  gContour1SigmaBootstrapMinChi2->Write();
  outFileContour.Close();
}
