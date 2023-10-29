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
#include <TROOT.h>

#include "./functions.h"

// singlet = 3/2
// definition of shared parameter
// same charge function
const int nParsSC = 9;
int iparSC[nParsSC] = {
    0,  // first source radius
    1,  // second source radius
    2,  // normalisation of two contributions
    3,  // scattering length singlet
    4,  // effective range singlet
    7,  // QS
    8,  // charge product
    10,  // reduced mass
    11,  // overall normalization
};

// opposite charge function
const int nParsOC = 11;
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
    10,  // reduced mass
    12,  // overall normalization
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
                                       bool scalableLL = false, bool deleteFile = false,
                                       std::string outName = "LLSimFit_Dpi",
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

    double massD = TDatabasePDG::Instance()->GetParticle(413)->Mass() * 1000;
    double massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass() * 1000;
    double redMass = massD * massPi / (massD + massPi);

    // now perform global fit
    auto fSCSim = new TF1("fSCSim", ScalableGeneralCoulombLednickyTwoRadii, 0., maxRange, nParsSC);
    auto fOCSim = new TF1("fOCSim", ScalableGeneralCoulombLednickySecondTwoRadii, 0., maxRange, nParsOC);

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

    const int nPars = 13;
    double par0[nPars] = {radiusFirst, radiusSecond, weightFirst, -0.08, 0., -0.08, 0., 0., 1., -1., redMass, 1, 1};

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
    if (scalableLL) {
        fitter.Config().ParSettings(11).SetLowerLimit(0.95);
        fitter.Config().ParSettings(11).SetUpperLimit(1.05);
        fitter.Config().ParSettings(12).SetLowerLimit(0.95);
        fitter.Config().ParSettings(12).SetUpperLimit(1.05);
    } else {
        fitter.Config().ParSettings(11).Fix();
        fitter.Config().ParSettings(12).Fix();
    }

    fitter.Config().MinimizerOptions().SetPrintLevel(1);
    fitter.Config().SetMinimizer("Minuit2", "Migrad");

    // // fit FCN function directly
    // // (specify optionally data size and flag to indicate that is a chi2 fit)
    // // repeat with only stat unc
    fitter.FitFCN(nPars, globalChi2, 0, dataSC.Size() + dataOC.Size(), true);
    ROOT::Fit::FitResult result = fitter.Result();

    return result;
}

void fitSimultaneousLLDpi(const char *inFileName, const char *oFileName, int nIter, double sourceRelUnc = -1,
                          bool scalableLL = false) {
    gROOT->SetBatch(true);

    TFile *inFile = new TFile(inFileName);
    TFile *oFile = new TFile(Form("%s.root", oFileName), "create");
    if (oFile->IsZombie()) {  // file already exists
        return;
    }

    for (auto uncType : {"tot", "stat"}) {
        oFile->mkdir(Form("%s/iters", uncType));
        oFile->cd(Form("%s/iters", uncType));
        TNtuple *tResults = new TNtuple("tResults", "", "a0sin:a0tri:chi2ndf:r1:r2:w1:normsin:normtri");

        TCanvas *cCombFit = new TCanvas("cCombFit", "", 1200, 600);
        cCombFit->Divide(3, 1);
        cCombFit->SaveAs(Form("%s.pdf[", oFileName));
            // Form("/home/daniel/an/DstarPi/20_luuksel/GenCFCorr_nopc_kStarBW50MeV_bs%d%s.pdf[", nIter, uncType));

        std::vector<double> weights1 = {0.66, 0.69, 0.64};
        std::vector<double> radii1 = {0.97};
        std::vector<double> radii2 = {2.52};

        if (sourceRelUnc < 0) {  // use normal uncertainties extracted from the source fit
            radii1.push_back(1.06);
            radii1.push_back(0.89);

            radii2.push_back(2.88);
            radii2.push_back(2.32);
        } else {  // use source from extrap of r_core(mult)
            radii1.push_back(radii1[0] * (1 - sourceRelUnc));
            radii1.push_back(radii1[0] * (1 + sourceRelUnc));

            radii2.push_back(radii2[0] * (1 - sourceRelUnc));
            radii2.push_back(radii2[0] * (1 + sourceRelUnc));
        }

        // std::vector<std::vector<double>> fitRanges = {{10, 450}, {10, 400}, {10, 500}};
        std::vector<std::vector<double>> fitRanges = {{10, 350}, {10, 300}};

        for (int iIter = 0; iIter < nIter; iIter++) {
            TGraphAsymmErrors *gSC =
                reinterpret_cast<TGraphAsymmErrors *>(inFile->Get(Form("sc/%s/gCFGen%d", uncType, iIter)));
            gSC->SetName(Form("gSC%s%d", uncType, iIter));
            gSC->Write();
            TGraphAsymmErrors *gOC =
                reinterpret_cast<TGraphAsymmErrors *>(inFile->Get(Form("oc/%s/gCFGen%d", uncType, iIter)));
            gOC->SetName(Form("gOC%s%d", uncType, iIter));
            gOC->Write();

            int iRadius = 0;
            int iFitRange = 0;
            if (uncType == "tot") {
                iRadius = gRandom->Integer(3);
                iFitRange = gRandom->Integer(fitRanges.size());
            }

            ROOT::Fit::FitResult result = fitSimultaneousLL(gSC, gOC, radii1[iRadius], radii2[iRadius],
                                                            weights1[iRadius], fitRanges[iFitRange][1], scalableLL);
            const double *params = result.GetParams();
            const double *parErrors = result.GetErrors();

            double a0singlet = params[3];
            double a0triplet = params[5];
            double a0singletUnc = parErrors[3];
            double a0tripletUnc = parErrors[5];

            double chi2ndf = result.Chi2() / result.Ndf();
            if (std::abs(a0singlet) < 1 && std::abs(a0triplet) < 1 && chi2ndf < 10)
                tResults->Fill(a0singlet, a0triplet, chi2ndf, radii1[iRadius], radii2[iRadius], weights1[iRadius],
                               params[11], params[12]);

            // Draw
            TVirtualPad *pad1 = cCombFit->cd(1);
            pad1->DrawFrame(0, 0.8, 300, 1.4);
            auto fSCSim = new TF1("fSCSim", ScalableGeneralCoulombLednickyTwoRadii, 0., 300, nParsSC);

            fSCSim->FixParameter(0, params[0]);
            fSCSim->FixParameter(1, params[1]);
            fSCSim->FixParameter(2, params[2]);
            fSCSim->FixParameter(3, params[3]);
            fSCSim->FixParameter(4, params[4]);
            fSCSim->FixParameter(5, params[7]);
            fSCSim->FixParameter(6, params[8]);
            fSCSim->FixParameter(7, params[10]);
            fSCSim->FixParameter(8, params[11]);

            TVirtualPad *pad2 = cCombFit->cd(2);
            pad2->DrawFrame(0, 0.8, 300, 1.4);
            auto fOCSim = new TF1("fSCSim", ScalableGeneralCoulombLednickySecondTwoRadii, 0., 300, nParsOC);
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
            fOCSim->FixParameter(10, params[12]);
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
            cCombFit->SaveAs(Form("%s.pdf", oFileName));
            delete gSC;
            delete gOC;
        }
        cCombFit->SaveAs(Form("%s.pdf]", oFileName));
        oFile->cd(uncType);
        tResults->Write();
        delete tResults;
        delete cCombFit;
    }
    oFile->Close();
    printf(Form("Output saved in %s.root\n", oFileName));
}
