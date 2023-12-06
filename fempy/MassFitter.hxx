#ifndef FEMPY_MASSFITTER_HXX_
#define FEMPY_MASSFITTER_HXX_

#include <map>
#include <tuple>
#include <string>

#include "../../../alice/AliPhysics/PWG/Tools/yaml-cpp/include/yaml-cpp/yaml.h"


double Gaus(double *x, double *par) {
    double norm = 1. / TMath::Sqrt((2. * TMath::Pi())) / par[2];
    return norm * par[0] * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2]);
}

double Hat(double *x, double *par) {
    // p0: total yield
    // p1: mean
    // p2: sigma of thin gaussian
    // p3: fraction of narrow gaussian yield
    // p4: wide/narrow gaussian width

    double normThin = 1. / TMath::Sqrt((2. * TMath::Pi())) / par[2];
    double gThin = normThin * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / par[2] / par[2]);

    double wideSigma = par[2] * par[4];
    double normWide = 1. / TMath::Sqrt((2. * TMath::Pi())) / wideSigma;
    double gWide = normWide * TMath::Exp(-(x[0] - par[1]) * (x[0] - par[1]) / 2. / wideSigma / wideSigma);

    return par[0] * (par[3] * gThin + (1 - par[3]) * gWide);
}

double PowEx(double *x, double *par) {
    // p0: total yield
    // p1: slope
    Double_t mpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

    if (x[0] < mpi) return 0;
    return par[0] * TMath::Sqrt(x[0] - mpi) * TMath::Exp(-1. * par[1] * (x[0] - mpi));
}

double Pol1(double *x, double *par) { return par[0] + par[1] * x[0]; }

double Pol5(double *x, double *par) { return par[0] + 
                                             par[1] * (x[0]-par[6]) + 
                                             par[2] * (x[0]-par[6]) * (x[0]-par[6]) +
                                             par[3] * (x[0]-par[6]) * (x[0]-par[6]) * (x[0]-par[6]) +
                                             par[4] * (x[0]-par[6]) * (x[0]-par[6]) * (x[0]-par[6]) * (x[0]-par[6]) + 
                                             par[5] * (x[0]-par[6]) * (x[0]-par[6]) * (x[0]-par[6]) * (x[0]-par[6]) * (x[0]-par[6]); }

double Exp(double *x, double *par) {
    // p0: total yield
    // p1: slope
    return par[0] * TMath::Exp(par[1] * x[0]);
}

double Spline5(double *x, double *par){
    int numKnots = 8;
    Double_t xKnots[numKnots];
    Double_t yKnots[numKnots];
    for(int iKnot=0; iKnot<numKnots; iKnot++){
        xKnots[iKnot] = par[iKnot];
        yKnots[iKnot] = par[numKnots+iKnot];
    }
    TSpline5* sp5 = new TSpline5("sp5", xKnots, yKnots, numKnots, "");
    return sp5->Eval(x[0]);
}

class MassFitter {
 public:
    std::map<std::string, int> nPars = {{"ngaus", 3}, {"gaus", 3}, {"hat", 5}, 
                                        {"pol1", 2}, {"pol5", 7}, {"powex", 2}, {"exp", 2}, {"spline5", 16}}; 

    enum SgnFuncs { kGaus = 0 }; // not used, delete?
    enum BkgFuncs { kPol1 = 0 }; // not used, delete?
    MassFitter(TH1 *hist, std::string sgnFuncName, std::string bkgFuncName, double fitRangeMin, double fitRangeMax,
               std::string configPath) {
        this->fHist = reinterpret_cast<TH1 *>(hist->Clone());
        this->fFitRangeMin = fitRangeMin;
        this->fFitRangeMax = fitRangeMax;
        this->fSgnFuncName = sgnFuncName;
        this->fBkgFuncName = bkgFuncName;

        this->fNSgnPars = nPars[this->fSgnFuncName];
        this->fNBkgPars = nPars[this->fBkgFuncName];
        int nTotPars = this->fNSgnPars + this->fNBkgPars; 

        this->fConfigPath = configPath.data();
        this->fCfgFile = YAML::LoadFile(configPath.data());

        if (fSgnFuncName == "gaus") {
            this->fSgnFunc = Gaus;
        } else if (fSgnFuncName == "hat") {
            this->fSgnFunc = Hat;
        } else {
            printf("Function not implemented\n");
            exit(1);
        }

        if (fBkgFuncName == "powex") {
            this->fBkgFunc = PowEx;
        } else if (fBkgFuncName == "exp") {
            this->fBkgFunc = Exp;
        } else if (fBkgFuncName == "pol1") {
            this->fBkgFunc = Pol1;
        } else if (fBkgFuncName == "pol5") {
            this->fBkgFunc = Pol5;
        } else if (fBkgFuncName == "spline5") {
            this->fBkgFunc = Spline5;    
        } else {
            printf("Function not implemented\n");
            exit(1);
        }

        this->fFit = new TF1(
            "fTot",
            [&, this](double *x, double *pars) {
                return this->fSgnFunc(x, pars) + this->fBkgFunc(x, &pars[this->fNSgnPars]);
            },
            fFitRangeMin, fFitRangeMax, nTotPars);
        fFit->SetNpx(300);

        this->fPrefit = new TF1(
            "fPrefit",
            [&, this, bkgFuncName](double *x, double *pars) -> double {
                if (fBkgFuncName == "powex" && std::abs(x[0] - 0.1455) < 0.003) {
                    TF1::RejectPoint();
                    return this->fBkgFunc(x, pars);
                } else if (fBkgFuncName == "exp" && std::abs(x[0] - 1.87) < 0.05) {
                    TF1::RejectPoint();
                    return this->fBkgFunc(x, pars);
                } else {
                    return this->fBkgFunc(x, pars);
                }
            },
            this->fFitRangeMin, fFitRangeMax, fNBkgPars);
    }

    /*
    Set the fit parameters. Signal parameters first, then background.
    format: {index_of_parameter, {init_value, lower_lim, upper_lim}}

    Parameters:
        - name: unique identifier of the set of parameters. Format: pdg_signal-func_bkg-func_suffix
    */
    void SetFitSettings(std::string name, int npart = 0, int version = 0) {
        if (name == "421_gaus_exp_alice3_siblings") {
            fFitPars = {
                // gaussian
                {0, {"norm", 1000, 0, 1.e6}},
                {1, {"mean", 1.86484, 1.86484*0.98, 1.86484*1.02}},
                {2, {"sigma", 10, 0.005, 0.05}},

                // exp
                {3, {"norm", 10, 0.01, 1.e4}},
                {4, {"slope", -0.5, -10, 10}},
            };
        } else if (name == "431_gaus_powex_DstarFemto") {
            fFitPars = {
                // gaussian
                {0, {"norm", 0.5, 0, 50}},
                {1, {"mean", 0.145, 0.144, 0.146}},
                {2, {"sigma", 0.0006, 0.0003, 0.0009}},

                // powex
                {3, {"norm", 0.5, 200, 3000}},
                {4, {"slope", 0.1, 0, 100}},
            };
        } else if (name == "411_gaus_exp_DstarFemto") {
            fFitPars = { {0, {"norm", 0.5, 0, 50}},
                {0, {"mean", 1.87, 1.84, 1.9}},
                {2, {"sigma", 0.006, 0.0001, 0.01}},

                // exp
                {3, {"norm", 10, 0, 1.e6}},
                {4, {"slope", -0.5, -10, 0}},
            };
        } else if (name == "411_hat_exp_DstarFemto") {
            fFitPars = {    // gauss1
                {0, {"norm", 0.1, 0, 50}},
                {1, {"mean", 0.145, 0.144, 0.146}},
                {2, {"sigma", 0.001, 0.0002, 0.002}},

                // gauss2
                {3, {"Yfrac", 0.5, 0.1, 1}},
                {4, {"sigmaFrac", 1.5, 1, 5}},

                // exp
                {5, {"norm", 1, 0, 1e6}},
                {6, {"slope", 0.1, -100, 100}},
            };
        } else if ("Purity") {
            const YAML::Node& parSett = fCfgFile["neutrpart"][npart]["fitsettings"][version];
            if(this->fBkgFuncName == "spline5"){
                for (int iPar = 0; iPar < this->fNSgnPars; iPar++){
                    fFitPars.insert({iPar, {parSett["p" + std::to_string(iPar)][0].as<std::string>(),
                                        parSett["p" + std::to_string(iPar)][1].as<double>(),
                                        parSett["p" + std::to_string(iPar)][2].as<double>(),
                                        parSett["p" + std::to_string(iPar)][3].as<double>()}});
                }
                for (int iPar = 0; iPar < this->fNBkgPars/2; iPar++){
                    double xKnotFix = parSett["xknots"][iPar].as<double>();
                    double yKnotSet = this->fHist->GetBinContent(this->fHist->FindBin(xKnotFix));
                    fFitPars.insert({iPar + this->fNSgnPars, {"pXKn" + std::to_string(iPar), 
                                     xKnotFix, xKnotFix, xKnotFix}});
                    fFitPars.insert({iPar + this->fNSgnPars + this->fNBkgPars/2, {"pYKn" + std::to_string(iPar), 
                                     yKnotSet, yKnotSet - 1e3, yKnotSet + 1e3}});
                }
            }
            else {
                for (int iPar = 0; iPar < this->fNSgnPars + this->fNBkgPars; iPar++){
                    fFitPars.insert({iPar, {parSett["p" + std::to_string(iPar)][0].as<std::string>(),
                                            parSett["p" + std::to_string(iPar)][1].as<double>(),
                                            parSett["p" + std::to_string(iPar)][2].as<double>(),
                                            parSett["p" + std::to_string(iPar)][3].as<double>()}});
                }
            }
        } else {
            std::cout << "The set of parameters '" << name << "' is not valid. Exit!" << std::endl;
            exit(1);
        }
    }

    int Fit() {
        // printf("---> %s\n", this->bkg.data());
        // prefit
        // this->fPrefit = new TF1("fPrefit", [&, this](double *x, double * par) {
        //     if (std::abs(x[0] - 0.145) < 0.001)
        //         TF1::RejectPoint();
        //         return 0;
        //     return this->fBkgFunc.;EvalPar(x, par);cc

        // }, this->fFitRangeMin, fFitRangeMax, fNBkgPars);
        // printf("\n\n\nPerfomring the prefit to the background:\n");
        // hist->Fit(this->fPrefit, "QMR0+", "");

        for (int iPar = 0; iPar < this->fNSgnPars + this->fNBkgPars; iPar++) {
            auto pars = fFitPars[iPar];
            std::cout << iPar << " " 
                      << std::get<0>(pars)
                      << " " << std::get<1>(pars)
                      << " " << std::get<2>(pars)
                      << " " << std::get<3>(pars)
                      << std::endl;

            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            this->fFit->SetParameter(iPar, std::get<1>(pars));
            this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
        }

        // printf("\n\n\nPerfomring the full fit:\n");
        int status = fHist->Fit(this->fFit, "SMRL+0", "")->Status();

        // decompose the fit function in its contributions
        this->fBkg = new TF1(fBkgFuncName.data(), fBkgFunc, fFitRangeMin, fFitRangeMax, fNBkgPars);
        for(int iBkgPar = 0; iBkgPar < fNBkgPars; iBkgPar++){
            this->fBkg->SetParameter(iBkgPar, this->fFit->GetParameter(this->fNSgnPars + iBkgPar));
        }

        if (this->fSgnFuncName == "hat") {
            this->fHatThin = new TF1("fHatThin", Gaus, fFitRangeMin, fFitRangeMax, 3);
            this->fHatThin->SetParameter(0, this->fFit->GetParameter(0) * this->fFit->GetParameter(3));
            this->fHatThin->SetParameter(1, this->fFit->GetParameter(1));
            this->fHatThin->SetParameter(2, this->fFit->GetParameter(2));

            this->fHatWide = new TF1("fHatThin", Gaus, fFitRangeMin, fFitRangeMax, 3);
            this->fHatWide->SetParameter(0, this->fFit->GetParameter(0) * (1 - this->fFit->GetParameter(3)));
            this->fHatWide->SetParameter(1, this->fFit->GetParameter(1));
            this->fHatWide->SetParameter(2, this->fFit->GetParameter(2) * this->fFit->GetParameter(4));

            this->fSgn = new TF1("fSgn", Hat, fFitRangeMin, fFitRangeMax, 5);
            this->fSgn->SetParameter(0, this->fFit->GetParameter(0));
            this->fSgn->SetParameter(1, this->fFit->GetParameter(1));
            this->fSgn->SetParameter(2, this->fFit->GetParameter(2));
            this->fSgn->SetParameter(3, this->fFit->GetParameter(3));
            this->fSgn->SetParameter(4, this->fFit->GetParameter(4));

        } else {
            this->fSgn = new TF1(fSgnFuncName.data(), fSgnFunc, fFitRangeMin, fFitRangeMax, fNSgnPars);
            for(int iSgnPar = 0; iSgnPar < fNSgnPars; iSgnPar++){
                this->fSgn->SetParameter(iSgnPar, this->fFit->GetParameter(iSgnPar));
            }
        }
        return status;
    }

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad, std::string method) {
        pad->cd();
        fHist->GetYaxis()->SetRangeUser(0, 1.3 * fHist->GetMaximum());
        gPad->DrawFrame(fFitRangeMin, 0, fFitRangeMax, 1.3 * fHist->GetMaximum(),
                        Form("%s;%s;%s", this->fHist->GetTitle(), this->fHist->GetXaxis()->GetTitle(),
                             this->fHist->GetYaxis()->GetTitle()));
        if (this->fBkgFuncName == "pol1") {
            this->fBkg->SetNpx(300);
            this->fBkg->SetLineColor(kGray + 2);
            this->fBkg->Draw("same");
        } else if (this->fBkgFuncName == "pol5") {
            this->fBkg->SetNpx(300);
            this->fBkg->SetLineColor(kGray + 2);
            this->fBkg->Draw("same");
        } else if (this->fBkgFuncName == "powex" || this->fBkgFuncName == "exp") {
            this->fBkg->SetNpx(300);
            this->fBkg->SetLineColor(kGray + 2);
            this->fBkg->Draw("same");
        } 
        else if (this->fBkgFuncName == "spline5") {
            this->fBkg->SetNpx(300);
            this->fBkg->SetLineColor(kGray + 2);
            this->fBkg->Draw("same");
        } 

        if (this->fSgnFuncName == "hat") {
            this->fHatThin->SetLineColor(kMagenta + 3);
            this->fHatThin->SetNpx(300);
            this->fHatThin->Draw("same");

            this->fHatWide->SetNpx(300);
            this->fHatWide->SetLineColor(kAzure + 2);
            this->fHatWide->Draw("same");

            this->fSgn->SetNpx(300);
            this->fSgn->SetLineColor(kBlue + 2);
            this->fSgn->Draw("same");
        } else if (this->fSgnFuncName == "gaus") {
            this->fSgn->SetNpx(300);
            this->fSgn->SetLineColor(kBlue + 2);
            this->fSgn->Draw("same");
        }

        fPrefit->SetLineStyle(9);
        fPrefit->SetLineColor(kGray + 2);
        fPrefit->Draw("same");

        fFit->SetLineColor(kRed);
        fFit->Draw("same");

        fHist->SetMarkerSize(1);
        fHist->SetMarkerStyle(20);
        fHist->SetMarkerColor(kBlack);
        fHist->SetLineColor(kBlack);
        fHist->SetLineWidth(2);
        fHist->Draw("same pe");

        TLatex tl;
        tl.SetTextSize(0.035);
        tl.SetTextFont(42);
        double nSigma = 2;
        double step = 0.05;
        int iStep = 0;
        tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("#chi^{2}/NDF = %.2f", fFit->GetChisquare() / fFit->GetNDF()));
        tl.DrawLatexNDC(.15, .85 - step * iStep++,
                        Form("S(%.2f#sigma) = %.2f #pm %.2f", nSigma, this->GetSignal(nSigma, method),
                             this->GetSignalUnc(nSigma, method)));

        tl.DrawLatexNDC(
            .15, .85 - step * iStep++,
            Form("B(%.2f#sigma) = %.2f #pm %.2f", nSigma, this->GetBackground(nSigma), this->GetBackgroundUnc(nSigma)));
        
        tl.DrawLatexNDC(.15, .85 - step * iStep++, Form("Counts = %.2f", this->GetCounts()));
        pad->Update();
    }

    double GetMean() {
        if (!fFit) return -1;
        if (this->fSgnFuncName == "gaus" || this->fSgnFuncName == "hat") return fFit->GetParameter(1);
        return -1;
    }

    double GetMeanUnc() {
        if (!fFit) return -1;
        if (this->fSgnFuncName == "gaus" || this->fSgnFuncName == "hat") return fFit->GetParError(1);
        return -1;
    }

    double GetWidth() {
        if (!fFit) return -1;

        double probs[2] = {TMath::Freq(-1), TMath::Freq(1)};
        double quantiles[2];
        fSgn->GetQuantiles(2, quantiles, probs);
        return (quantiles[1] - quantiles[0]) / 2;
    }

    double GetWidthUnc() {
        if (!fFit) return -1;

        if (this->fSgnFuncName == "gaus")
            return fFit->GetParError(2);
        else if (this->fSgnFuncName == "hat")
            return fFit->GetParError(2) / fFit->GetParameter(2) * GetWidth();  // assume rel error is the same
        return -1;
    }

    double GetCounts() {
        if (!fFit) return -1;
        int firstBin = fHist->GetXaxis()->FindBin(fFitRangeMin * 1.0001);
        int lastBin = fHist->GetXaxis()->FindBin(fFitRangeMax * 0.9999);

        double totCounts = fHist->Integral(firstBin, lastBin);
        double bkgCounts =
            fBkg->Integral(fHist->GetXaxis()->GetBinLowEdge(firstBin), fHist->GetXaxis()->GetBinLowEdge(lastBin + 1)) /
            this->fHist->GetBinWidth(1);
        return totCounts - bkgCounts;
    }

    double GetSignal(double nSigma, std::string method) {
        if (!fFit) return -1;

        // convert nSigma to quantiles
        double probs[2] = {TMath::Freq(-nSigma), TMath::Freq(nSigma)};
        double quantiles[2];
        fSgn->GetQuantiles(2, quantiles, probs);  // returns the limits in which the sgn func should be integrated

        if (method == "sgn_int") {
            return fSgn->Integral(quantiles[0], quantiles[1]) / this->fHist->GetBinWidth(1);
        } else if (method == "data_minus_bkg") {
            int firstBin = fHist->GetXaxis()->FindBin(quantiles[0] * 1.0001);
            int lastBin = fHist->GetXaxis()->FindBin(quantiles[1] * 0.9999);

            double data = this->fHist->Integral(firstBin, lastBin);
            double bkg = GetBackground(nSigma);
            return data - bkg;
        } else {
            printf("'%s' is an invalid signal extraction method. Exit!\n", method.data());
            exit(1);
        }
        return -1;
    }
    double GetChi2Ndf() { return fFit->GetChisquare() / fFit->GetNDF(); }

    double GetBackground(double nSigma) {
        if (!fFit) return -1;
        double probs[2] = {TMath::Freq(-nSigma), TMath::Freq(nSigma)};
        double quantiles[2];
        fSgn->GetQuantiles(2, quantiles, probs);  // returns the limits in which the sgn func should be integrated

        return fBkg->Integral(quantiles[0], quantiles[1]) / this->fHist->GetBinWidth(1);
    }

    double GetBackgroundUnc(double nSigma) {
        Int_t leftBand = this->fHist->FindBin(this->GetMean() - 6 * this->GetWidth());
        Int_t rightBand = this->fHist->FindBin(this->GetMean() + 6 * this->GetWidth());

        int start = this->fHist->FindBin(this->fFitRangeMin * 1.0001);
        int end = this->fHist->FindBin(this->fFitRangeMax * 0.9999);
        double SidebandBkg =
            this->fHist->Integral(1, leftBand) + this->fHist->Integral(rightBand, this->fHist->GetNbinsX());

        double sum2 = 0;
        for (Int_t i = 1; i <= leftBand; i++) {
            sum2 += this->fHist->GetBinError(i) * this->fHist->GetBinError(i);
        }
        for (Int_t i = rightBand; i <= this->fHist->GetNbinsX(); i++) {
            sum2 += this->fHist->GetBinError(i) * this->fHist->GetBinError(i);
        }

        return TMath::Sqrt(sum2) / SidebandBkg * this->GetBackground(nSigma);
    }

    double GetSignalUnc(double nSigma, std::string method) {
        if (!fFit) return -1;
        if (this->GetSignal(nSigma, method) <= 0) return 0;

        return fFit->GetParError(0) / fFit->GetParameter(0) * this->GetSignal(nSigma, method);
    }

    TF1* GetSgnFunc() {
        //if (!fFit) return -1;
        return fSgn;
    }

    TF1* GetBkgFunc() {
        //if (!fFit) return -1;
        return fBkg;
    }


 private:
    YAML::Node fCfgFile;

    TH1 *fHist = nullptr;
    TF1 *fSgn = nullptr;
    TF1 *fHatThin = nullptr;
    TF1 *fHatWide = nullptr;
    TF1 *fBkg = nullptr;
    TF1 *fFit = nullptr;
    TF1 *fPrefit = nullptr;

    std::string fSgnFuncName;
    std::string fBkgFuncName;
    const char* fConfigPath;

    double (*fSgnFunc)(double *x, double *par);
    double (*fBkgFunc)(double *x, double *par);

    int fNSgnPars;
    int fNBkgPars;
    double fFitRangeMin;
    double fFitRangeMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;
};

#endif  // FEMPY_MASSFITTER_HXX_
