#ifndef FEMPY_CORRELATIONFITTER_HXX_
#define FEMPY_CORRELATIONFITTER_HXX_

#include <map>
#include <string>
#include <tuple>

#include "TF1.h"
#include "TFitResult.h"
#include "TH1.h"
#include "TMath.h"

double Pol0(double *x, double *par) { return par[0]; }

double BreitWigner(double *x, double *par) {
    double kstar = x[0];

    double yield = par[0];
    double mean = par[1];
    double gamma = par[2];

    return yield * gamma / TMath::Pi() / (gamma * gamma + (kstar - mean) * (kstar - mean));
}

class CorrelationFitter {
   public:
    CorrelationFitter(TH1 *hist, double fitRangeMin, double fitRangeMax) {
        this->hist = reinterpret_cast<TH1 *>(hist->Clone());
        this->fitRangeMin = fitRangeMin;
        this->fitRangeMax = fitRangeMax;
        this->nPars = {0};  // The first parameter has index zero
    }

    /*
    Add a term to the CF model. Available options:
        - pol0
        - bw (Breit-Wigner)

    In pyroot, the fit parameters (pars) are specified as a list of tuples.

    The syntax is: `(name, initialValue, lowerLim, upperLim)`. To fix a parameter set `lowerLim > upperLim`, e.g
    `fitter.Add('bw',
    [('p0', 0.9, 1, -1)])`

    Example:
    ```fitter.Add('bw', [
        ('yield', 0.03, 0.001, 0.1),
        ('mean', 0.126, 0.125, 0.127),
        ('gamma', 0.004, 0.003, 0.005),
    ])```

    It is your responsibility to give the correct number of parameters. All fit parameters must be initialized.
    */
    void Add(std::string name, std::vector<std::tuple<std::string, double, double, double>> pars) {
        if (name == "pol0") {  // Constant function
            this->fitFunc.push_back(Pol0);
        } else if (name == "bw") {  // Breit-Wigner function
            this->fitFunc.push_back(BreitWigner);
        } else {
            printf("Error: function '%s' is not implemented. Exit!", name.data());
            exit(1);
        }

        this->nPars.push_back(pars.size());

        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    // Perform the fit
    int Fit() {
        // Build the fit function
        this->fFit = new TF1(
            "fFit",
            [&, this](double *x, double *pars) {
                double result = 0;
                for (int iTerm = 0; iTerm < this->fitFunc.size(); iTerm++) {
                    auto func = this->fitFunc[iTerm];
                    result += func(x, &pars[this->nPars[iTerm]]);  // Shift the index of the parameters
                }
                return result;
            },
            fitRangeMin, fitRangeMax, this->fFitPars.size());

        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];

            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            if (std::get<2>(pars) > std::get<3>(pars)) {
                this->fFit->FixParameter(iPar, std::get<1>(pars));
            } else {
                this->fFit->SetParameter(iPar, std::get<1>(pars));
                this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
            }
        }

        int status = hist->Fit(this->fFit, "SMR+0", "")->Status();

        return status;
    }

    /*
    Define a canvas before calling this function and pass gPad as TVirtualPad
    */
    void Draw(TVirtualPad *pad) {
        pad->cd();
        hist->GetYaxis()->SetRangeUser(0, 1.3 * hist->GetMaximum());
        gPad->DrawFrame(fitRangeMin, 0, fitRangeMax, 1.3 * hist->GetMaximum(),
                        Form("%s;%s;%s", this->hist->GetTitle(), this->hist->GetXaxis()->GetTitle(),
                             this->hist->GetYaxis()->GetTitle()));
        this->fFit->SetNpx(300);
        this->fFit->SetLineColor(kGray + 2);
        this->fFit->SetLineColor(kRed);
        this->fFit->Draw("same");

        hist->SetMarkerSize(1);
        hist->SetMarkerStyle(20);
        hist->SetMarkerColor(kBlack);
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(2);
        hist->Draw("same pe");

        pad->Update();
    }

   private:
    TH1 *hist = nullptr;
    TF1 *fFit = nullptr;

    std::vector<double (*)(double *x, double *par)> fitFunc;  // List of function describing each term of the CF model
    std::vector<int> nPars;                                   // Keeps track of how many parameters each function has

    double fitRangeMin;
    double fitRangeMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;  // List of fit parameters
};

#endif  // FEMPY_CORRELATIONFITTER_HXX_
