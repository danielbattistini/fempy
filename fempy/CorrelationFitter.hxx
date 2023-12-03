#ifndef FEMPY_CORRELATIONFITTER_HXX_
#define FEMPY_CORRELATIONFITTER_HXX_

#include <map>
#include <string>
#include <tuple>

#include "TH1.h"
#include "TFitResult.h" 
#include "TF1.h"

double Pol0(double *x, double *par) { return par[0]; }

class CorrelationFitter {
   public:
    CorrelationFitter(TH1 *hist, double fitRangeMin, double fitRangeMax) {
        this->hist = reinterpret_cast<TH1 *>(hist->Clone());
        this->fitRangeMin = fitRangeMin;
        this->fitRangeMax = fitRangeMax;
    }

    void Add(std::string name, std::vector<std::tuple<std::string, double, double, double>> pars) {
        if (name == "pol0") {  // Constant function
            this->fitFunc.push_back(Pol0);
        } else {
            printf("Error: function '%s' is not implemented. Exit!", name.data());
            exit(1);
        }

        // Save fit settings
        for (const auto &par : pars) {
            this->fFitPars.insert({this->fFitPars.size(), par});
        }
    }

    int Fit() {
        // Build the fit function
        this->fFit = new TF1(
            "fFit",
            [&, this](double *x, double *pars) {
                double result = 0;
                for (const auto &func : this->fitFunc) {
                    result += func(x, pars);
                }
                return result;
            },
            fitRangeMin, fitRangeMax, this->fFitPars.size());

        for (size_t iPar = 0; iPar < this->fFitPars.size(); iPar++) {
            auto pars = this->fFitPars[iPar];
            std::cout << iPar << " " << std::get<0>(pars) << " " << std::get<1>(pars) << " " << std::get<2>(pars) << " "
                      << std::get<3>(pars) << std::endl;

            this->fFit->SetParName(iPar, std::get<0>(pars).data());
            this->fFit->SetParameter(iPar, std::get<1>(pars));
            this->fFit->SetParLimits(iPar, std::get<2>(pars), std::get<3>(pars));
        }

        int status = hist->Fit(this->fFit, "SMR+0", "")->Status();

        return status;
    }

   private:
    TH1 *hist = nullptr;
    TF1 *fFit = nullptr;

    std::vector<double (*)(double *x, double *par)> fitFunc;

    double fitRangeMin;
    double fitRangeMax;
    std::map<int, std::tuple<std::string, double, double, double>> fFitPars;
};

#endif  // FEMPY_CORRELATIONFITTER_HXX_
