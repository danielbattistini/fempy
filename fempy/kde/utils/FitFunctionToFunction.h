#ifndef FEMPY_KDE_UTILS_FITFUNCTIONTOFUNCTION_H_
#define FEMPY_KDE_UTILS_FITFUNCTIONTOFUNCTION_H_

#include <iostream>
#include <string>
#include <map>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom2.h"

// using namespace std::placeholders; // for _1, _2, _3...
struct variableSettings {
    std::string fName;
    Double_t fStart;
    Double_t fLimLow;
    Double_t fLimUp;
    Double_t fStep;
};

class FitFunctionToFunction : public TObject {
 public:
    FitFunctionToFunction() {}
    FitFunctionToFunction(TF1 *f1, TF1 *f2, TF1 *f2UpperUnc, TF1 *f2LowerUnc, Int_t nPar)
        : fFunc1(f1), fFunc2(f2), fNPar(nPar), fFunc1UpperUnc(f2UpperUnc), fFunc2LowerUnc(f2LowerUnc) {
        // fFunc1->SetName("fFunc1");
        // fFunc2->SetName("fFunc2");
        // fFunc2UpperUnc->SetName("fFunc2UpperUnc");
        // fFunc2LowerUnc->SetName("fFunc2LowerUnc");
        this->fFunc1->AddToGlobalList();
        this->fFunc2->AddToGlobalList();
        this->fFunc1UpperUnc->AddToGlobalList();
        this->fFunc2LowerUnc->AddToGlobalList();
    }
    ~FitFunctionToFunction() {}

    double customChi2(const double *fitpar);
    double GetChi2() const { return this->fChi2; }
    int Fit();

    void SetFitRange(Double_t xMin, Double_t xMax) {
        this->fRangeLow = xMin;
        this->fRangeUp = xMax;
    }
    void SetVariableSettings(Int_t iVar, variableSettings settings) { this->fVarSettings.insert({iVar, settings}); }

    // double customChi2Formula(double *xx, double *par) {
    // const Double_t p0 = xx[0];
    // const Double_t p1 = xx[1];

    // fTheo->SetParameter(0, p0);
    // fTheo->SetParameter(1, p1);
    // // std::cout<<"p0: " << p0 <<"p1: " << p1<< "    int theo: "<<fExp->Integral(0, 1)<< "    int exp:
    // "<<fTheo->Integral(0, 1)<<std::endl; double chi2 = (fExp->Integral(0, 1) - fTheo->Integral(0,
    // 1))*(fExp->Integral(0, 1) - fTheo->Integral(0, 1));

    //     return customChi2(xx);
    // }

    // void printChi2ToFile() {
    //     TCanvas *cc = new TCanvas("c2", "c", 600, 600);

    //     auto ptr_to_print_sum = std::mem_fn(&FitFunctionToFunction::customChi2Formula);
    //     auto f4 = std::bind(ptr_to_print_sum, this, std::placeholders::_1, std::placeholders::_1);

    //     TF2 *fChi2 = new TF2("fChi2", f4, -2, 2, -2, 2);
    //     fChi2->Draw("colz");
    //     cc->SaveAs("plotChi2.pdf");
    // }

 private:
    TF1 *fFunc1;
    TF1 *fFunc2;
    TF1 *fFunc1UpperUnc;
    TF1 *fFunc2LowerUnc;

    Int_t fNPar;
    Double_t fRangeLow;
    Double_t fRangeUp;
    Double_t fChi2;

    ROOT::Math::Minimizer *fMinimizer;

    std::map<Int_t, variableSettings> fVarSettings;
}

#endif  // FEMPY_KDE_UTILS_FITFUNCTIONTOFUNCTION_H_
