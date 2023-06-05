#include "FitFunctionToFunction.h"

#include <iostream>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom2.h"

// using namespace std::placeholders; // for _1, _2, _3...

double FitFunctionToFunction::customChi2(const double *fitpar) {
    for (auto &[iPar, iSettings] : this->fVarSettings) {
        fFunc2->SetParameter(iPar, fitpar[iPar]);
    }

    TString f1n = this->fFunc1->GetName();
    TString f2n = this->fFunc2->GetName();
    TString f1un = this->fFunc1UpperUnc->GetName();
    TString f2un = this->fFunc2LowerUnc->GetName();

    TString deltaFormula = Form("4/(%f-%f)*((%s-%s)/(%s - %s))**2", this->fRangeUp, this->fRangeLow, f1n.Data(),
                                f2n.Data(), f1un.Data(), f2un.Data());
    TF1 *fDelta = new TF1("f", deltaFormula.Data(), this->fRangeLow, this->fRangeUp);

    // std::cout << "p0: " << p0 << " \tp1: " << p1 << " \tchi2: " <<fDelta->Integral(0, 10)<< std::endl;
    return fDelta->Integral(0, 10);
}

int FitFunctionToFunction::Fit() {
    this->fMinimizer = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
    this->fMinimizer->SetMaxFunctionCalls(1000000);  // for Minuit/Minuit2
    this->fMinimizer->SetMaxIterations(10000);       // for GSL
    this->fMinimizer->SetTolerance(0.001);
    this->fMinimizer->SetPrintLevel(1);

    auto ptr_to_print_sum = std::mem_fn(&FitFunctionToFunction::customChi2);
    auto f4 = std::bind(ptr_to_print_sum, this, std::placeholders::_1);

    ROOT::Math::Functor f(f4, this->fNPar);
    this->fMinimizer->SetFunction(f);

    for (auto &[iPar, iSettings] : this->fVarSettings) {
        std::cout << "iPar: " << iPar << "   fName: " << iSettings.fName.data() << "   fStart: " << iSettings.fStart
                  << "   fLimLow: " << iSettings.fLimLow << "   fLimUp: " << iSettings.fLimUp
                  << "   fStep: " << iSettings.fStep << std::endl;
        this->fMinimizer->SetVariable(iPar, iSettings.fName.data(), iSettings.fStart, iSettings.fStep);
        this->fMinimizer->SetVariableLimits(iPar, iSettings.fLimLow, iSettings.fLimUp);
    }

    // do the minimization
    this->fMinimizer->Minimize();

    this->fChi2 = this->fMinimizer->MinValue();

    return 0;
}
