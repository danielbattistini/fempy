#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TCanvas.h"
#include "TError.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom2.h"
#include <iostream>

// using namespace std::placeholders;  // for _1, _2, _3...

class FitFunctionToFunction : public TObject {
  public:
    FitFunctionToFunction(){};
    FitFunctionToFunction(TF1 *f1, TF1 *f2) : fFunc1(f1), fFunc2(f2) {
        printf("constr\n");
        fFunc1->SetName("fFunc1");
        fFunc2->SetName("fFunc2");
    };
    ~FitFunctionToFunction(){};

    double customChi2(const double *fitpar) {
        const Double_t p0 = fitpar[0];
        const Double_t p1 = fitpar[1];

        fFunc2->SetParameter(0, p0);
        fFunc2->SetParameter(1, p1);
        
        TF1 *fDelta = new TF1("f", "(fFunc1-fFunc2)**2", 0, 1);
        
        return fDelta->Integral(0, 1);
    };
    int Fit() {
        ROOT::Math::Minimizer *minimum = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
        minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
        minimum->SetMaxIterations(10000);      // for GSL
        minimum->SetTolerance(0.001);
        minimum->SetPrintLevel(2);

        // typedef  double (FitFunctionToFunction::*F2FMem)(const double *);
        // F2FMem chisq = &FitFunctionToFunction::customChi2;

        auto ptr_to_print_sum = std::mem_fn(&FitFunctionToFunction::customChi2);
        auto f4 = std::bind(ptr_to_print_sum, this, std::placeholders::_1);
        // auto chisq = [&](double *fitpar) { return customChi2(fitpar); };
        // auto chisq = [&]
        ROOT::Math::Functor f(f4, 2);
        double step[2] = {0.01, 0.01};
        double variable[2] = {0.5, 0.5};

        minimum->SetFunction(f);
        minimum->SetVariable(0, "p0", variable[0], step[0]);
        minimum->SetVariable(1, "p1", variable[1], step[1]);

        // do the minimization
        minimum->Minimize();
        const double *xs = minimum->X();

        // minimum->Dump();

        return 0;
    };

  private:
    TF1 *fFunc1;
    TF1 *fFunc2;
};

double 
void test_class() {
    TF1 *f1 = new TF1("f1", "x+1", 0, 10);
    TF1 *f2 = new TF1("f2", "[0]*x*x + [1]", 0, 10);

    FitFunctionToFunction fitter(f1, f2);
    fitter.Fit();
}