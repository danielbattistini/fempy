#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "TError.h"
#include "TF1.h"
#include "TF2.h"
#include "TRandom2.h"
#include "TCanvas.h"
#include "TObject.h"
#include <iostream>

double RosenBrock(const double *xx) {
    const Double_t x = xx[0];
    const Double_t y = xx[1];
    const Double_t tmp1 = y - x * x;
    const Double_t tmp2 = 1 - x;
    return 100 * tmp1 * tmp1 + tmp2 * tmp2;
}
TF1 *fExp = new TF1("fExp", "x + x*x", 0, 1);
TF1 *fTheo = new TF1("fTheo", "[0]+[1]*x", 0, 1);

double customChi2(const double *fitpar) {
    const Double_t p0 = fitpar[0];
    const Double_t p1 = fitpar[1];

    fTheo->SetParameter(0, p0);
    fTheo->SetParameter(1, p1);
    // std::cout<<"p0: " << p0 <<"p1: " << p1<< "    int theo: "<<fExp->Integral(0, 1)<< "    int exp: "<<fTheo->Integral(0, 1)<<std::endl;
    TF1 *fDelta = new TF1("f", "(fExp-fTheo)**2", 0, 1);
    // double diff = (fExp->Integral(0, 1) - fTheo->Integral(0, 1));
    return fDelta->Integral(0, 1);
}

double customChi2Formula(double *xx, double *par) {
    const Double_t p0 = xx[0];
    const Double_t p1 = xx[1];

    // fTheo->SetParameter(0, p0);
    // fTheo->SetParameter(1, p1);
    // // std::cout<<"p0: " << p0 <<"p1: " << p1<< "    int theo: "<<fExp->Integral(0, 1)<< "    int exp: "<<fTheo->Integral(0, 1)<<std::endl;
    // double chi2 = (fExp->Integral(0, 1) - fTheo->Integral(0, 1))*(fExp->Integral(0, 1) - fTheo->Integral(0, 1));
    return customChi2(xx);
}

int testCustomChi2Minimization(char *minName = "Minuit2", const char *algoName = "", int randomSeed = -1) {
    // create minimizer giving a name and a name (optionally) for the specific
    // algorithm
    // possible choices are:
    //     minName                  algoName
    // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
    //  Minuit2                     Fumili2
    //  Fumili
    //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
    //                              BFGS2, SteepestDescent
    //  GSLMultiFit
    //   GSLSimAn
    //   Genetic
    // ROOT::Math::Minimizer *minimum = ROOT::Math::Factory::CreateMinimizer(minName, algoName);
    ROOT::Math::Minimizer *minimum = ROOT::Math::Factory::CreateMinimizer("Minuit", algoName);
    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(10000);      // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(2);
    // create function wrapper for minimizer
    // a IMultiGenFunction type
    // ROOT::Math::Functor f(&RosenBrock, 2);
    // ROOT::Math::Functor f(&customChi2, 1);
    ROOT::Math::Functor f(&customChi2, 2);
    // double step[1] = {0.01};
    double step[2] = {0.01, 0.01};
    // starting point
    // double variable[1] = {0.5};
    double variable[2] = {0.5, 0.5};

    minimum->SetFunction(f);
    // Set the free variables to be minimized !
    minimum->SetVariable(0, "p0", variable[0], step[0]);
    minimum->SetVariable(1, "p1", variable[1], step[1]);
    // do the minimization
    minimum->Minimize();
    const double *xs = minimum->X();
    // std::cout << "fit p0: " << xs[0]<< "\nfit p1: " << xs[1]  << std::endl;
    // // expected minimum is 0
    // if (minimum->MinValue() < 1.E-4 && f(xs) < 1.E-4)
    //     std::cout << "Minimizer " << minName << " - " << algoName << "   converged to the right minimum" << std::endl;
    // else {
    //     std::cout << "Minimizer " << minName << " - " << algoName << "   failed to converge !!!" << std::endl;
    //     Error("NumericalMinimization", "fail to converge");
    // }

    TCanvas *c = new TCanvas("c", "c", 600, 600);
    fExp->Draw();
    fTheo->Draw("same");


    TCanvas *cc = new TCanvas("c2", "c", 600, 600);
    cc->SetRightMargin(0.17);
    cc->DrawFrame(-2, -2, 3, 5, ";p0;p1;#chi^{2}");

    TF2 *fChi2 = new TF2("fChi2", customChi2Formula, -2, 3, -2, 5);
    fChi2->GetZaxis()->SetTitle("#chi^{2}");
    fChi2->Draw("same colz");

    cc->SaveAs("correlation_param.pdf");
    return 0;
}