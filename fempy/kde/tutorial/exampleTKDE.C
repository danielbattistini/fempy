/// \file
/// \ingroup tutorial_math
/// \notebook
/// Example of using the TKDE class (kernel density estimator)
///
/// \macro_image
/// \macro_code
///
/// \authors Lorenzo Moneta, Bartolomeu Rabacal (Dec 2010)

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TKDE.h"
/*#include "TStopwatch.h"*/
#include "Math/DistFunc.h"
#include "TF1.h"
#include "TLegend.h"
#include "TRandom.h"

// test TKDE

// struct GlobalChi2 {
//     GlobalChi2(TF1 *f1, TF1 *f2) : fTheo(f1), fExp(f2) {}

//     // parameter vector is first background (in common 1 and 2)
//     // and then is signal (only in 2)
//     double operator()(const double *par) const {
//         double p1[2];
//         for (int i = 0; i < 2; ++i) p1[i] = par[iparB[i]];

//         double p2[5];
//         for (int i = 0; i < 5; ++i) p2[i] = par[iparSB[i]];

//        return (fTheo)(p1) + (*fExp)(p2);
//     }
//     // const?
//     const TF1 *fTheo;
//     const TF1 *fExp;
// };

double mychi2impl(double x) { return x * x - x; }
struct customChi2 {
    customChi2(TF1 *f1, TF1 *f2) : fTheo(f1), fExp(f2) {}

    // parameter vector is first background (in common 1 and 2)
    // and then is signal (only in 2)
    double operator()(const double *par) const {
        fTheo->SetParameter(0, par[0]);

        return fTheo->Integral(0, 10) - fExp->Integral(0, 10);
    }

    TF1 *fTheo;
    TF1 *fExp;
};

void exampleTKDE(int n = 1000) {
    // generate some gaussian points

    int nbin = 100;
    double xmin = 0;
    double xmax = 10;

    TH1D *h1 = new TH1D("h1", "h1", nbin, xmin, xmax);

    // generate some points with bi- gaussian distribution

    std::vector<double> data(n);
    for (int i = 0; i < n; ++i) {
        if (i < 0.4 * n) {
            data[i] = gRandom->Gaus(2, 1);
            h1->Fill(data[i]);
        } else {
            data[i] = gRandom->Gaus(7, 1.5);
            h1->Fill(data[i]);
        }
    }

    // scale histogram
    h1->Scale(1. / h1->Integral(), "width");
    h1->SetStats(false);
    h1->SetTitle("Bi-Gaussian");
    h1->Draw();

    // drawn true normalized density
    TF1 *f1 = new TF1("f1", "0.4*ROOT::Math::normal_pdf(x,1,2)+0.6*ROOT::Math::normal_pdf(x,1.5,7)", xmin, xmax);
    TF1 *f1Fit =
        new TF1("f1Fit", "[0]*ROOT::Math::normal_pdf(x,1,2)+(1-[0])*ROOT::Math::normal_pdf(x,1.5,7)", xmin, xmax);
    f1->SetLineColor(kGreen + 2);
    f1->Draw("SAME");

    // create TKDE class
    double rho = 1.0;  // default value
    TKDE *kde = new TKDE(n, &data[0], xmin, xmax, "", rho);
    // kde->Draw("ConfidenceInterval@95 Same");
    // kde->Draw("SAME");
    // kde->GetDrawnLowerFunction()->Draw("same");
    // kde->GetDrawnUpperFunction()->Draw("same");

    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.95);
    legend->AddEntry(f1, "True function");
    // legend->AddEntry(kde->GetDrawnFunction(),"TKDE");

    // kde->GetGraphWithErrors(100, 0, 10)->Draw("same");
    // kde->GetUpperFunction(0.68, 100, 0, 10)->Draw("same");
    // legend->AddEntry(kde->GetDrawnLowerFunction(),"TKDE - #sigma");
    // legend->AddEntry(kde->GetDrawnUpperFunction(),"TKDE + #sigma");
    // legend->Draw();

    TF1 *fKde = kde->GetFunction(100, 0, 10);
    fKde->AddToGlobalList();
    fKde->SetName("ciccio");
    TF1 *f2 = new TF1("f", "ciccio", 0, 10);
    f2->SetLineColor(kMagenta);
    f2->Draw("same");
    // fit kde
    // ROOT::Fit::Fitter fitter;
    // customChi2 myChi2(f1Fit, fKde);
    // fitter.FitFCN(1, mychi2impl, 0, n, true);
    // fitter.FitFCN(1, myChi2, 0, n, true);
    // fitter.Config().SetParamsSettings(0, 0.5);

    // fitter.Result();

    // fKde->Draw("same");
    // legend->AddEntry(fKde,"TKDE estim");
}
