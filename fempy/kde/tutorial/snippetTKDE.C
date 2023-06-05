#include "Math/DistFunc.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TKDE.h"
#include "TLegend.h"
#include "TRandom.h"
// #include "TSpline3.h"


double ratioSplines(TSpline3 *s1, TSpline3 *s2, double x){
    return s1->Eval(x) / s2->Eval(x);
}

TF1 *fKde1; 
TF1 *fKde2; 

double ratioFunc(TF1 *f1, TF1 *f2, double x){
    return f1->Eval(x) / f2->Eval(x);
}


double ratioFuncRoot(double *x, double *par){
    return ratioFunc(fKde1, fKde2, x[0]);
}


void snippetTKDE() {
    int n = 1000;
    int nbin = 100;
    double xmin = 0;
    double xmax = 10;
    std::vector<double> data1(n);
    std::vector<double> data2(n);

    for (int i = 0; i < n; ++i) {
        data1[i] = gRandom->Gaus(1, 1);
        data2[i] = gRandom->Gaus(2, 1);
    }
    TKDE *kde1 = new TKDE(n, &data1[0], xmin, xmax);
    TKDE *kde2 = new TKDE(n, &data2[0], xmin, xmax);
    
    fKde1 = kde1->GetFunction(100, 0, 3);
    fKde2 = kde2->GetFunction(100, 0, 3);

    // fKde1->SetName("f1");
    // fKde2->SetName("f2");
    

    // TSpline3 *spl1 = new TSpline3("spl1", 0, 3, fKde1, 1000, 0, fKde1->Eval(0), fKde1->Eval(3));
    // TSpline3 *spl2 = new TSpline3("spl2", 0, 3, fKde2, 1000, 0, fKde2->Eval(0), fKde2->Eval(3));

    // TF1 * cf = TF1("ff", "[&]")
    // TF1 * cf  = new TF1("cf1", ,0,10,0);

    // auto mylambda = ratioSplines
    // auto mylambda = [&](double *x, double *p){ return spl1->Eval(x[0]) + spl2->Eval(x[0]); };
    // TF1 * cf  = new TF1("cf1",[&](double *x, double *p){ return spl1->Eval(x[0]) + spl2->Eval(x[0]); },0,10,0);
    TF1 * cf  = new TF1("cf1",ratioFuncRoot,0,3);
    
    cf->Draw();
    // // TF1 *f1f2 = new TF1("f1f2", "f1 / f2", xmin, xmax);


    // fKde1->Draw("same");
    // spl1->Draw("same");




    

    // std::cout<<fKde1->GetFormula()<<std::endl;
    // TF1 *ratio = new TF1("ratio", Form("%s/%s", fKde1->GetName(), fKde2->GetName()), 0, 10);
    // ratio->Draw("");

    // /////

    // TF1 *g1 = new TF1("g1", "x*x+1", 0, 10);
    // TF1 *g2 = new TF1("g2", "x*x+2", 0, 10);

    // std::cout<<g1->GetFormula()<<std::endl;
    // // TF1 *g1g2 = new TF1("g1g2", "g1/g2", 0, 10);
    // // g1g2->Draw();

}
