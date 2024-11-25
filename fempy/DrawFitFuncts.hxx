#ifndef FEMPY_DRAWFITFUNCTS_HXX_
#define FEMPY_DRAWFITFUNCTS_HXX_

#include <map>
#include <string>
#include <tuple>
#include <stdexcept>
#include <numeric>

#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "THashList.h"

#if LOG_LEVEL_DRAW
#ifndef DEBUG
#define DEBUG(msg) std::cout << __FUNCTION__ << "  " << msg << std::endl
#else
#define DEBUG(msg)
#endif
#endif

class DrawFitFuncts {
   public:
    DrawFitFuncts(TH1 *fithist, double drawRangeMin, double drawRangeMax, 
                  bool globNorm=0, int basIdx=-1) {
    
        this->hFitHist = fithist;
        this->fBasIdx = basIdx; 
        this->fMult = false; 
        this->fDrawRangeMin = drawRangeMin;
        this->fDrawRangeMax = drawRangeMax;
        this->fGlobNorm = globNorm; 
    }

    void SetParHist(TH1 *parHist) {
        this->hParameters = parHist;
    }

    void SetTotalFitFunc(TF1 *totalFitFunc) {
        this->fFit = totalFitFunc;
    }

    void SetBasIdx(int basIdx, bool doMultiply) {
        this->fBasIdx = basIdx;
        this->fMult = doMultiply;
    }

    void SetGlobNorm(bool setGlobNorm) {
        this->fGlobNorm = setGlobNorm;
    }

    void AddFitCompName(TString fitFuncComp) {
        cout << "Adding " << fitFuncComp << endl;
        this->fFitFuncNames.push_back(fitFuncComp);
    }

    void AddSplineHisto(TH1 *splinehisto) {
        TH1D *splineHisto = static_cast<TH1D*>(splinehisto);
        TSpline3* spline = new TSpline3(splinehisto);
        DEBUG("Adding histo " << splineHisto->GetName() << endl;
        cout << "Bin content for 3rd bin " << splineHisto->GetBinContent(3) << endl;
        cout << "Spline at 10 MeV/c " << spline->Eval(10));
        this->fSplines.push_back(spline);
    }

    void AddSplineHisto(TGraph *splinegraph) {
        TSpline3* spline = new TSpline3(splinegraph->GetName(), splinegraph);
        DEBUG("Adding TGraph " << splinegraph->GetName() << endl;
        cout << "Spline at 10 MeV/c " << spline->Eval(10));
        this->fSplines.push_back(spline);
    }

    void EvaluateToBeDrawnComponents(std::vector<bool> onBaseline, std::vector<bool> multNorm, std::vector<bool> multGlobNorm, 
                            std::vector<double> funcshifts, int basIdx=-1, std::vector<TString> addComps = {""}) {

        // Warnings that prevent the evaluation from being successful
        if(basIdx == -1){
            std::cerr << "Warning: Baseline is not fixed!" << std::endl;
        }
        if(onBaseline.size() != this->fFitFuncNames.size()){
            std::cerr << "Warning: onbaseline status not defined for all components!" << std::endl;
        }
        if(multNorm.size() != this->fFitFuncNames.size()){
            std::cerr << "Warning: multnorm status not defined for all components!" << std::endl;
        }
        if(multGlobNorm.size() != this->fFitFuncNames.size()){
            std::cerr << "Warning: multglobnorm status not defined for all components!" << std::endl;
        }

        DEBUG("Total number of parameters with subcomponents: " << this->hParameters->GetNbinsX());
        DEBUG("Number of components to be drawn: " << this->fFitFuncNames.size());
        for(int iFunc=0; iFunc<this->fFitFuncNames.size(); iFunc++) {
            DEBUG("Component " << this->fFitFuncNames[iFunc]); 
        }
        // Evaluate the single fit components alone
        int startPar=0;
        int iSpline=0;
        std::vector<TF1 *> rawComps;
        std::vector<int> nParsComps;
        for(int iFunc=0; iFunc<this->fFitFuncNames.size(); iFunc++) {
            std::cout << endl;

            DEBUG("Processing " << this->fFitFuncNames[iFunc]);

            if(this->fFitFuncNames[iFunc].Contains("spline")) { // Build shifted spline
                startPar += 1;
                int iPar = startPar + iFunc + 1;

                auto parName = this->hParameters->GetXaxis()->GetLabels()->At(iPar)->GetName();
                double splineShift = this->hParameters->GetBinContent(iPar + 1);
                DEBUG("Loaded shift parameter from histogram. Bin = " << iPar + 1 << "  " << parName << "=" << splineShift);

                // Build shifted spline
                TF1 *shiftedSpline = new TF1(this->fFitFuncNames[iFunc],
                    [&, this, iSpline] (double *x, double *pars) {
                        return this->fSplines[iSpline]->Eval(x[0] - pars[0]);
                    }, this->fDrawRangeMin, this->fDrawRangeMax, 1);

                nParsComps.push_back(1);
                rawComps.push_back(shiftedSpline);
                rawComps.back()->FixParameter(0, splineShift);
                iSpline++;
            } else {
                nParsComps.push_back(std::get<1>(functions[this->fFitFuncNames[iFunc]]));
                rawComps.push_back(new TF1(this->fFitFuncNames[iFunc], std::get<0>(functions[this->fFitFuncNames[iFunc]]), 
                                           fDrawRangeMin, fDrawRangeMax, std::get<1>(functions[this->fFitFuncNames[iFunc]])));
                DEBUG("Set pars of comp " << iFunc << ", named " << this->fFitFuncNames[iFunc] << ", having " << rawComps.back()->GetNpar() << " parameters" << endl; 
                      cout << "StartPar: " << startPar);
                for(int iPar=0; iPar<rawComps.back()->GetNpar(); iPar++) {
                    DEBUG("Set par n. " << iPar << " to " << this->hParameters->GetXaxis()->GetLabels()->At(startPar+iFunc+iPar+1+this->fGlobNorm)->GetName(); 
                          cout << ", bin no. " << startPar+iFunc+iPar+1 << " bin content: " << this->hParameters->GetBinContent(startPar+iFunc+iPar+2+this->fGlobNorm));
                    rawComps.back()->FixParameter(iPar, this->hParameters->GetBinContent(startPar+iFunc+iPar+2+this->fGlobNorm));
                }
                startPar += std::get<1>(functions[this->fFitFuncNames[iFunc]]);
            }
            DEBUG("Evaluating " << this->fFitFuncNames[iFunc] << " at 200 MeV/c: " << rawComps.back()->Eval(200));
        }
        DEBUG("Number of raw components pre-sum: " << rawComps.size()); 

        // save the normalization constant for which each component has to be multiplied when drawing
        std::vector<double> norms;
        DEBUG("Compute normalization -----------------------------------\\");
        std::cout << std::showpos;
        cout.precision(4);
        std::cout << std::scientific;
        for(int iFunc=0; iFunc<this->fFitFuncNames.size(); iFunc++) {
            if(multNorm[iFunc]) {
                int normIdx = accumulate(nParsComps.begin(), std::next(nParsComps.begin(), iFunc), 0) + iFunc + this->fGlobNorm;
                DEBUG("Set component " + std::to_string(iFunc) + " norm to: ";
                cout << this->hParameters->GetXaxis()->GetLabels()->At(normIdx)->GetName() << ", val: ";
                cout << this->hParameters->GetBinContent(normIdx+1) << ", norm idx: ";
                cout << std::to_string(normIdx));
                norms.push_back(this->hParameters->GetBinContent(normIdx+1));
            } else {
                DEBUG("Set component norm to 1");
                norms.push_back(1.);
            }
        } 
        DEBUG("Compute normalization END --------------------------------/");
        std::cout << std::noshowpos;

        // append to the raw components vector the functions that are sum of more than one component
        DEBUG("Number of raw components: " << rawComps.size());
        if(addComps[0] != "") {
            for(int iAddComp=0; iAddComp<addComps.size(); iAddComp++) {
                DEBUG("Evaluating sum of components " << addComps[iAddComp]);
                DEBUG("Function name " << "SumComp_" + addComps[iAddComp]);

                // push back the components multiplied by their norm
                this->fFitFuncNames.push_back("SumComp_" + addComps[iAddComp]);
                TF1 *sumComps = new TF1("SumComp_" + addComps[iAddComp], 
                        [&, this, rawComps, norms, addComps, iAddComp, onBaseline]
                        (double *x, double *pars) {
                        double sum=0.;
                        for(int iFunc=0; iFunc<onBaseline.size(); iFunc++) {
                            if(addComps[iAddComp].Contains(std::to_string(iFunc))) {
                                sum += norms[iFunc] * rawComps[iFunc]->Eval(x[0]);
                            }
                        }
                        return sum;}, this->fDrawRangeMin, this->fDrawRangeMax, 0);
                rawComps.push_back(sumComps);
                DEBUG("Number of raw components: " << rawComps.size());
                DEBUG("Eval last component: " << rawComps.back()->Eval(200));
                // determine whether, when drawing, the newly added component has to be drawn on 
                // the baseline and has to be multiplied for the global normalization constant
                bool addBaseline = true;
                bool addMultGlobNorm = true;
                for(int iFunc=0; iFunc<onBaseline.size(); iFunc++) {
                    if(addComps[iAddComp].Contains(std::to_string(iFunc))) {

                        // all the components should have the property set to true for the sum component
                        // to also have the same feature
                        if(!onBaseline[iFunc]) {
                            cout << "Not all components are to be drawn on the baseline, their sum will not be ";
                            cout << "drawn on the baseline!" << endl;
                            addBaseline = false;
                        }
                        if(!multGlobNorm[iFunc]) {
                            cout << "Not all components are to be multiplied for the global normalization constant, their sum will not be ";
                            cout << "multiplied!" << endl;
                            addMultGlobNorm = false;
                        }
                    }
                }

                if(addBaseline) {
                    onBaseline.push_back(1);
                } else {
                    onBaseline.push_back(0);
                }
                if(addMultGlobNorm) {
                    multGlobNorm.push_back(1);
                } else {
                    multGlobNorm.push_back(0);
                }
     
                norms.push_back(1);
                // kill the components to be summed by setting the normalization constants to zero 
                for(int iFunc=0; iFunc<this->fFitFuncNames.size(); iFunc++) {
                    if(addComps[iAddComp].Contains(std::to_string(iFunc))) {
                        onBaseline[iFunc] = false;
                        multNorm[iFunc] = 0.0000;
                        multGlobNorm[iFunc] = 0.0000;
                    }
                }
                DEBUG("Finished sum of components " << addComps[iAddComp]);
            }
        }

        // Define the baseline with its norm, if not indicated it is set to 1
        double baselineNorm;
        TF1 *bas = nullptr;
        DEBUG("--------------------------------");
        if(basIdx != -1) {
            int previousCompsPars = accumulate(nParsComps.begin(), std::next(nParsComps.begin(), basIdx), 0) + basIdx;
            baselineNorm = this->hParameters->GetBinContent(previousCompsPars + this->fGlobNorm + 1);
            DEBUG("Set baseline norm for function to: " << baselineNorm);
            bas = new TF1(this->fFitFuncNames[basIdx],
                [&, this, rawComps, multNorm, multGlobNorm, basIdx]
                    (double *x, double *pars) {
                       return rawComps[basIdx]->Eval(x[0]);
                    }, this->fDrawRangeMin, this->fDrawRangeMax, 0);
        } else {
            baselineNorm = 1.000000;
            bas = new TF1("fBas", "1", this->fDrawRangeMin, this->fDrawRangeMax, 0);
        }
        DEBUG("--------------------------------");

        // Define the weight of the baseline for each component, if it specified not to draw the component
        // on the baseline it will be set to zero 
        std::vector<double> onBasNorms;
        DEBUG("--------------------------------");
        for(int iFunc=0; iFunc<rawComps.size(); iFunc++) {
            if(onBaseline[iFunc]) {
                DEBUG("Set norm of the baseline for function " << iFunc << " to: " << baselineNorm);
                onBasNorms.push_back(baselineNorm);
            } else {
                DEBUG("Set norm of the baseline for function " << iFunc << " to: " << static_cast<double>(int(0)));
                onBasNorms.push_back(0.00000);
            }
        } 
        DEBUG("--------------------------------");

        // Define the global normalization constant for which every component will be multiplied, 1 if we want 
        // to draw the component as not multiplied
        std::vector<double> globNorms;
        DEBUG("--------------------------------");
        DEBUG("Number of global norms " << multGlobNorm.size());
        for(int iFunc=0; iFunc<rawComps.size(); iFunc++) {
            if(multGlobNorm[iFunc]) {
                DEBUG("Set global norm for function " << iFunc << " to: " << this->hParameters->GetBinContent(1));
                globNorms.push_back(this->hParameters->GetBinContent(1));
            } else {
                DEBUG("Set global norm for function " << iFunc << " to 1");
                globNorms.push_back(1.);
            }
        } 
        DEBUG("--------------------------------");
        DEBUG("Global norm mult or add: " << this->fMult); 
        // Define the final functions that will be drawn on the canvas
        DEBUG("Number of raw components: " << rawComps.size()); 
        for(int iRawComp=0; iRawComp<rawComps.size(); iRawComp++) {
            DEBUG("Global norm of the component: " << globNorms[iRawComp]);

            double globNorm = globNorms[iRawComp];
            double norm = norms[iRawComp];

            this->fFuncToBeDrawn.push_back(new TF1(this->fFitFuncNames[iRawComp],
                [&, this, globNorm, iRawComp, norms, funcshifts, rawComps, onBasNorms, bas] (double *x, double *pars) {
                    double term1 = onBasNorms[iRawComp] * norms[this->fBasIdx] * bas->Eval(x[0]);
                    double term2 = norms[iRawComp] * rawComps[iRawComp]->Eval(x[0]) + funcshifts[iRawComp];
                    if(iRawComp != this->fBasIdx) {
                        if(fMult) {
                            return globNorm * term1 + term2;
                        } else {
                            return globNorm * (term1 + term2);
                        }
                    } else {
                        return globNorm * term2;
                    }
                }, this->fDrawRangeMin, this->fDrawRangeMax, 0));
            DEBUG("Evaluate component: " <<  this->fFitFuncNames[iRawComp] << " value@200MeV: " << this->fFuncToBeDrawn.back()->Eval(200));
        }

        DEBUG("Raw components defined!");
    }

    /*
    */
    void Draw(std::vector<TString> legLabels, std::vector<int> colors, std::vector<double> legCoords, int linesThickness, 
              double lowRangeUser=0.0, double uppRangeUser=1.05, std::string title=";k* (MeV/c);C(k*)") {
        std::cout << endl;

        gPad->DrawFrame(fDrawRangeMin, lowRangeUser, fDrawRangeMax, uppRangeUser, title.data());

        hFitHist->GetYaxis()->SetRangeUser(lowRangeUser, uppRangeUser); 
        hFitHist->SetMarkerSize(0.1);
        hFitHist->SetMarkerStyle(24);
        hFitHist->SetMarkerColor(kBlack);
        hFitHist->SetLineColor(kBlack);
        hFitHist->SetLineWidth(3);
        hFitHist->Draw("same pe");

        std::vector<TF1 *> gaussians;
        DEBUG("Drawing " << fFuncToBeDrawn.size() << " components:");
        for(int iFuncEval=0; iFuncEval<fFuncToBeDrawn.size(); iFuncEval++) {
            this->fFuncToBeDrawn[iFuncEval]->SetNpx(1000);
            this->fFuncToBeDrawn[iFuncEval]->SetLineColor(colors[iFuncEval]);
            this->fFuncToBeDrawn[iFuncEval]->SetLineWidth(linesThickness);
            this->fFuncToBeDrawn[iFuncEval]->DrawF1(fDrawRangeMin+1,fDrawRangeMax,"same");
            DEBUG("idx: " << iFuncEval << " value: " << this->fFuncToBeDrawn[iFuncEval]->Eval(200) << " name: " << fFitFuncNames[iFuncEval]);
        }

        this->fFit->SetNpx(1000);
        this->fFit->SetLineColor(kRed);
        this->fFit->SetLineWidth(linesThickness);
        this->fFit->DrawF1(fDrawRangeMin+1,fDrawRangeMax,"same");

        // Build legend
        TLegend *legend = new TLegend(legCoords[0], legCoords[1], legCoords[2], legCoords[3]);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.045);
        legend->AddEntry(this->hFitHist, legLabels[0].Data(), "lp");
        legend->AddEntry(this->fFit, legLabels[1].Data(), "l");
        for(int iFuncEval=0; iFuncEval<fFuncToBeDrawn.size(); iFuncEval++) {
            legend->AddEntry(this->fFuncToBeDrawn[iFuncEval], legLabels[iFuncEval+2].Data(), "l");
        }
        legend->Draw("same");

        gPad->Update();
    }

    std::pair<std::vector<TF1*>, std::vector<TSpline3*>> GetFitComponents() const {
        return std::pair<std::vector<TF1*>, std::vector<TSpline3*>>({this->fFuncToBeDrawn, this->fSplines});
    };

   private:

    TH1 *hFitHist = nullptr; // Histogram to be fitted
    TH1 *hParameters = nullptr; // Histogram containing all the fit parameters
    TF1 *fFit = nullptr;
    bool fGlobNorm; 
    int fBasIdx; 
    bool fMult;         // True if the baseline is multiplicative, false if it is additive
    double fDrawRangeMin;
    double fDrawRangeMax;

    std::vector<TString> fFitFuncNames;     // Function names of fit components
    std::vector<TF1*> fFuncToBeDrawn;           // Final functions with all scalings/corrections applied. Are the ones drawn
    std::vector<TF1*> fRawFuncs;            // Fit components evaluated after the fitting
    std::vector<TSpline3*> fSplines;        // Fit components evaluated after the fitting

};

#endif  // FEMPY_DRAWFITFUNCTS_HXX_
