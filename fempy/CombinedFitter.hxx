#ifndef FEMPY_COMBINEDFITTER_HXX_
#define FEMPY_COMBINEDFITTER_HXX_

#include <map>
#include <string>
#include <tuple>
#include <stdexcept>
#include <algorithm> 

#include "TF1.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TFitResultPtr.h"
// #include "FitFunctions.cxx"
#include "CorrelationFitter.hxx"

#include "Math/WrappedMultiTF1.h"
#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "Fit/Fitter.h"
#include "HFitInterface.h"

#if LOG_LEVEL_COMBFIT
#define DEBUG(msg) std::cout << msg << std::endl
#else
#define DEBUG(msg)
#endif

class GlobalChi2 {
    public:
    std::vector<ROOT::Math::IMultiGenFunction *> fChi2Functions; 
    std::vector<std::vector<int> > fModelsSharedParsIdxs;
    std::vector<std::vector<int>> fPars;

	GlobalChi2(std::vector<std::vector<int>> pars, std::vector<std::vector<int>> sharedidxs) {
	    fPars = pars;
        fModelsSharedParsIdxs = sharedidxs;
        DEBUG("fPars size: " << fPars.size());
    }

    void AddChi2Function(ROOT::Math::IMultiGenFunction &func) {
        fChi2Functions.push_back(&func);
    }

	double operator() (const double *par) const {
        // Set parameters of the single functions, if the number
        // of parameters is different in the two cases two loops
        // are required
        double totalChi2 = 0.;
        std::vector<double> parModel[fPars.size()];
        for(int iModel=0; iModel<fPars.size(); iModel++) {
            DEBUG("");
            DEBUG("MODEL: " << iModel);
            DEBUG("Number of parameters: " << fPars[iModel].size());
            for (int i=0; i<fPars[iModel].size(); i++) {
                DEBUG("iPar" << this->fPars[iModel][i] << ": " << par[this->fPars[iModel][i]]);
	    	    parModel[iModel].push_back(par[this->fPars[iModel][i]]);
	        }
            DEBUG("");
            totalChi2 += (*fChi2Functions[iModel])(parModel[iModel].data());
        }
        DEBUG("Total chi2: " << totalChi2);
        return totalChi2;

	}
};

class CombinedFitter {
   public:
    CombinedFitter(std::vector<CorrelationFitter > models, std::vector<std::vector<int> > singleModelsSharedPars) {
        DEBUG("CIAO COMBINED");
        this->fModels = models;
        this->fSingleModelsSharedPars = singleModelsSharedPars; 
        this->fNSharedPars = singleModelsSharedPars[0].size(); 
        this->fTotalPars = 0;

        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            fModelFuncts.push_back(fModels[iModel].GetFitFunction()); 
        }

        SetTotalParameters();
        SetUniqueParametersVectors();
        SetFitParameters();

    }

    void CombinedFit() {
        DEBUG("Filling models");
        std::vector<ROOT::Math::WrappedMultiTF1> wrappedFitFuncts;
        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            wrappedFitFuncts.push_back(ROOT::Math::WrappedMultiTF1(*fModelFuncts[iModel], 1));
            // wrappedFitFuncts.push_back(ROOT::Math::WrappedMultiTF1(*fModels[iModel].GetFitFunction(), 1));
        }

        // Set fit range
        DEBUG("Fit range");
        ROOT::Fit::DataOptions opt;
        std::vector<ROOT::Fit::DataRange> fitRanges;
        std::vector<ROOT::Fit::BinData> fitBinData;
        for(int iFitRange=0; iFitRange<this->fModels.size(); iFitRange++) {
            fitRanges.push_back(ROOT::Fit::DataRange());
            fitRanges.back().SetRange(fModels[iFitRange].GetLowFitRange(), fModels[iFitRange].GetUppFitRange());
            fitBinData.push_back(ROOT::Fit::BinData(opt, fitRanges.back()));
            DEBUG("GetNbins histo: " << fModels[iFitRange].GetFitHisto()->GetNbinsX());
            ROOT::Fit::FillData(fitBinData.back(), fModels[iFitRange].GetFitHisto());
        }

        DEBUG("Chi squared");
        std::vector<ROOT::Fit::Chi2Function> chi2Functions;
        for(int iFit=0; iFit<this->fModels.size(); iFit++) {
            chi2Functions.push_back(ROOT::Fit::Chi2Function(fitBinData[iFit], wrappedFitFuncts[iFit]));
        }


        DEBUG("Global chi2");
        GlobalChi2 global_chi2(fSingleModelsUniquePars, fSingleModelsSharedPars);
        for(int iChi2Fcn=0; iChi2Fcn<chi2Functions.size(); iChi2Fcn++) {
            global_chi2.AddChi2Function(chi2Functions[iChi2Fcn]);
        }

        DEBUG("Fitter");
        ROOT::Fit::Fitter fitter;
        SetupGlobalFitter(&fitter);
        
        int nPar = 0;
        for(int iModel=0; iModel<fModels.size(); iModel++) {
            nPar += fitBinData[iModel].Size();
        }

        DEBUG("Setup fit");
        DEBUG("Fit");
        DEBUG("nPar: " << nPar);
        DEBUG("Total pars: " << this->fTotalPars);
        fitter.FitFCN(this->fTotalPars, global_chi2, nullptr, nPar, 1);
        ROOT::Fit::FitResult result = fitter.Result();
        DEBUG("Combined fit Chi2 = " << result.Chi2());
        // separate fit results
        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            fModelFuncts[iModel]->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
            // fModels[iModel].GetFitFunction()->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
        }
    }

    void SetupGlobalFitter(ROOT::Fit::Fitter *fitter) { 
        DEBUG("Number of parameters to be initialized: " << this->fInitPars.size());
        std::vector<double> dummyInit(this->fTotalPars, 0.0);
        fitter->Config().SetParamsSettings(this->fTotalPars, dummyInit.data());  // number of total parameters, list of init values
        for(int iInitPar=0; iInitPar<this->fInitPars.size(); iInitPar++) {
            DEBUG("InitPar" << iInitPar << ": [" << std::get<0>(fInitPars[iInitPar]) << ", " << std::get<1>(fInitPars[iInitPar]) << ", " <<
                  std::get<2>(fInitPars[iInitPar]) << ", " << std::get<3>(fInitPars[iInitPar]) << "]");
            fitter->Config().ParSettings(iInitPar).SetName(std::get<0>(fInitPars[iInitPar]));
            fitter->Config().ParSettings(iInitPar).SetValue(std::get<1>(fInitPars[iInitPar]));
            if(std::get<2>(fInitPars[iInitPar]) >= std::get<3>(fInitPars[iInitPar])) {
                fitter->Config().ParSettings(iInitPar).Fix();
            } else {
                fitter->Config().ParSettings(iInitPar).SetValue(std::get<1>(fInitPars[iInitPar]));
                fitter->Config().ParSettings(iInitPar).SetLimits(std::get<2>(fInitPars[iInitPar]), std::get<3>(fInitPars[iInitPar]));
            }
        }

        fitter->Config().MinimizerOptions().SetPrintLevel(1);
        fitter->Config().SetMinimizer("Minuit2", "Migrad");
    
    }

   private:

    void SetTotalParameters() {
        for(int iModel=0; iModel<fSingleModelsSharedPars.size(); iModel++) {
            this->fModels[iModel].BuildFitFunction();
            DEBUG("Parameters of the model: " << this->fModels[iModel].GetFitFunction()->GetNpar());
            DEBUG("Shared parameters: " << fSingleModelsSharedPars[iModel].size());
            this->fTotalPars += this->fModels[iModel].GetFitFunction()->GetNpar() - 
                                fSingleModelsSharedPars[iModel].size();
        }
        this->fTotalPars += fNSharedPars;
        DEBUG("Total parameters: " << this->fTotalPars);
    }

    void SetUniqueParametersVectors() {
        int previousModelsUniquePars = 0;
        DEBUG("Previous models unique parameters: " << previousModelsUniquePars);
        for(int iModel=0; iModel<fSingleModelsSharedPars.size(); iModel++) {
            DEBUG("Model " << iModel);
            std::vector<int> modelUniquePars;
            int nModelUniquePars = this->fModels[iModel].GetFitFunction()->GetNpar() - 
                                   fSingleModelsSharedPars[iModel].size(); 
            int iSharedPar=0;
            for(int iPar=0; iPar<this->fModels[iModel].GetFitFunction()->GetNpar(); iPar++) {
                DEBUG("iPar " << iPar);
                if(iPar == this->fSingleModelsSharedPars[iModel][iSharedPar]-1) {
                    DEBUG("Shared, push back n par " << this->fTotalPars + iSharedPar - this->fNSharedPars);
                    modelUniquePars.push_back(this->fTotalPars + iSharedPar - this->fNSharedPars);
                    iSharedPar++;
                } else {
                    DEBUG("Adding unique parameter: " << iPar + previousModelsUniquePars - iSharedPar << " to model " << iModel);
                    DEBUG("Unique, push back n par " << iPar + previousModelsUniquePars - iSharedPar);
                    modelUniquePars.push_back(iPar + previousModelsUniquePars - iSharedPar);
                }
            }
            fSingleModelsUniquePars.push_back(modelUniquePars);
            previousModelsUniquePars += nModelUniquePars;
            DEBUG("Previous models unique parameters: " << previousModelsUniquePars);
            DEBUG("");
        }

        for(int iModel=0; iModel<fSingleModelsSharedPars.size(); iModel++) {
            DEBUG("Model: " << iModel);
            DEBUG("[ ");
            for(int iIdx=0; iIdx<fSingleModelsUniquePars[iModel].size(); iIdx++) {
                DEBUG(fSingleModelsUniquePars[iModel][iIdx] << ", ");
            }
            DEBUG("]");
            DEBUG("");
        }
    }

    void SetFitParameters() {

        std::tuple<std::string, double, double, double> initSharedPars[fNSharedPars];
        double lowerLim, upperLim;
        for(int iModel=0; iModel<fModels.size(); iModel++) { 
            DEBUG("iModel: " << iModel); 
            int iSharedPar = 0;
            for(int iPar=0; iPar<this->fModels[iModel].GetFitFunction()->GetNpar(); iPar++) {
                this->fModels[iModel].GetFitFunction()->GetParLimits(iPar, lowerLim, upperLim);
                std::tuple<std::string, double, double, double> initPar = 
                    {this->fModels[iModel].GetFitFunction()->GetParName(iPar),
                     this->fModels[iModel].GetFitFunction()->GetParameter(iPar),
                     lowerLim, upperLim};    
                DEBUG("iPar" << iPar << ": "; 
                      cout << this->fModels[iModel].GetFitFunction()->GetParName(iPar) << ", "; 
                      cout << this->fModels[iModel].GetFitFunction()->GetParameter(iPar) << ", "
                           << lowerLim << ", " << upperLim); 
                if(std::find(fSingleModelsSharedPars[iModel].begin(), fSingleModelsSharedPars[iModel].end(), iPar+1)
                   != fSingleModelsSharedPars[iModel].end()) {
                    DEBUG("iPar" << iPar << " is shared!");
                    initSharedPars[iSharedPar] = initPar;
                    iSharedPar++;
                } else {
                    fInitPars.push_back(initPar);
                }
            }
        }
        for(int iSharedPar=0; iSharedPar<fNSharedPars; iSharedPar++) {
            fInitPars.push_back(initSharedPars[iSharedPar]);
        }
    }

    void Debug() {
    }

    std::vector<CorrelationFitter > fModels;        // vector storing the models used for the fitting
    int fTotalPars;                                 // sum of unique pars of single fits and shared ones 
    int fNSharedPars;                               // number of shared parameters 
    std::vector<TF1 *> fModelFuncts;
    std::vector<std::vector<int> > fSingleModelsSharedPars;
    std::vector<std::vector<int> > fSingleModelsUniquePars;
    std::vector<std::tuple<std::string, double, double, double> > fInitPars;
};

#endif  // FEMPY_COMBINEDFITTER_HXX_
