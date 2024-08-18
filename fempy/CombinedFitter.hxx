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
    
    // Class that defines the global chi2 for a combined fit

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
        // Set parameters of the single functions exploiting the 
        // correspondance between each fit function parameter and
        // the global fit parameters, contained in the fPars data
        // member 
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
            
            // to test the set of parameters, we evaluate the chi2 of each function taking 
            // care to maintain the correspondance between the parameters, i.e. if the 5th
            // parameter of the second fit function is the 17th parameter of the global fit
            // we'll have this information in the second vector of fPars 
            totalChi2 += (*fChi2Functions[iModel])(parModel[iModel].data());
        }
        DEBUG("Total chi2: " << totalChi2);
        return totalChi2;

	}
};

class CombinedFitter {
   public:
    CombinedFitter(std::vector<CorrelationFitter > models, std::vector<std::vector<int> > singleModelsSharedPars) {
        this->fModels = models;
        this->fSingleModelsSharedPars = singleModelsSharedPars; 
        this->fNSharedPars = singleModelsSharedPars[0].size(); 
        this->fTotalPars = 0;

        // It would be ideal not to use this data member and retrieve the fit functions
        // directly from the CorrelationFitter objects, but I was not able to do so
        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            fModelFuncts.push_back(fModels[iModel].GetFitFunction()); 
        }

        SetTotalParameters();
        SetUniqueParametersVectors();
        SetFitParameters();

    }

    TH1D *SampledHisto(TH1D *histo, double upprange, int itry) {
        int nBinsFitCF = upprange/histo->GetBinWidth(1);
        TH1D *sampledHisto = new TH1D(Form("hTry_%i", itry), "", nBinsFitCF, 
                                      histo->GetBinLowEdge(1), upprange);
        for(int iSampledBin=0; iSampledBin<sampledHisto->GetNbinsX(); iSampledBin++) {
            double binValue = histo->GetBinContent(iSampledBin+1);
            double binError = histo->GetBinError(iSampledBin+1);
            sampledHisto->SetBinContent(iSampledBin+1, gRandom->Gaus(binValue, binError));
            sampledHisto->SetBinError(iSampledBin+1, binError);
        }
        return sampledHisto;
    }

    std::vector <TH1D *> CombinedFitBootstrap(int ntries, bool difference) {
        
        DEBUG("Performing fit with non-sampled histograms!");
        ROOT::Fit::FitResult originalFitResult = CombinedFit(); 

        // Can go in DM
        DEBUG("Filling models");
        std::vector<ROOT::Math::WrappedMultiTF1> wrappedFitFuncts;
        ROOT::Fit::DataOptions opt;
        std::vector<ROOT::Fit::DataRange> fitRanges;
        std::vector<ROOT::Fit::BinData> fitBinData;

        // Set fit range
        DEBUG("Fit range");
        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            // wrappedFitFuncts.push_back(ROOT::Math::WrappedMultiTF1(*fModels[iModel].GetFitFunction(), 1));
            wrappedFitFuncts.push_back(ROOT::Math::WrappedMultiTF1(*fModelFuncts[iModel], 1));
            fitRanges.push_back(ROOT::Fit::DataRange());
            fitRanges.back().SetRange(fModels[iModel].GetLowFitRange(), fModels[iModel].GetUppFitRange());
            fitBinData.push_back(ROOT::Fit::BinData(opt, fitRanges.back()));
        }

        // Initialize histogram for saving fit parameters
        std::vector <TH1D *> hBTHistos;
        for(int iPar=0; iPar<this->fInitPars.size(); iPar++) {
            double histoMean = originalFitResult.Parameter(iPar);
            if(originalFitResult.IsParameterFixed(iPar)) {
                std::string parName = originalFitResult.ParName(iPar) + "_fixed";
                double histoBound = 0.1;
                hBTHistos.push_back(new TH1D(parName.c_str(), parName.c_str(), 
                                          1000, histoMean - 5*histoBound, histoMean + 5*histoBound));
            } else {
                std::string parName = originalFitResult.ParName(iPar);
                double histoBound = originalFitResult.ParError(iPar);
                hBTHistos.push_back(new TH1D(parName.c_str(), parName.c_str(), 
                                          1000, histoMean - 5*histoBound, histoMean + 5*histoBound));
            }
        }

        if(difference) {
            for(int iModel=0; iModel<this->fModels.size(); iModel++) {
                double binWidth = this->fModels[iModel].GetFitHisto()->GetBinWidth(1);
                int nBinsFitCF = fModels[iModel].GetUppFitRange()/binWidth;
                for(int iBin=0; iBin<nBinsFitCF; iBin++) {
                    int binCenter = static_cast<int>((iBin * binWidth) + (binWidth / 2));
                    double fitDifference = this->fModels[iModel].GetFitHisto()->GetBinContent(iBin+1) - 
                                           this->fModels[iModel].GetFitFunction()->Eval(this->fModels[iModel].GetFitHisto()->GetBinCenter(iBin+1));
                    hBTHistos.push_back(new TH1D(Form("Model%i_Subtraction_%iMeV", iModel, binCenter), 
                                                 Form("Model%i_Subtraction_%iMeV", iModel, binCenter), 
                                                 10000, fitDifference - 5*fitDifference, fitDifference + 5*fitDifference));
                }
            }
        }

        // hBTHistos.push_back(new TH1D("Chi2/DOF", "Chi2/DOF", 1000, 0, 100));

        for(int iTry=0; iTry<ntries; iTry++) {
            std::vector<TH1D *> sampledHistos;
            std::vector<ROOT::Fit::Chi2Function> chi2Functions;
            GlobalChi2 global_chi2(fSingleModelsUniquePars, fSingleModelsSharedPars);
            int nPar = 0;
            for(int iModel=0; iModel<this->fModels.size(); iModel++) {
        
                sampledHistos.push_back(SampledHisto(static_cast<TH1D *>(fModels[iModel].GetFitHisto()), 
                                                     fModels[iModel].GetUppFitRange(), iTry));
                DEBUG("GetNbins histo: " << sampledHistos.back()->GetNbinsX());
                ROOT::Fit::FillData(fitBinData[iModel], sampledHistos.back());
                nPar += fitBinData[iModel].Size();

                DEBUG("Global chi2");
                chi2Functions.push_back(ROOT::Fit::Chi2Function(fitBinData[iModel], wrappedFitFuncts[iModel]));

            }

            for(int iModel=0; iModel<this->fModels.size(); iModel++) {
                global_chi2.AddChi2Function(chi2Functions[iModel]);
            }
                
            DEBUG("Fitter");
            ROOT::Fit::Fitter fitter;
            SetupGlobalFitter(&fitter);

            DEBUG("Setup fit");
            DEBUG("Fit");
            DEBUG("nPar: " << nPar);
            DEBUG("Total pars: " << this->fTotalPars);
            fitter.FitFCN(this->fTotalPars, global_chi2, nullptr, nPar, 1);
            ROOT::Fit::FitResult result = fitter.Result();

            for(int iPar=0; iPar<this->fInitPars.size(); iPar++) {
                hBTHistos[iPar]->Fill(result.Parameter(iPar));
            }

            // separate fit results
            for(int iModel=0; iModel<this->fModels.size(); iModel++) {
                fModelFuncts[iModel]->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
                // fModels[iModel].GetFitFunction()->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
            }

            if(difference) {
                for(int iModel=0; iModel<this->fModels.size(); iModel++) {
                    double binWidth = this->fModels[iModel].GetFitHisto()->GetBinWidth(1);
                    int nBinsFitCF = fModels[iModel].GetUppFitRange()/binWidth;
                    for(int iBin=0; iBin<nBinsFitCF; iBin++) {
                        hBTHistos[this->fInitPars.size() + iBin + iModel*nBinsFitCF]->Fill(
                            this->fModels[iModel].GetFitHisto()->GetBinContent(iBin+1) - 
                            fModelFuncts[iModel]->Eval(this->fModels[iModel].GetFitHisto()->GetBinCenter(iBin+1)));
                    }
                }
            }
        
            sampledHistos.clear();
            chi2Functions.clear();
        }

        if(difference) {
            for(int iModel=0; iModel<this->fModels.size(); iModel++) {
                double binWidth = this->fModels[iModel].GetFitHisto()->GetBinWidth(1);
                int nBinsFitCF = fModels[iModel].GetUppFitRange()/binWidth;

                // Restore original fit results to evaluate the subtraction
                fModelFuncts[iModel]->SetFitResult(originalFitResult, fSingleModelsUniquePars[iModel].data());

                TH1D *hSubtraction = new TH1D(Form("hDifference_Model%i", iModel), Form("hDifference_Model%i", iModel), nBinsFitCF, 0, this->fModels[iModel].GetUppFitRange());
                TH1D *hYields = new TH1D(Form("hYields_Model%i_stat", iModel), Form("hYields_Model%i_stat", iModel), nBinsFitCF, 0, this->fModels[iModel].GetUppFitRange());
                TH1D *hMeans = new TH1D(Form("hMeans_Model%i_stat", iModel), Form("hMeans_Model%i_stat", iModel), nBinsFitCF, 0, this->fModels[iModel].GetUppFitRange());
                TH1D *hStdDevs = new TH1D(Form("hStdDevs_Model%i_stat", iModel), Form("hStdDevs_Model%i_stat", iModel), nBinsFitCF, 0, this->fModels[iModel].GetUppFitRange());
                TH1D *hRelUnc = new TH1D(Form("hRelUnc_Model%i_stat", iModel), Form("hRelUnc_Model%i_stat", iModel), nBinsFitCF, 0, this->fModels[iModel].GetUppFitRange());
                
                for(int iBin=0; iBin<nBinsFitCF; iBin++) {
                    int iHistoIdx = this->fTotalPars + iModel * nBinsFitCF + iBin + 1;
                    TF1 *gaus = new TF1("gaus", "gaus", hBTHistos[iHistoIdx]->GetBinLowEdge(1),
                                        hBTHistos[iHistoIdx]->GetBinLowEdge(hBTHistos[iHistoIdx]->GetNbinsX()) + hBTHistos[iHistoIdx]->GetBinWidth(1));
                    gaus->SetParameter(1, hBTHistos[iHistoIdx]->GetBinCenter(hBTHistos[iHistoIdx]->GetMaximumBin()) );
                    gaus->SetParameter(2, hBTHistos[iHistoIdx]->GetRMS() );
                    hBTHistos[iBin]->Fit(gaus, "SMRL+", "");
                    hSubtraction->SetBinContent(iBin+1,  this->fModels[iModel].GetFitHisto()->GetBinContent(iBin+1) - 
                                                         fModelFuncts[iModel]->Eval(this->fModels[iModel].GetFitHisto()->GetBinCenter(iBin+1)));
                    hYields->SetBinContent(iBin+1, gaus->GetParameter(0));
                    hMeans->SetBinContent(iBin+1, gaus->GetParameter(1));
                    hStdDevs->SetBinContent(iBin+1, gaus->GetParameter(2));
                    hRelUnc->SetBinContent(iBin+1, hStdDevs->GetBinContent(iBin+1) / hSubtraction->GetBinContent(iBin+1));
                    hSubtraction->SetBinError(iBin+1, hStdDevs->GetBinContent(iBin+1) );
                }
                hBTHistos.push_back(hYields);
                hBTHistos.push_back(hMeans);
                hBTHistos.push_back(hStdDevs);
                hBTHistos.push_back(hRelUnc);
                hBTHistos.push_back(hSubtraction);
            }
        }

        return hBTHistos;

    }

    ROOT::Fit::FitResult CombinedFit() {

        // Can go in DM
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

        DEBUG("Global chi2");
        std::vector<ROOT::Fit::Chi2Function> chi2Functions;
        for(int iFit=0; iFit<this->fModels.size(); iFit++) {
            chi2Functions.push_back(ROOT::Fit::Chi2Function(fitBinData[iFit], wrappedFitFuncts[iFit]));
        }
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

        // separate fit results
        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            fModelFuncts[iModel]->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
            // fModels[iModel].GetFitFunction()->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
        }

        return result;
        // // TODO: define Chi2
    }

    std::vector<TH1D *> CombinedFitDifference() {

        // Can go in DM
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

        DEBUG("Global chi2");
        std::vector<ROOT::Fit::Chi2Function> chi2Functions;
        for(int iFit=0; iFit<this->fModels.size(); iFit++) {
            chi2Functions.push_back(ROOT::Fit::Chi2Function(fitBinData[iFit], wrappedFitFuncts[iFit]));
        }
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

        // separate fit results
        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            fModelFuncts[iModel]->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
            // fModels[iModel].GetFitFunction()->SetFitResult(result, fSingleModelsUniquePars[iModel].data());
        }

        std::vector<TH1D *> hDifferences;
        for(int iModel=0; iModel<this->fModels.size(); iModel++) {
            double binWidth = this->fModels[iModel].GetFitHisto()->GetBinWidth(1);
            int nBinsFitCF = fModels[iModel].GetUppFitRange()/binWidth;
            TH1D *hSubtraction = new TH1D(Form("hSubtraction_Model%i", iModel), 
                                          Form("hSubtraction_Model%i", iModel), nBinsFitCF, 0, 
                                          this->fModels[iModel].GetUppFitRange());
            for(int iBin=0; iBin<nBinsFitCF; iBin++) {
                hSubtraction->SetBinContent(iBin+1,  this->fModels[iModel].GetFitHisto()->GetBinContent(iBin+1) - 
                                                     fModelFuncts[iModel]->Eval(this->fModels[iModel].GetFitHisto()->GetBinCenter(iBin+1)));
            }
            hDifferences.push_back(hSubtraction);
        }

        return hDifferences;
        // // TODO: define Chi2
    }

    void SetupGlobalFitter(ROOT::Fit::Fitter *fitter) { 
        DEBUG("Number of parameters to be initialized: " << this->fInitPars.size());
        std::vector<double> dummyInit(this->fTotalPars, 0.0);

        // trivial initialization with dummyInit, which is a vector of zeros, 
        // then setup each parameter specifically according to fInitPar
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

    int GetModelUniqueParams(int imodel) {
        return this->fSingleModelsUniquePars[imodel].size() - this->fNSharedPars;
    }

    int GetSharedParams() {
        return this->fNSharedPars;
    }

    int GetTotalParams() {
        return this->fTotalPars;
    }

   private:

    void SetTotalParameters() {

        // For each model we sum the number of non-shared parameters, 
        // then at the end we sum once the number of the shared ones
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

        // To implement a combined fit of two functions, we need a vector for each single fit 
        // which stores the indices of the parameters in the global combined fit associated 
        // to the fit function parameters

        // i.e. if we execute a combined fit of two pol3 sharing the linear coefficient, the 
        // total number of parameters will be 7 = 3(unique)*2 + 1(shared). The global fit 
        // parameters will be 7, thus it is needed to specify which of the 7 parameters corre-
        // sponds to the p0, p1, p2, p3 of the first and second pol

        // In this class, the vector storing the global fit parameters is arranged in the 
        // following way: fInitPars: {(unique pars 1st comp), ... ,(unique pars nth comp), (shared pars)}  

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

        // The shared pars are put after the unique pars of 
        // all models which are considered for the combined fit
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
    std::vector<std::vector<int> > fSingleModelsSharedPars;     // vector of vectors where the latter contain
                                                                // the indices of the parameters that are shared
                                                                // for each specific model 
    std::vector<std::vector<int> > fSingleModelsUniquePars;     // vector of vectors where the latter contain
                                                                // the indices of the model parameter in the 
                                                                // global list of parameter for the combined fit
    std::vector<std::tuple<std::string, double, double, double> > fInitPars;    // initialization of the combined
                                                                                // fit parameters, {"par_name", init,
                                                                                // low_lim, upp_lim}
};

#endif  // FEMPY_COMBINEDFITTER_HXX_
