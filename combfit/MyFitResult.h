#ifndef COMBFIT_MYFITRESULT_H_
#define COMBFIT_MYFITRESULT_H_

#include <Fit/BinData.h>
#include <Fit/Chi2FCN.h>
#include <Fit/FitResult.h>
#include <Fit/Fitter.h>
#include <HFitInterface.h>
#include <Math/WrappedMultiTF1.h>

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <utility>


class MyFitResult : public ROOT::Fit::FitResult {
 public:
    MyFitResult(const ROOT::Fit::FitResult &fitresult, TF1 *fit, int *idx = NULL) : ROOT::Fit::FitResult(fitresult) {
        std::shared_ptr<IModelFunction> fitFCN(new ROOT::Math::WrappedMultiTF1(*fit, 1));
        SetFitFunction(fitFCN);
        UseParameters(idx);
    }
    void SetFitFunction(std::shared_ptr<IModelFunction> func) { fFitFunc = func; }
    void UseParameters(int *idx) {
        if (idx == NULL) return;

        int npars = fFitFunc->NPar();

        // replace the member vectors and maps with new ones that contain only the given parameters
        // parameter indices as indices of the overall Fitter!
        std::vector<double> params;
        std::vector<double> errors;
        std::vector<string> par_names;
        std::vector<pair<double, double>> param_bounds;
        std::vector<double> cov_matrix;
        std::vector<double> global_cc;
        std::map<unsigned int, bool> fixed_params;
        std::map<unsigned int, unsigned int> bound_params;
        std::map<unsigned int, pair<double, double>> minos_errors;

        cov_matrix.reserve(npars * (npars + 1) / 2);
        for (int jpar = 0; jpar < npars; ++jpar) {
            int ipar = idx[jpar];

            params.push_back(fParams[ipar]);
            errors.push_back(fErrors[ipar]);
            par_names.push_back(fParNames[ipar]);
            global_cc.push_back(fGlobalCC[ipar]);
            if (fFixedParams.find(ipar) != fFixedParams.end())
                fixed_params.insert(pair<unsigned int, bool>(jpar, fFixedParams[ipar]));
            if (fBoundParams.find(ipar) != fBoundParams.end()) {
                bound_params.insert(pair<unsigned int, unsigned int>(jpar, param_bounds.size()));
                int kpar = fBoundParams[ipar];
                param_bounds.push_back(fParamBounds[kpar]);
            }
            if (fMinosErrors.find(ipar) != fMinosErrors.end())
                minos_errors.insert(pair<unsigned int, pair<double, double>>(jpar, fMinosErrors[ipar]));
        }

        for (int ipar = 0; ipar < npars; ++ipar) {
            int kpar = idx[ipar];
            for (int jpar = 0; jpar <= ipar; ++jpar) {
                int lpar = idx[jpar];
                cov_matrix.push_back(CovMatrix(kpar, lpar));
            }
        }

        fParams = params;
        fErrors = errors;
        fParNames = par_names;
        fParamBounds = param_bounds;
        fCovMatrix = cov_matrix;
        fGlobalCC = global_cc;
        fFixedParams = fixed_params;
        fBoundParams = bound_params;
        fMinosErrors = minos_errors;
    }
};

#endif  // COMBFIT_MYFITRESULT_H_
