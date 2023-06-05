'''
This macro only works for ROOT>=6.26. Kept for reference/consultation

'''

from cgi import test
from ROOT import Math, TF1, TF2, TRandom2, TCanvas, addressof


fExp = TF1("fExp", "(x-0.5)**2", 0, 1)
fTheo = TF1("fTheo", "[0]+[1]*x", 0, 1)


def customChi2(fitpar):
    p0 = fitpar[0]
    p1 = fitpar[1]

    fTheo.SetParameter(0, p0)
    fTheo.SetParameter(1, p1)
    # print(<<"p0: " << p0 <<"p1: " << p1<< "    int theo: "<<fExp.Integral(0, 1)<< "    int exp: "<<fTheo.Integral(0, 1))
    fDelta = TF1("f", "(fExp-fTheo)**2", 0, 1)
    # double diff = (fExp.Integral(0, 1) - fTheo.Integral(0, 1))
    return fDelta.Integral(0, 1)


def customChi2Formula(xx, par):
    # p0 = xx[0]
    # p1 = xx[1]

    # fTheo.SetParameter(0, p0)
    # fTheo.SetParameter(1, p1)
    # # print(<<"p0: " << p0 <<"p1: " << p1<< "    int theo: "<<fExp.Integral(0, 1)<< "    int exp: "<<fTheo.Integral(0, 1))
    # double chi2 = (fExp.Integral(0, 1) - fTheo.Integral(0, 1))*(fExp.Integral(0, 1) - fTheo.Integral(0, 1))
    return customChi2(xx)


def testCustomChi2Minimization(minName="Minuit2",algoName="", randomSeed=42):
    minimum = Math.Factory.CreateMinimizer("Minuit", algoName)
    # set tolerance , etc...
    minimum.SetMaxFunctionCalls(1000000)  # for Minuit/Minuit2
    minimum.SetMaxIterations(10000)      # for GSL
    minimum.SetTolerance(0.001)
    minimum.SetPrintLevel(2)
    # create function wrapper for minimizer
    # a IMultiGenFunction type
    # Math.Functor f(&RosenBrock, 2)
    # Math.Functor f(&customChi2, 1)

    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    f = Math.Functor(customChi2, 2)
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    ### !!! DOESNT WORK WITH ROOT <= 6.24
    
    # f = Math.Functor(addressof(customChi2), 2)
    # double step[1] = {0.01}
    step = [0.01, 0.01]
    # starting point
    # double variable[1] = {0.5}
    variable = [0.5, 0.5]

    minimum.SetFunction(f)
    # Set the free variables to be minimized !
    minimum.SetVariable(0, "p0", variable[0], step[0])
    minimum.SetVariable(1, "p1", variable[1], step[1])
    # do the minimization
    minimum.Minimize()
    xs = minimum.X()

    c = TCanvas("c", "c", 600, 600)
    fExp.Draw("")
    fTheo.Draw("same")

    cc = TCanvas("c2", "c", 600, 600)
    fChi2 = TF2("fChi2", customChi2Formula, -2, 2, -2, 2)
    fChi2.Draw("colz")

    return 0

if __name__ == '__main__':
    testCustomChi2Minimization()