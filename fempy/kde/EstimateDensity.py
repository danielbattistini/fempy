'''
Script to create a tree dataset from binned dataset
'''
import numpy as np
import sys
import itertools

from utils.io import GetObjectFromFile
from utils.handle import GetFemtoDreamPairId
from utils.correlation_handler import CorrelationHandler

from ROOT import TFile, TTree, gRandom, TH1D, TKDE, RDataFrame, TF1, SetOwnership, gROOT, TCanvas, nullptr, EnableImplicitMT, TList


def ComputeMultWeights(hhSE, hhME):
    weights = []

    nMultBinsSE = hhSE.GetNbinsY()
    nMultBinsME = hhME.GetNbinsY()

    if nMultBinsSE != nMultBinsME:
        print('not same mult bins. Exit!')
        sys.exit()

    nMultBins = nMultBinsSE

    nSEEntries = hhSE.GetEntries()
    nMEEntries = hhME.GetEntries()
    for iMultBin in range(nMultBins):
        hSEmult = hhSE.ProjectionX(f'{hhSE.GetName()}_projy_{iMultBin}', iMultBin, iMultBin)
        hMEmult = hhME.ProjectionX(f'{hhME.GetName()}_projy_{iMultBin}', iMultBin, iMultBin)

        nSEProjEntries = hSEmult.GetEntries()
        nMEProjEntries = hMEmult.GetEntries()

        if nMEProjEntries == 0:
            weights.append(0)
        else:
            # relAbSE = nSEProjEntries / nSEEntries
            # relAbME = nMEProjEntries / nMEEntries
            # weights.append(relAbSE/relAbME) # gets up to 50

            weights.append(nSEProjEntries / nMEProjEntries)  # the GF way

    return weights


def EstimateDensity(
        kernelSetup,
        inFileName='/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults_merged.root',
        oFileName='/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/densities/Densities_test2.root'):
    ''' Estimate the PDF of the SE/ME '''

    '''KDE options:
    Binning: kUnbinned, kRelaxedBinning, kForcedBinning
    Iteration: kAdaptive, kFixed
    KernelType: kGaussian, kEpanechnikov, kBiweight, kCosineArch, kUserDefined, kTotalKernels
    Mirror: kNoMirror, kMirrorLeft, kMirrorRight, kMirrorBoth, kMirrorAsymLeft, kMirrorRightAsymLeft, kMirrorAsymRight, kMirrorLeftAsymRight, kMirrorAsymBoth
    '''

    EnableImplicitMT()

    kStarMax = 3000
    pairIds = ['sc']
    eventIds = ['SE', 'ME']
    regionIds = ['sgn', 'sbl', 'sbr']
    fracToKeep = 1
    kdeSetupStr = f'KernelType:{kernelSetup["kernel"]};' \
                  f'Iteration:{kernelSetup["iter"]};' \
                  f'Mirror:{kernelSetup["mirror"]};' \
                  f'Binning:{kernelSetup["binning"]}'
    oFileName = oFileName[:-5] + f'_kern-{kernelSetup["kernel"]}' \
                                 f'_iter-{kernelSetup["iter"]}' \
                                 f'_mirror-{kernelSetup["mirror"]}' \
                                 f'_binning-{kernelSetup["binning"]}.root'
    if fracToKeep < 0.99:
        oFileName = oFileName[:-5] + f'_red-{fracToKeep}.root'

    kdeEstimates = {}

    inFile = TFile(inFileName)
    oFile = TFile(oFileName, 'recreate')

    for (pId, rId) in itertools.product(pairIds, regionIds):
        dId = f'SE_{pId}_{rId}'

        def cfRelUncUpperFcn(x, par):
            upper = kdeEstimates[f'{dId}_upper'].Eval(x[0])
            centr = kdeEstimates[f'{dId}_centr'].Eval(x[0])
            return (upper - centr) / centr

        def cfRelUncLowerFcn(x, par):
            centr = kdeEstimates[f'{dId}_centr'].Eval(x[0])
            lower = kdeEstimates[f'{dId}_lower'].Eval(x[0])
            return (centr - lower) / centr

        print(dId)

        # same-event
        rdfSE = RDataFrame(f'treeDpi_{dId}', inFileName)
        if fracToKeep < 0.99:
            rdfSE = rdfSE.Filter(f'gRandom->Rndm() < {fracToKeep}')  # todo: remove
        kStarSE = rdfSE.AsNumpy(['kStar'])['kStar']
        nPairsSE = len(kStarSE)

        kdeSE = TKDE(nPairsSE, kStarSE, 0.0, kStarMax, kdeSetupStr)

        kdeEstimates[f'{dId}_centr'] = kdeSE.GetFunction(1000, 0, kStarMax).Clone()
        kdeEstimates[f'{dId}_centr'].SetTitle(" ;k* (MeV/#it{c});Density")
        kdeEstimates[f'{dId}_centr'].SetName(f'f{dId}_centr')
        kdeEstimates[f'{dId}_centr'].Write()

        kdeEstimates[f'{dId}_upper'] = kdeSE.GetUpperFunction(0.683, 1000, 0, kStarMax)
        kdeEstimates[f'{dId}_upper'].SetTitle(" ;k* (MeV/#it{c});Density")
        kdeEstimates[f'{dId}_upper'].SetName(f'f{dId}_upper')
        kdeEstimates[f'{dId}_upper'].Write()

        kdeEstimates[f'{dId}_lower'] = kdeSE.GetLowerFunction(0.683, 1000, 0, kStarMax)
        kdeEstimates[f'{dId}_lower'].SetTitle(" ;k* (MeV/#it{c});Density")
        kdeEstimates[f'{dId}_lower'].SetName(f'f{dId}_lower')
        kdeEstimates[f'{dId}_lower'].Write()

        kdeEstimates[f'{dId}_relunc_upper'] = TF1(f'f{dId}_relunc_upper', cfRelUncUpperFcn, 0, kStarMax, 0)
        kdeEstimates[f'{dId}_relunc_upper'].SetTitle(" ;k* (MeV/#it{c});Relative error")
        kdeEstimates[f'{dId}_relunc_upper'].Write()

        kdeEstimates[f'{dId}_relunc_lower'] = TF1(f'f{dId}_relunc_lower', cfRelUncLowerFcn, 0, kStarMax, 0)
        kdeEstimates[f'{dId}_relunc_lower'].SetTitle(" ;k* (MeV/#it{c});Relative error")
        kdeEstimates[f'{dId}_relunc_lower'].Write()
        # auto cfRelUncLowerFcn = [&](double *x, double *) {
        #     double lower = mapDensities[TString(Form("SE_%s_lower", baseName.Data()))]->EvalPar(x, nullptr);
        #     double centr = mapDensities[TString(Form("SE_%s_centr", baseName.Data()))]->EvalPar(x, nullptr);
        #     return (centr - lower) / centr;
        # };
        # auto fff = new TF1(Form("%s_relunc_upper", baseName.Data()), cfRelUncUpperFcn, kStarMin, kStarMax, 0)
        # # auto fffL = new TF1(Form("%s_relunc_lower", baseName.Data()), cfRelUncLowerFcn, kStarMin, kStarMax, 0);
        # // // fff -> Draw()
        # # // mapCorrelationsUnc[mapName]->SetTitle(";#it{k}* (MeV/#it{c});#it{C}");
        # fff -> SetNpx(1000)
        # fff -> Write()
        # // mapCorrelationsUnc[mapName]->SetTitle(";#it{k}* (MeV/#it{c});#it{C}");
        # fffL->SetNpx(1000);
        # fffL->Write();

        ###################################
        ###################################
        ###################################
        ###################################

        dId = f'ME_{pId}_{rId}'

        print(dId)

        # mixed-event

        rdfME = RDataFrame(f'treeDpi_{dId}', inFileName)
        if fracToKeep < 0.99:
            rdfME = rdfME.Filter(f'gRandom->Rndm() < {fracToKeep}')  # todo: remove
        kStarME = rdfME.AsNumpy(['kStar'])['kStar']
        nPairsME = len(kStarME)

        kdeME = TKDE(nPairsME, kStarME, 0.0, kStarMax, kdeSetupStr)

        kdeEstimates[f'{dId}_centr'] = kdeME.GetFunction(1000, 0, kStarMax).Clone()
        kdeEstimates[f'{dId}_centr'].SetTitle(" ;k* (MeV/#it{c});Density")
        kdeEstimates[f'{dId}_centr'].SetName(f'f{dId}_centr')
        kdeEstimates[f'{dId}_centr'].Write()

        kdeEstimates[f'{dId}_upper'] = kdeME.GetUpperFunction(0.683, 1000, 0, kStarMax)
        kdeEstimates[f'{dId}_upper'].SetTitle(" ;k* (MeV/#it{c});Density")
        kdeEstimates[f'{dId}_upper'].SetName(f'f{dId}_upper')
        kdeEstimates[f'{dId}_upper'].Write()

        kdeEstimates[f'{dId}_lower'] = kdeME.GetLowerFunction(0.683, 1000, 0, kStarMax)
        kdeEstimates[f'{dId}_lower'].SetTitle(" ;k* (MeV/#it{c});Density")
        kdeEstimates[f'{dId}_lower'].SetName(f'f{dId}_lower')
        kdeEstimates[f'{dId}_lower'].Write()

        kdeEstimates[f'{dId}_relunc_upper'] = TF1(f'f{dId}_relunc_upper', cfRelUncUpperFcn, 0, kStarMax, 0)
        kdeEstimates[f'{dId}_relunc_upper'].SetTitle(" ;k* (MeV/#it{c});Relative error")
        kdeEstimates[f'{dId}_relunc_upper'].Write()

        kdeEstimates[f'{dId}_relunc_lower'] = TF1(f'f{dId}_relunc_lower', cfRelUncLowerFcn, 0, kStarMax, 0)
        kdeEstimates[f'{dId}_relunc_lower'].SetTitle(" ;k* (MeV/#it{c});Relative error")
        kdeEstimates[f'{dId}_relunc_lower'].Write()

        # dId = f'ME_{pId}_{rId}'

        # rdfME = RDataFrame(f'treeDpi_ME_{pId}_{rId}', inFileName)
        # dataME = rdfME.AsNumpy(['kStar'])

        # kStarME = dataME['kStar']
        # nPairsME = len(kStarME)

        # kdeME = TKDE(nPairsME, kStarME, 0.0, kStarMax, kdeSetupStr)

        # kdeEstimates[f'f{dId}centr'] = kdeME.GetFunction(1000, 0, kStarMax).Clone()
        # kdeEstimates[f'f{dId}centr'].SetTitle(" ;k* (MeV/#it{c});Density")
        # kdeEstimates[f'f{dId}centr'].SetName(f'ME{dId}centr')
        # kdeEstimates[f'f{dId}centr'].Write(f'fME{dId}centr')

        # kdeEstimates[f'f{dId}upper'] = kdeME.GetUpperFunction(0.683, 1000, 0, kStarMax)
        # kdeEstimates[f'f{dId}upper'].SetTitle(" ;k* (MeV/#it{c});Density")
        # kdeEstimates[f'f{dId}upper'].SetName(f'ME{dId}upper')
        # kdeEstimates[f'f{dId}upper'].Write(f'fME{dId}upper')

        # kdeEstimates[f'f{dId}lower'] = kdeME.GetLowerFunction(0.683, 1000, 0, kStarMax)
        # kdeEstimates[f'f{dId}lower'].SetTitle(" ;k* (MeV/#it{c});Density")
        # kdeEstimates[f'f{dId}lower'].SetName(f'ME{dId}lower')
        # kdeEstimates[f'f{dId}lower'].Write(f'fME{dId}lower')
        # weightsGF for mult
        hhSE = inFile.Get(f'hhDpi_Mult_kStar_SE_{pId}_{rId}')
        hhME = inFile.Get(f'hhDpi_Mult_kStar_ME_{pId}_{rId}')

        hSEProj = hhSE.ProjectionY('hSEProj', 1, hhSE.GetXaxis().GetNbins())
        hMEProj = hhME.ProjectionY('hMEProj', 1, hhME.GetXaxis().GetNbins())

        hSEProj.Write()
        hSEProj.Scale(1. / hSEProj.GetEntries())
        hSEProj.Write('hSEProjNorm')

        hMEProj.Write()
        hMEProj.Scale(1. / hMEProj.GetEntries())
        hMEProj.Write('hMEProjNorm')

        kdeMESlices = []
        kdeMESlicesUpperUnc = []
        kdeMESlicesLowerUnc = []

        kdeMERew = []

        print("Performing the multiplicity reweighting...")
        nBinsMult = hhSE.GetYaxis().GetNbins()
        # for iMultBin in range(10, 12):
        for iMultBin in range(nBinsMult):
            multMin = hhSE.GetYaxis().GetBinLowEdge(iMultBin)
            multMax = hhSE.GetYaxis().GetBinUpEdge(iMultBin)
            if iMultBin < 10 or iMultBin >= 11:
                slice = TF1("", "0", 0, kStarMax)
                slice.Write(f'fME_mult_{multMin:.0f}_{multMax:.0f}_{pId}_{rId}_centr')
                sliceUpperUnc = slice.Clone()
                sliceLowerUnc = slice.Clone()
                
                kdeMESlices.append(slice)
                kdeMESlicesUpperUnc.append(slice)
                kdeMESlicesLowerUnc.append(slice)
                continue

            print(f"Mult bin {iMultBin}: {multMin:.1f} < mult < {multMax:.1f} ")

            dataME_filter = rdfME.Filter(f'{multMin} < mult && mult < {multMax}').AsNumpy(['kStar'])['kStar']

            # if len(dataME_filter) > 0:
            if len(dataME_filter) > nPairsME / 100 and len(dataME_filter) > 1000:
                kdeME = TKDE(len(dataME_filter), dataME_filter, 0.0, kStarMax, kdeSetupStr)
                slice = kdeME.GetFunction(1000, 0, kStarMax).Clone()
                sliceUpperUnc = kdeME.GetUpperFunction(0.683, 1000, 0, kStarMax)
                sliceLowerUnc = kdeME.GetLowerFunction(0.683, 1000, 0, kStarMax)
                
                # print(sliceUpperUnc.Eval(2000))
                kdeMESlices.append(slice)
            else:
                slice = TF1("", "0", 0, kStarMax)
                sliceUpperUnc = slice.Clone()
                sliceLowerUnc = slice.Clone()
                kdeMESlices.append(slice)

            # print(sliceUpperUnc.Eval(2000))
            kdeMESlicesUpperUnc.append(sliceUpperUnc)
            kdeMESlicesLowerUnc.append(sliceLowerUnc)
            slice.Write(f'fME_mult_{multMin:.0f}_{multMax:.0f}_{pId}_{rId}_centr')
        # compute weights
        weights = ComputeMultWeights(hhSE, hhME)
        hWeights = TH1D(f'hWeights_{pId}_{rId}',
                        f'Weights {pId} {rId};mult;weight', nBinsMult, 0, nBinsMult)
        for iW, w in enumerate(weights):
            hWeights.SetBinContent(iW + 1, w)
        hWeights.Write()

        norm = sum([hWeights.GetBinContent(iSlice+1) * kdeMESlices[iSlice].Integral(0, kStarMax, 1e-4) for iSlice in range(nBinsMult)])
        def cfRewFcn(x, par):
            return sum([hWeights.GetBinContent(iSlice+1) * kdeMESlices[iSlice].Eval(x[0]) for iSlice in range(nBinsMult)])/norm
        
        cfMERew = TF1(f'fMERew_{pId}_{rId}_centr', cfRewFcn, 0, kStarMax)
        cfMERew.Write()

        def cfRewRelUncUpperFcn(x, par):
            [print(me.Eval(x[0])) for me in kdeMESlicesUpperUnc]
            sigmas = [kdeMESlicesUpperUnc[iSlice].Eval(x[0]) - kdeMESlices[iSlice].Eval(x[0]) for iSlice in range(nBinsMult)]
            print(sigmas)
            print(sigmas)
            error = np.sqrt(sum([w*w*s*s for (w, s) in zip(weights, sigmas)]))
            return error/cfMERew.Eval(x[0])/norm
        print(cfRewRelUncUpperFcn([2000], None))
        print(cfRewRelUncUpperFcn([2000], None))

        def cfRewRelUncLowerFcn(x, par):
            sigmas = [kdeMESlices[iSlice].Eval(x[0]) - kdeMESlicesLowerUnc[iSlice].Eval(x[0]) for iSlice in range(nBinsMult)]
            error = np.sqrt(sum([w*w*s*s for (w, s) in zip(weights, sigmas)]))            
            return error/cfMERew.Eval(x[0])/norm
        
        def cfRewUncUpperFcn(x, par):
            return cfMERew.Eval(x[0]) * (1 + cfRewRelUncUpperFcn(x, par))/norm
            
        def cfRewUncLowerFcn(x, par):
            return cfMERew.Eval(x[0]) * (1 - cfRewRelUncLowerFcn(x, par))/norm

        cfMERewRelUpperUnc = TF1(f'fMERew_{pId}_{rId}_relunc_upper', cfRewRelUncUpperFcn, 0, kStarMax)
        cfMERewRelLowerUnc = TF1(f'fMERew_{pId}_{rId}_relunc_lower', cfRewRelUncLowerFcn, 0, kStarMax)
        cfMERewRelUpperUnc.Write()
        cfMERewRelLowerUnc.Write()

        cfMERewUpperUnc = TF1(f'fMERew_{pId}_{rId}_upper', cfRewUncUpperFcn, 0, kStarMax)
        cfMERewLowerUnc = TF1(f'fMERew_{pId}_{rId}_lower', cfRewUncLowerFcn, 0, kStarMax)
        cfMERewUpperUnc.Write()
        cfMERewLowerUnc.Write()
        break
    oFile.Close()
    print('Densities were saved in:', oFileName)


if __name__ == '__main__':
    # parser = argparse.ArgumentParser(description='Arguments to pass')
    # parser.add_argument('inFileName', metavar='text',
    #                     help='AnalysisResults.root file')
    # parser.add_argument('oFileName', metavar='text',
    #                     help='AnalysisResults.root file')

    # args = parser.parse_args()

    # EstimateDensity(inFileName=args.inFileName, oFileName=args.oFileName)

    # kernelTypes = ['Gaussian', 'Epanechnikov', 'Biweight', 'CosineArch']
    # iterTypes = ['Adaptive', 'Fixed']
    # mirrorTypes = ['NoMirror', 'MirrorBoth', 'MirrorAsymBoth']
    # binningTypes = ['Unbinned', 'RelaxedBinning', 'ForcedBinning']
    # EstimateDensity()

    kernelSetup = {}

    # default
    kernelSetup['kernel'] = 'Gaussian'
    kernelSetup['iter'] = 'Adaptive'
    kernelSetup['mirror'] = 'NoMirror'
    kernelSetup['binning'] = 'RelaxedBinning'
    print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    EstimateDensity(kernelSetup)

    # # variation on kernel
    # kernelSetup['kernel'] = 'Epanechnikov'
    # kernelSetup['iter'] = 'Adaptive'
    # kernelSetup['mirror'] = 'NoMirror'
    # kernelSetup['binning'] = 'RelaxedBinning'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)
    # kernelSetup['kernel'] = 'Biweight'
    # kernelSetup['iter'] = 'Adaptive'
    # kernelSetup['mirror'] = 'NoMirror'
    # kernelSetup['binning'] = 'RelaxedBinning'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)
    # kernelSetup['kernel'] = 'CosineArch'
    # kernelSetup['iter'] = 'Adaptive'
    # kernelSetup['mirror'] = 'NoMirror'
    # kernelSetup['binning'] = 'RelaxedBinning'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)

    # variation on iter
    # kernelSetup['kernel'] = 'Gaussian'
    # kernelSetup['iter'] = 'Fixed'
    # kernelSetup['mirror'] = 'NoMirror'
    # kernelSetup['binning'] = 'RelaxedBinning'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)

    # # # variation on mirror
    # kernelSetup['kernel'] = 'Gaussian'
    # kernelSetup['iter'] = 'Fixed'
    # kernelSetup['mirror'] = 'MirrorBoth'
    # kernelSetup['binning'] = 'RelaxedBinning'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)
    # kernelSetup['kernel'] = 'Gaussian'
    # kernelSetup['iter'] = 'Fixed'
    # kernelSetup['mirror'] = 'MirrorAsymBoth'
    # kernelSetup['binning'] = 'RelaxedBinning'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)

    # # # variation on binning
    # kernelSetup['kernel'] = 'Gaussian'
    # kernelSetup['iter'] = 'Fixed'
    # kernelSetup['mirror'] = 'MirrorBoth'
    # kernelSetup['binning'] = 'Unbinned'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)
