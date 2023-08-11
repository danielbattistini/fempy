from ROOT import TCanvas, TFile

eventIds = ["SE", "ME"]
pairIds = ["sc"]
regionIds = ["sgn", "sbl", "sbr"]


def CompareKernelSettings(inFileNames, oFileName):
    TCanvas('cVariations', 'cVariations', 600, 600)
    for inFileName in inFileNames:
        inFile = TFile(inFileName)

        densityIds = [f'{eId}_{pId}_{rId}' for eId in eventIds for pId in pairIds for rId in regionIds]

        cf = {}
        for dId in densityIds:
            cf['centr'] = inFile.Get(f'{dId}_central')
            cf['upper'] = inFile.Get(f'{dId}_upper')
            cf['lower'] = inFile.Get(f'{dId}_lower')

        cf['centr'].Draw()
        cf['upper'].Draw('same')
        cf['lower'].Draw('same')

    oFile = TFile(oFileName)


if __name__ == '__main__':
    kernelIds = ['Gaussian', 'Epanechnikov', 'Biweight', 'CosineArch']
    iterIds = ['Adaptive', 'Fixed']
    mirrorIds = ['NoMirror', 'MirrorBoth', 'MirrorAsymBoth']
    binningIds = ['Unbinned', 'RelaxedBinning', 'ForcedBinning']

    # # default
    # kernelSetup['kernel'] = 'Gaussian'
    # kernelSetup['iter'] = 'Adaptive'
    # kernelSetup['mirror'] = 'NoMirror'
    # kernelSetup['binning'] = 'RelaxedBinning'
    # print('Performing kde estimation using the setting:', list(kernelSetup.values()))
    # EstimateDensity(kernelSetup)
    kernelVariationList = [
        f'Densities_kern-{k}_iter-Adaptive_mirror-NoMirror_binning-RelaxedBinning.root' for k in kernelIds]
    CompareKernelSettings(kernelVariationList, "test_variations.root")
