'''Macro to merge same/opposite charge SE/ME distributions.'''

from posixpath import ismount
from ROOT import TFile, TList, TTree
from sys import exit


def MergePairs():
    # merge data
    isMC = False
    combos = 'pap'

    if isMC:
        inFileName = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/mcgp/AnalysisResults.root'
        oFileName = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/mcgp/AnalysisResults_merged.root'
    else:
        inFileName = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults.root'
        oFileName = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults_merged.root'

    if combos == 'pap':  # pair anti-pair
        mergeComb = {
            'sc_sgn': ['sgn_Particle0_Particle2', 'sgn_Particle1_Particle3'],
            'sc_sbl': ['sbl_Particle0_Particle2', 'sbl_Particle1_Particle3'],
            'sc_sbr': ['sbr_Particle0_Particle2', 'sbr_Particle1_Particle3'],
        }

    inFile = TFile(inFileName)
    oFile = TFile(oFileName, 'recreate')

    inFile.ls()

    print('Merging trees into', oFileName, '...')
    for mergingId, mergingList in mergeComb.items():
        # merge trees
        for eventType in ['SE', 'ME']:
            cue = TList()
            for name in mergingList:
                cue.Add(inFile.Get(f'treeDpi_{eventType}_{name}'))
            mergedTree = TTree.MergeTrees(cue)
            mergedTree.Write(f'treeDpi_{eventType}_{mergingId}')
            print('TTree', f'treeDpi_{eventType}_{mergingId}', 'was written.')

            hToSum = []
            for name in mergingList:
                print(f'--- hhDpi_Mult_kStar_{eventType}_{name}')
                hToSum.append(inFile.Get(f'hhDpi_Mult_kStar_{eventType}_{name}'))

            hSummed = hToSum[0]
            for h in hToSum[1:]:
                hSummed.Add(h)
            hSummed.Write(f'hhDpi_Mult_kStar_{eventType}_{mergingId}')

    # merge MC trees
    if isMC:
        for eventType in ['SE', 'ME']:
            for mcType in ['Gen', 'Rec']:
                cue = TList()

                for name in mergeComb['sc_sgn']:
                    cue.Add(inFile.Get(f'treeDpi_{eventType}_MC{mcType}_{name}'))
                mergedTree = TTree.MergeTrees(cue)
                mergedTree.Write(f'treeDpi_{eventType}_MC{mcType}_sc_sgn')
                print('TTree', f'treeDpi_{eventType}_MC{mcType}_sc_sgn', 'was written.')

    oFile.Close()


MergePairs()
