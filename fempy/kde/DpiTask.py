import luigi
import sys
import os

from MergePairs import MergePairs
# from ExtractGenCF import ExtractGenCF
from EstimateDensity import EstimateDensity

analysisBasePath = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/'


class CreateDatasetDirTask(luigi.Task):
    path = luigi.Parameter()

    def run(self):
        os.makedirs(self.path, exist_ok=True)

    def output(self):
        return luigi.LocalTarget(self.path)


class CreateDatasetTask(luigi.Task):
    inFileName = luigi.Parameter()
    oFileName = luigi.Parameter()
    scaleFactor = luigi.Parameter()
    isMCstr = luigi.Parameter(default='false')

    def run(self):
        macroName = "/home/daniel/alice/femtokde/CreateFakeDataSet.C"
        command = f'root -q \'{macroName}("{self.inFileName}", "{self.oFileName}", {self.scaleFactor}, {self.isMCstr})\''
        print('\033[31mRunning:', command, '\033[0m')
        os.system(command)

    def output(self):
        return luigi.LocalTarget(self.oFileName)

    def requires(self):
        return [CreateDatasetDirTask(path=os.path.dirname(self.oFileName))]


class MergePairsTask(luigi.Task):
    # inFileName = luigi.Parameter()
    # oFileName = luigi.Parameter()
    # mergingCombinations = luigi.Parameter()

    inFileName = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults.root',
    oFileName = '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults_merged.root',
    mergingCombinations = 'pap'

    def run(self):
        MergePairs(inFileName=self.inFileName, oFileName=self.oFileName, combos=self.mergingCombinations)

    def output(self):
        return luigi.LocalTarget(self.oFileName)

    def requires(self):
        return [CreateDatasetTask(
            inFileName='/home/daniel/alice/CharmingAnalyses/DKDpi/oton_mctruth/data/AnalysisResults_all.root',
            oFileName='/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults.root',
            scaleFactor=0.001,
            isMCstr="false"
        )]


class EstimateDensityTask(luigi.Task):
    inFile = luigi.Parameter()
    oFile = luigi.Parameter()

    def run(self):

        EstimateDensity(inFileName=self.inFile, oFileName=self.oFile)

    def output(self):
        return luigi.LocalTarget(self.oFile)

    def requires(self):
        return [MergePairsTask(
                inFileName='/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults.root',
                oFileName='/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults_merged.root',
                mergingCombinations='pap'
                )]


class ComputeRawCFTask(luigi.Task):
    inFile = luigi.Parameter()
    oFile = luigi.Parameter()

    def run(self):
        macroName = "/home/daniel/alice/femtokde/ComputeRawCF.C"
        command = f'root \'{macroName}("{self.inFile}", "{self.oFile}")\''

        print('Running:', command)

        os.system(f'root {macroName}("{self.inFile}", "{self.oFile}")')

    def output(self):
        return luigi.LocalTarget(self.oFile)

    def requires(self):
        return [EstimateDensityTask(
                '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults_merged.root',
                '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/Densities.root',
                )]


class DpiTask(luigi.Task):

    def run(self):
        print('Done!')

    def requires(self):
        return [
            CreateDatasetTask(  # data
                inFileName='/home/daniel/alice/CharmingAnalyses/DKDpi/oton_mctruth/data/AnalysisResults_all.root',
                oFileName='/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults.root',
                scaleFactor=0.001,
                isMCstr="false"
            ),
            CreateDatasetTask(  # mcgp
                inFileName='/home/daniel/alice/CharmingAnalyses/DKDpi/oton_mctruth/mcgp/AnalysisResults_all.root',
                oFileName='/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/mcgp/AnalysisResults.root',
                scaleFactor=0.001,
                isMCstr="true"
            ),
            MergePairsTask(),
            # EstimateDensityTask(
            #     '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/data/data/AnalysisResults_merged.root',
            #     '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/Densities.root',
            # ),
            # ComputeRawCFTask(
            #     '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/Densities.root',
            #     '/home/daniel/alice/CharmingAnalyses/Dpi_HMpp13TeV_quick/RawCF.root',
            # )
        ]

    # def output(self):
    #     return luigi.LocalTarget('hello_world.txt')


if __name__ == '__main__':
    luigi.run()
