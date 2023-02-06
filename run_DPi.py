import luigi
import os

analysisBasePath = '/home/daniel/an/DPi'

datasets = {
    "10_tree_syst": {
        "data": [
            ("/data/DPi/10_tree_syst/data/AnalysisResults_4623.root", "4623"),
            ("/data/DPi/10_tree_syst/data/AnalysisResults_4624.root", "4624"),
            ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/005/AnalysisResults.root", "4625_005"),
            ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/001/AnalysisResults.root", "4625_001"),
            ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/003/AnalysisResults.root", "4625_003"),
            ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/004/AnalysisResults.root", "4625_004"),
            ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/002/AnalysisResults.root", "4625_002"),
            ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/006/AnalysisResults.root", "4625_006"),
            ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/007/AnalysisResults.root", "4625_007"),
            ("/data/DPi/10_tree_syst/data/AnalysisResults_4626.root", "4626"),
        ]
    }
}


class MakeDistrTask(luigi.Task):
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    kStarBW = luigi.FloatParameter()
    nJobs = luigi.IntParameter()
    targets = []
    treeDir = "HM_CharmFemto_DplusPion_Trees0"

    cfg = '/home/daniel/an/DPi/cfg_selection_nopc_centr.yml'
    oDir = f"{analysisBasePath}/{data_version}/distr"
    uniq ='false'
    macro = '/home/daniel/phsw/fempy/MakeDistr.C'

    def run(self):
        print("ok")
        # print(f'Running with {self.data_version}: Done!')

        # macro = '/home/daniel/phsw/fempy/MakeDistr.C'

        # commands = []
        # self.targets = []
        # uniq = 'false'

        # for (file, datasetSuffix) in datasets[self.data_version][self.sample]:
        #     suffix = f"{self.sample}_kStarBW{self.kStarBW}MeV_{datasetSuffix}"
        #     params = f'"{file}", "{cfg}", "{oDir}", "DPi", "{suffix}", {self.kStarBW}, "{self.treeDir}", {uniq}, {self.nJobs}'

        #     commands.append(f'root -l -b -q \'{macro}+({params})\'')
        #     self.targets.append(f'{oDir}/Distr_{suffix}.root')

        #     print('comma', commands[-1], self.targets[-1])

        #     os.system(commands[-1])

    def requires(self):
        suffixes =[]
        inFiles = []
        for (inFile, datasetSuffix) in datasets[self.data_version][self.sample]:
            suffixes.append(f'{self.sample}_kStarBW{self.kStarBW}MeV_{datasetSuffix}')
            inFiles.append(inFile)
        # print("lallalal", datasets[self.data_version][self.sample][])
        # print("lallalal", list(zip(datasets[self.data_version][self.sample], suffixes)))
        return [MakeOneDistrTask(inFile=f, suffix=s, data_version=self.data_version, sample=self.sample, kStarBW=self.kStarBW, nJobs=self.nJobs) for (f, s) in zip(inFiles, suffixes)]

    def output(self):
        oFiles = []
        for (_, datasetSuffix) in datasets[self.data_version][self.sample]:
            fn = f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.sample}_kStarBW{self.kStarBW}MeV_{datasetSuffix}.root'
            print("auauauauauau ", fn)
            oFiles.append(fn)
        return [luigi.LocalTarget(f) for f in oFiles]


class MakeOneDistrTask(luigi.Task):
    inFile = luigi.Parameter()
    suffix = luigi.Parameter()
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    kStarBW = luigi.FloatParameter()
    nJobs = luigi.IntParameter()
    targets = []
    treeDir = "HM_CharmFemto_DplusPion_Trees0"

    cfg = '/home/daniel/an/DPi/cfg_selection_nopc_centr.yml'
    uniq ='false'
    macro = '/home/daniel/phsw/fempy/MakeDistr.C'
    def run(self):
        oDir = f"{analysisBasePath}/{self.data_version}/distr"
        params = f'"{self.inFile}", "{self.cfg}", "{oDir}", "DPi", "{self.suffix}", {self.kStarBW}, "{self.treeDir}", {self.uniq}, {self.nJobs}'
        print(f'lalalalal root -l -b -q \'{self.macro}+({params})')
        os.system(f'root -l -b -q \'{self.macro}+({params})\'')

    def output(self):
        target = f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.suffix}.root'
        return luigi.LocalTarget(target)


#not tested
class MergeDistrTask(luigi.Task):
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    kStarBW = luigi.Parameter()

    oFile = f"~/an/DPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV.root"
    inFiles = f"~/an/DPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV_*.root"

    def run(self):
        os.system(f"hadd {self.oFile} {self.inFiles}")

    def output(self):
        return [luigi.LocalTarget(self.oFile)]

    # def requires(self):
    #     return super().requires()


# class ComputeCFTaskBroken(luigi.Task):
#     data_version = luigi.Parameter()
#     sample = luigi.Parameter()
#     kStarBW = luigi.Parameter()
#     sgnFitFunc = luigi.Parameter()
#     charmMassBW = luigi.FloatParameter()

#     inFile = f"~/an/DPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV.root"
#     oFile = f"~/an/DPi/{data_version}/cf/RawCF_{sample}_kStarBW{kStarBW}MeV_{sgnFitFunc}.root"

#     def run(self):
#         command = f"python3 -i ComputeDCFTrick.py --pair DPi {self.inFile} {analysisBasePath}/cf/   --suffix {self.sample}_{self.sgnFitFunc} --sgnFitFunc {self.sgnFitFunc} --charmMassBW={self.charmMassBW}"
#         print("running command: ", command)
#         os.system(command)

#     def output(self):
#         print("checking for file: " + self.oFile)
#         return [luigi.LocalTarget(self.oFile)]


# class DPiTaskBroken(luigi.Task):
#     data_version = luigi.Parameter()
#     kStarBW = luigi.FloatParameter()

#     # print(f"EEEE {data_version}   {kStarBW}", data_version)

#     def run(self):
#         print("Done")

#     def requires(self):
#         return [
#             # MakeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, nJobs=16),
#             # MergeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW),
#             # ComputeCFTaskBroken(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, sgnFitFunc="gaus", charmMassBW=0.2),
#         ]

class ComputeCFTask(luigi.Task):
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    kStarBW = luigi.FloatParameter()
    sgnFitFunc = luigi.Parameter()
    bkgFitFunc = luigi.Parameter()
    aliFitter = luigi.BoolParameter()
    charmMassBW = luigi.FloatParameter()

    # inFile = f"~/an/DPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV.root"
    # oFile = f"~/an/DPi/{data_version}/cf/RawCF_{sample}_kStarBW{kStarBW}MeV_{sgnFitFunc}.root"

    def run(self):
        inFile = f"~/an/DPi/{self.data_version}/distr/Distr_{self.sample}_kStarBW{self.kStarBW}MeV.root"

        parameters = '--pair DPi '
        parameters += f'{inFile} '
        parameters += f'{analysisBasePath}/{self.data_version}/cf/ '
        parameters += '--fitMass  '
        parameters += f'--sgnFitFunc {self.sgnFitFunc} '
        parameters += f'--bkgFitFunc {self.bkgFitFunc} '
        parameters += f'--charmMassBW {self.charmMassBW} '
        if self.aliFitter:
            parameters += f'--suffix {self.sample}_kStarBW{self.kStarBW}MeV_{self.sgnFitFunc}_charmMassBW{self.charmMassBW}_aliFitter '
            parameters += '--aliFitter '
        else:
            parameters += f'--suffix {self.sample}_kStarBW{self.kStarBW}MeV_{self.sgnFitFunc}_charmMassBW{self.charmMassBW} '

        command = f"python3 ComputeDstarCFTrick.py {parameters}"
        print("Running: ", command)
        os.system(command)


class DPiTask(luigi.Task):
    data_version = luigi.Parameter()
    kStarBW = luigi.FloatParameter()

    def run(self):
        print("ok")

    def requires(self):
        return [
            MakeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, nJobs=16),
        ]


if __name__ == '__main__':
    luigi.build([DPiTask(data_version='10_tree_syst', kStarBW=10)], local_scheduler=True)
