import luigi
import os

analysisBasePath = '/home/daniel/an/DstarPi'

datasets = {
    "18_fixmix": {
        "data": [
            "/data/DstarPi/18_fixmix/data/AnalysisResults_4654.root",
            "/data/DstarPi/18_fixmix/data/AnalysisResults_4655.root",
            "/data/DstarPi/18_fixmix/data/AnalysisResults_4656.root",
            "/data/DstarPi/18_fixmix/data/AnalysisResults_4657.root",
        ]
    }
}


class MakeOneDistrTask(luigi.Task):
    command = luigi.Parameter()
    target = luigi.Parameter()

    def run(self):
        os.system(self.command)

    def output(self):
        return luigi.LocalTarget(self.target)


class MakeDistrTask(luigi.Task):
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    kStarBW = luigi.FloatParameter()
    nJobs = luigi.IntParameter()
    targets = []

    def run(self):
        print(f'Running with {self.data_version}: Done!')

        macro = '/home/daniel/phsw/fempy/MakeDistr.C'
        cfg = '/home/daniel/an/DstarPi/cfg_selection_nopc_syst.yml'
        treeDir = "HM_CharmFemto_DstarPion_Trees0"

        oDir = f"{analysisBasePath}/{self.data_version}/distr"
        commands = []
        self.targets = []

        for file in datasets[self.data_version][self.sample]:
            rn = file[-9:-5]
            suffix = f"{self.sample}_kStarBW{self.kStarBW}MeV_{rn}"
            params = f'"{file}", "{cfg}", "{oDir}", "DstarPi", "{suffix}", {self.kStarBW}, "{treeDir}", {self.nJobs}'

            commands.append(f'root -l -b -q \'{macro}+({params})\'')
            self.targets.append(f'{oDir}/Distr_{suffix}.root')

            print('comma', commands[-1], self.targets[-1])

            os.system(commands[-1])

    def output(self):
        return [luigi.LocalTarget("/home/daniel/ciccioc")]


#not tested
class MergeDistrTask(luigi.Task):
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    kStarBW = luigi.Parameter()

    oFile = f"~/an/DstarPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV.root"
    inFiles = f"~/an/DstarPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV_*.root"

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

#     inFile = f"~/an/DstarPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV.root"
#     oFile = f"~/an/DstarPi/{data_version}/cf/RawCF_{sample}_kStarBW{kStarBW}MeV_{sgnFitFunc}.root"

#     def run(self):
#         command = f"python3 -i ComputeDstarCFTrick.py --pair DstarPi {self.inFile} {analysisBasePath}/cf/   --suffix {self.sample}_{self.sgnFitFunc} --sgnFitFunc {self.sgnFitFunc} --charmMassBW={self.charmMassBW}"
#         print("running command: ", command)
#         os.system(command)

#     def output(self):
#         print("checking for file: " + self.oFile)
#         return [luigi.LocalTarget(self.oFile)]


# class DstarPiTaskBroken(luigi.Task):
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
    data_version=luigi.Parameter()
    sample=luigi.Parameter()
    kStarBW=luigi.FloatParameter()
    sgnFitFunc=luigi.Parameter()
    bkgFitFunc=luigi.Parameter()
    aliFitter=luigi.BoolParameter()
    charmMassBW=luigi.FloatParameter()

    # inFile = f"~/an/DstarPi/{data_version}/distr/Distr_{sample}_kStarBW{kStarBW}MeV.root"
    # oFile = f"~/an/DstarPi/{data_version}/cf/RawCF_{sample}_kStarBW{kStarBW}MeV_{sgnFitFunc}.root"


    def run(self):
        inFile = f"~/an/DstarPi/{self.data_version}/distr/Distr_{self.sample}_kStarBW{self.kStarBW}MeV.root"

        parameters = '--pair DstarPi '
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



class DstarPiTask(luigi.Task):
    data_version = luigi.Parameter()
    kStarBW = luigi.FloatParameter()

    def run(self):
        print("ok")

    def requires(self):
        return [
            # MakeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, nJobs=16),
            # MergeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW),
            # ComputeCFTaskBroken(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, sgnFitFunc="gaus", charmMassBW=0.2),
            # ComputeCFTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, sgnFitFunc="gaus", charmMassBW=0.2),
            ComputeCFTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, sgnFitFunc="gaus", bkgFitFunc="powex", charmMassBW=0.5, aliFitter=True),
            ComputeCFTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, sgnFitFunc="gaus", bkgFitFunc="powex", charmMassBW=0.5, aliFitter=False),
            ComputeCFTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, sgnFitFunc="hat", bkgFitFunc="powex", charmMassBW=0.5, aliFitter=False),
        ]

if __name__ == '__main__':
    luigi.build([DstarPiTask(data_version='18_fixmix', kStarBW=10)], local_scheduler=True)
    # luigi.build([DstarPiTaskBroken(data_version='18_fixmix', kStarBW=10)], local_scheduler=True)
    # luigi.build([DstarPiTaskBroken(data_version="18_fixmix", kStarBW=15)], local_scheduler=True)
