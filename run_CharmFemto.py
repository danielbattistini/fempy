import argparse
import luigi
import sys
import os

analysisBasePath = ''
datasets = {}
pair = ''

def get_labels(version, sample):
    return [label for (_, label) in datasets[version][sample]]
        

class CharmFemtoTask(luigi.Task):
    data_version = luigi.Parameter()
    cleaner = luigi.Parameter()
    kStarBW = luigi.Parameter()

class MakeDirTreeTask(luigi.Task):
    data_version = luigi.Parameter()

    def dir_tree(self):
        return [
            os.path.join(analysisBasePath, self.data_version),
            os.path.join(analysisBasePath, self.data_version,  'distr'),
        ]

    def run(self):
        for dir in self.dir_tree():
            print("creating directory:", dir)
            os.mkdir(dir)

    def output(self):
        return [luigi.LocalTarget(dir) for dir in self.dir_tree()]


class MakeDistrTask(CharmFemtoTask):
    sample = luigi.Parameter()
    kStarBW = luigi.FloatParameter()
    nJobs = luigi.IntParameter()
    doSyst = luigi.IntParameter()
    doOnlyFD = luigi.Parameter()
    origin = luigi.Parameter(default='')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.origin = self.origin if self.origin == '' else f'_{self.origin}'

        
    def run(self):
        print("Merging the distributions")
        oFile = f"~/an/{pair}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV{self.origin}.root"
        print("ofileeee: ", oFile)
        inFiles = " ".join([f"~/an/{pair}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{label}{self.origin}.root" for label in get_labels(self.data_version, self.sample)])
        print(f"hadd {oFile} {inFiles}")
        os.system(f"hadd {oFile} {inFiles}")

        # sum charge combos
        os.system(f"python3 SumPairs.py /home/daniel/an/{pair}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV{self.origin}.root")

    def requires(self):
        print(self.origin, "--->", end='  ')
        print(self.origin)
        os.system('root -l -b -q -e ".L /home/daniel/phsw/fempy/MakeDistr.C+"')
        
        # suffixes = []
        # inFiles = []
        # for (inFile, datasetSuffix) in datasets[self.data_version][self.sample]:
        #     suffixes.append(f'{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{datasetSuffix}{self.origin}')
        #     inFiles.append(inFile)
        # return [MakeOneDistrTask(
        #     inFile=f,
        #     suffix=s,
        #     data_version=self.data_version,
        #     sample=self.sample,
        #     kStarBW=self.kStarBW,
        #     nJobs=self.nJobs,
        #     cleaner=self.cleaner,
        #     doSyst=self.doSyst,
        #     origin=self.origin,
        #     uniq="true" if self.doSyst == 'false' and self.doOnlyFD == 'false' else 'false',
        #     nSelToKeep=50) for (f, s) in zip(inFiles, suffixes)]


        # suffixes = []
        # inFiles = []
        # for (inFile, datasetSuffix) in datasets[self.data_version][self.sample]:
        #     suffixes.append(f'{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{datasetSuffix}{self.origin}')
        #     inFiles.append(inFile)





        print("kkkkkkkk", [(f, 
            f'{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{lab}{self.origin}', 
            self.data_version, 
            self.sample, 
            self.kStarBW, 
            self.nJobs, 
            self.cleaner, 
            self.doSyst, 
            self.origin, 
            "true",
            ) for (f, lab) in datasets[self.data_version][self.sample]])

        if len(datasets[self.data_version][self.sample]) == 0:
            return None

        return [MakeOneDistrTask(
            inFile=f,
            suffix=f'{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{lab}{self.origin}',
            data_version=self.data_version,
            sample=self.sample,
            kStarBW=self.kStarBW,
            nJobs=self.nJobs,
            cleaner=self.cleaner,
            doSyst=self.doSyst,
            origin=self.origin,
            uniq="true",
            # uniq="true" if self.doSyst == 'false' and self.doOnlyFD == 'false' else 'false',
            nSelToKeep=50) for (f, lab) in datasets[self.data_version][self.sample]]



    def output(self):
        oFiles = []
        for (_, datasetSuffix) in datasets[self.data_version][self.sample]:
            fn = f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{datasetSuffix}{self.origin}.root'
            oFiles.append(fn)
        return [luigi.LocalTarget(f) for f in oFiles] + [luigi.LocalTarget(f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV{self.origin}.root')]


class MakeOneDistrTask(luigi.Task):
    inFile = luigi.Parameter()
    suffix = luigi.Parameter()
    uniq = luigi.Parameter()
    cleaner = luigi.Parameter()
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    kStarBW = luigi.FloatParameter()
    nJobs = luigi.IntParameter()
    origin = luigi.Parameter()
    doSyst = luigi.Parameter(default='false')

    nSelToKeep = luigi.IntParameter(default=0)
    doOnlyFD = luigi.Parameter(default='false')

    macro = '/home/daniel/phsw/fempy/MakeDistr.C'

    def run(self):
        if pair == 'DPi':
            if self.sample == 'data':
                treeDir = "HM_CharmFemto_DplusPion_Trees0"
            else:
                treeDir = "HM_CharmFemto_Dpion_Trees0"
        elif pair == 'DK':
            treeDir = "HM_CharmFemto_Dkaon_Trees0"
        elif pair == 'DstarPi':
            treeDir = "HM_CharmFemto_DstarPion_Trees0"

        isMC = 'true' if self.sample in ['mcgp', 'mchf'] else 'false'
        cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_syst{self.origin}.yml'
        oDir = f"{analysisBasePath}/{self.data_version}/distr"
        params = f'"{self.inFile}", "{cfg}", {isMC}, {self.nSelToKeep}, "{oDir}", "{pair}", "{self.suffix}", {self.kStarBW}, "{treeDir}", {self.uniq}, {self.doOnlyFD}, {self.doSyst}, {self.nJobs}'
        print(f'lalalalal root -l -b -q \'{self.macro}+({params})')
        os.system(f'root -l -b -q \'{self.macro}+({params})\'')

    def output(self):
        target = f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.suffix}.root'
        print("Looking for file: ", target)
        return luigi.LocalTarget(target)

    def requires(self):
        if pair == 'DPi':
            if self.sample == 'data':
                treeDir = "HM_CharmFemto_DplusPion_Trees0"
            else:
                treeDir = "HM_CharmFemto_Dpion_Trees0"
        elif pair == 'DK':
            treeDir = "HM_CharmFemto_Dkaon_Trees0"
        elif pair == 'DstarPi':
            treeDir = "HM_CharmFemto_DstarPion_Trees0"

        isMC = 'true' if self.sample in ['mcgp', 'mchf'] else 'false'
        cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_syst{self.origin}.yml'
        oDir = f"{analysisBasePath}/{self.data_version}/distr"
        params = f'"{self.inFile}", "{cfg}", {isMC}, {self.nSelToKeep}, "{oDir}", "{pair}", "{self.suffix}", {self.kStarBW}, "{treeDir}", {self.uniq}, {self.doOnlyFD}, {self.doSyst}, {self.nJobs}'
        print(f'lalalalal root -l -b -q \'{self.macro}+({params})')
        
        
        return MakeDirTreeTask(data_version=self.data_version)


class ComputeCFTask(luigi.Task):
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    cleaner = luigi.Parameter()
    kStarBW = luigi.FloatParameter()
    sgnFitFunc = luigi.Parameter()
    bkgFitFunc = luigi.Parameter()
    aliFitter = luigi.BoolParameter()
    charmMassBW = luigi.FloatParameter()
    lowKStarBW = luigi.FloatParameter()
    lowKStarMax = luigi.FloatParameter()

    def make_suffix(self):
        suffix = f'{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{self.sgnFitFunc}_charmMassBW{self.charmMassBW}MeV'
        if self.lowKStarBW > 0 and self.lowKStarMax > 0:
            suffix += f'_lowKStarBW{self.lowKStarBW:.1f}MeV_until{self.lowKStarMax:.0f}MeV'
        if self.aliFitter:
            suffix += '_aliFitter'
        return suffix

    def run(self):
        inFile = f"~/an/{pair}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV.root"

        parameters = f'--pair {pair} '
        parameters += f'{inFile} '
        parameters += f'{analysisBasePath}/{self.data_version}/cf/ '
        parameters += '--fitMass  '
        parameters += f'--sgnFitFunc {self.sgnFitFunc} '
        parameters += f'--bkgFitFunc {self.bkgFitFunc} '
        parameters += f'--charmMassBW {self.charmMassBW} '
        parameters += f'--suffix {self.make_suffix()} '
        if self.aliFitter:
            parameters += '--aliFitter '
        if self.lowKStarBW > 0 and self.lowKStarMax > 0:
            parameters += f'--lowKStarCharmMassBW {self.lowKStarBW:.1f} {self.lowKStarMax:.0f}'

        command = f"python3 ComputeDstarCFTrick.py {parameters}"
        print("Running: ", command)
        os.system(command)

    def output(self):
        return luigi.LocalTarget(f'{analysisBasePath}/{self.data_version}/cf/RawCF_{self.make_suffix()}.root')


class DPiTask(CharmFemtoTask):
    def requires(self):
        yield MakeDirTreeTask(data_version=self.data_version)

        yield MakeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true')
        
        yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true')
        yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromc')
        yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromb')
        yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromc')
        yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromb')

class DKTask(CharmFemtoTask):
    def requires(self):
        yield MakeDirTreeTask(data_version=self.data_version)
        yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true')

        yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromc')
        yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromb')
        yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromc')
        yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false', doOnlyFD='true', origin='fromb')


class DstarPiTask(CharmFemtoTask):
    def requires(self):
        yield MakeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, nJobs=2, cleaner=self.cleaner, doSyst='false')
        # yield ComputeCFTask(
        #     data_version=self.data_version,
        #     sample="data",
        #     cleaner=self.cleaner,
        #     kStarBW=self.kStarBW,
        #     sgnFitFunc="gaus",
        #     bkgFitFunc="powex",
        #     charmMassBW=0.2,
        #     lowKStarBW=0.4,
        #     lowKStarMax=50,
        #     aliFitter=True)
        # yield ComputeCFTask(
        #     data_version=self.data_version,
        #     sample="data",
        #     cleaner=self.cleaner,
        #     kStarBW=self.kStarBW,
        #     sgnFitFunc="gaus",
        #     bkgFitFunc="powex",
        #     charmMassBW=0.2,
        #     lowKStarBW=0.4,
        #     lowKStarMax=50,
        #     aliFitter=False)
        # yield ComputeCFTask(
        #     data_version=self.data_version,
        #     sample="data",
        #     cleaner=self.cleaner,
        #     kStarBW=self.kStarBW,
        #     sgnFitFunc="hat",
        #     bkgFitFunc="powex",
        #     charmMassBW=0.2,
        #     lowKStarBW=0.4,
        #     lowKStarMax=50,
        #     aliFitter=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("pair", choices=('DstarPi', 'DPi', 'DK'))
    args = parser.parse_args()

    pair = args.pair

    analysisBasePath = f'/home/daniel/an/{pair}'
    if pair == 'DstarPi':

        datasets = {
            "18_fixmix": {
                "data": [
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4654.root", "4654"),
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4655.root", "4655"),
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4656.root", "4656"),
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4657.root", "4657"),
                ]
            }
        }

        # luigi.build([
        #     DstarPiTask(data_version='18_fixmix', kStarBW=15, cleaner='nopc'),
        # ], local_scheduler=True, workers=16)
    
    elif pair == 'DPi':
        datasets = {
            "10_tree_syst": {
                "data": [
                    ("/data/DPi/10_tree_syst/data/AnalysisResults_4624.root", "4624"),
                    ("/data/DPi/10_tree_syst/data/AnalysisResults_4626.root", "4626"),
                    ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/005/AnalysisResults.root", "4625_005"),
                    ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/001/AnalysisResults.root", "4625_001"),
                    ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/003/AnalysisResults.root", "4625_003"),
                    ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/004/AnalysisResults.root", "4625_004"),
                    ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/002/AnalysisResults.root", "4625_002"),
                    ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/006/AnalysisResults.root", "4625_006"),
                    ("/data/DPi/10_tree_syst/data/merge_4625/Stage_3/007/AnalysisResults.root", "4625_007"),
                    ("/data/DPi/10_tree_syst/data/AnalysisResults_4623.root", "4623"),
                ],
                'mcgp': [
                    ("/data/DPi/10_tree_syst/mcgp/AnalysisResults_4148.root", "4148"),
                    ("/data/DPi/10_tree_syst/mcgp/AnalysisResults_4149.root", "4149"),
                    ("/data/DPi/10_tree_syst/mcgp/AnalysisResults_4147.root", "4147"),
                ],
                'mchf': [
                    ("/data/DPi/10_tree_syst/mchf/AnalysisResults_4150.root", "4150"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/015/AnalysisResults.root", "4151_015"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/024/AnalysisResults.root", "4151_024"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/040/AnalysisResults.root", "4151_040"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/037/AnalysisResults.root", "4151_037"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/038/AnalysisResults.root", "4151_038"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/027/AnalysisResults.root", "4151_027"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/020/AnalysisResults.root", "4151_020"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/032/AnalysisResults.root", "4151_032"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/025/AnalysisResults.root", "4151_025"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/008/AnalysisResults.root", "4151_008"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/039/AnalysisResults.root", "4151_039"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/033/AnalysisResults.root", "4151_033"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/010/AnalysisResults.root", "4151_010"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/005/AnalysisResults.root", "4151_005"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/031/AnalysisResults.root", "4151_031"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/001/AnalysisResults.root", "4151_001"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/044/AnalysisResults.root", "4151_044"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/009/AnalysisResults.root", "4151_009"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/018/AnalysisResults.root", "4151_018"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/003/AnalysisResults.root", "4151_003"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/011/AnalysisResults.root", "4151_011"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/016/AnalysisResults.root", "4151_016"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/004/AnalysisResults.root", "4151_004"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/013/AnalysisResults.root", "4151_013"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/036/AnalysisResults.root", "4151_036"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/042/AnalysisResults.root", "4151_042"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/034/AnalysisResults.root", "4151_034"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/045/AnalysisResults.root", "4151_045"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/029/AnalysisResults.root", "4151_029"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/021/AnalysisResults.root", "4151_021"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/002/AnalysisResults.root", "4151_002"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/041/AnalysisResults.root", "4151_041"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/014/AnalysisResults.root", "4151_014"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/019/AnalysisResults.root", "4151_019"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/006/AnalysisResults.root", "4151_006"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/017/AnalysisResults.root", "4151_017"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/028/AnalysisResults.root", "4151_028"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/043/AnalysisResults.root", "4151_043"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/030/AnalysisResults.root", "4151_030"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/022/AnalysisResults.root", "4151_022"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/012/AnalysisResults.root", "4151_012"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/023/AnalysisResults.root", "4151_023"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/026/AnalysisResults.root", "4151_026"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/035/AnalysisResults.root", "4151_035"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4151/Stage_1/007/AnalysisResults.root", "4151_007"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/015/AnalysisResults.root", "4152_015"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/024/AnalysisResults.root", "4152_024"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/020/AnalysisResults.root", "4152_020"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/008/AnalysisResults.root", "4152_008"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/010/AnalysisResults.root", "4152_010"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/005/AnalysisResults.root", "4152_005"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/001/AnalysisResults.root", "4152_001"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/009/AnalysisResults.root", "4152_009"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/018/AnalysisResults.root", "4152_018"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/003/AnalysisResults.root", "4152_003"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/011/AnalysisResults.root", "4152_011"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/016/AnalysisResults.root", "4152_016"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/004/AnalysisResults.root", "4152_004"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/013/AnalysisResults.root", "4152_013"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/021/AnalysisResults.root", "4152_021"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/002/AnalysisResults.root", "4152_002"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/014/AnalysisResults.root", "4152_014"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/019/AnalysisResults.root", "4152_019"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/006/AnalysisResults.root", "4152_006"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/017/AnalysisResults.root", "4152_017"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/022/AnalysisResults.root", "4152_022"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/012/AnalysisResults.root", "4152_012"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/023/AnalysisResults.root", "4152_023"),
                    ("/data/DPi/10_tree_syst/mchf/merge_4152/Stage_1/007/AnalysisResults.root", "4152_007"),
                ],
            },
            "11_mix150": {
                "data": [
                    ("/data/DPi/11_mix150/data/sgn/AnalysisResults_4690.root", "4690"),
                    ("/data/DPi/11_mix150/data/sgn/AnalysisResults_4691.root", "4691"),
                    ("/data/DPi/11_mix150/data/sgn/merge_4692/Stage_3/001/AnalysisResults.root", "4692_001"),
                    ("/data/DPi/11_mix150/data/sgn/merge_4692/Stage_3/002/AnalysisResults.root", "4692_002"),
                    ("/data/DPi/11_mix150/data/sgn/merge_4692/Stage_3/003/AnalysisResults.root", "4692_003"),
                    ("/data/DPi/11_mix150/data/sgn/merge_4692/Stage_3/004/AnalysisResults.root", "4692_004"),
                    ("/data/DPi/11_mix150/data/sgn/merge_4692/Stage_3/005/AnalysisResults.root", "4692_005"),
                    ("/data/DPi/11_mix150/data/sgn/merge_4692/Stage_3/006/AnalysisResults.root", "4692_006"),
                    ("/data/DPi/11_mix150/data/sgn/merge_4692/Stage_3/007/AnalysisResults.root", "4692_007"),
                    ("/data/DPi/11_mix150/data/sgn/AnalysisResults_4693.root", "4693"),
                ],
                'mcgp': [ ],
                'mchf': [ ],
            }
        }


        luigi.build([
            DPiTask(data_version='10_tree_syst', kStarBW=15, cleaner='nopc'),
            DPiTask(data_version='10_tree_syst', kStarBW=15, cleaner='newpc'),
            DPiTask(data_version='10_tree_syst', kStarBW=15, cleaner='oldpc'),
            DPiTask(data_version='11_mix150', kStarBW=15, cleaner='nopc'),
            DPiTask(data_version='11_mix150', kStarBW=15, cleaner='newpc'),
            DPiTask(data_version='11_mix150', kStarBW=15, cleaner='oldpc'),
            # DPiTask(data_version='11_mix150', kStarBW=15, cleaner='newpc'),
            # DPiTask(data_version='11_mix150', kStarBW=15, cleaner='oldpc'),
        ], workers=10)

    elif pair == "DK":
        datasets = {
            "10_tree_syst": {
                "data": [
                    ("/data/DK/10_tree_syst/data/AnalysisResults_4620.root", "4620"),
                    ("/data/DK/10_tree_syst/data/AnalysisResults_4622.root", "4622"),
                    ("/data/DK/10_tree_syst/data/AnalysisResults_4621.root", "4621"),
                    ("/data/DK/10_tree_syst/data/AnalysisResults_4619.root", "4619"),
                ],
                'mcgp': [
                    ("/data/DK/10_tree_syst/mcgp/AnalysisResults_4153.root", "4153"),
                    ("/data/DK/10_tree_syst/mcgp/AnalysisResults_4154.root", "4154"),
                    ("/data/DK/10_tree_syst/mcgp/AnalysisResults_4155.root", "4155"),
                ],
                'mchf': [
                    ("/data/DK/10_tree_syst/mchf/AnalysisResults_4158.root", "4158"),
                    ("/data/DK/10_tree_syst/mchf/AnalysisResults_4156.root", "4156"),
                    ("/data/DK/10_tree_syst/mchf/AnalysisResults_4157.root", "4157"),
                ],
            }
        }
        luigi.build([
            DKTask(data_version='10_tree_syst', kStarBW=15, cleaner='nopc'),
            DKTask(data_version='10_tree_syst', kStarBW=15, cleaner='newpc'),
            DKTask(data_version='10_tree_syst', kStarBW=15, cleaner='oldpc'),
        ], workers=10)











# class CharmFemto(luigi.Task):
#     pair = luigi.Parameter(pair)

#     def requires(self):
#         if pair == "DPi":
#             yield DPiTask(data_version='10_tree_syst', kStarBW=15, cleaner='nopc')
