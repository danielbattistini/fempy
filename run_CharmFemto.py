import argparse
import luigi

from alive_progress import alive_bar

import sys
import os
import fempy
# from CompareGraphs import CompareGraphs

luigi.interface.core.log_level="WARNING"


def MakeSuffix(**kwargs):
    suffix = ''
    if kwargs.get('sample') is not None:
        suffix += '_' + kwargs['sample']
    if kwargs.get('cleaner') is not None:
        suffix += '_' + kwargs['cleaner']
    if kwargs.get('kStarBW') is not None:
        suffix += '_kStarBW' + str(kwargs['kStarBW']) + 'MeV'
    if kwargs.get('sgnFitFunc') is not None:
        suffix += '_' + kwargs['sgnFitFunc']
    if kwargs.get('charmMassBW') is not None:
        suffix += f'_charmMassBW{kwargs["charmMassBW"]:.1f}MeV'
    if kwargs.get('lowKStarCharmMassBW') is not None:
        suffix += f'_lowKStarCharmMassBW{kwargs["lowKStarCharmMassBW"]:.1f}MeV'
    if kwargs.get('until') is not None:
        suffix += f'_until{kwargs["until"]:.1f}MeV'
    if kwargs.get('origin') is not None:
        suffix += '_' + kwargs['origin']
    if kwargs.get('fitRange') is not None:
        suffix += f'_fitRange-{kwargs["fitRange"][0]}-{kwargs["fitRange"][1]}MeV'
    if kwargs.get('bl') is not None:
        suffix += '_' + kwargs['bl']
    if kwargs.get('bootstrap') is not None and kwargs.get('bootstrap') > 0:
        suffix += '_bootstrap' + str(kwargs['bootstrap'])
    if kwargs.get('syst') is not None:
        suffix += '_syst'
    return suffix


analysisBasePath = ''
datasets = {}
pair = ''

def Chunkyfy(array, nElemPerChunk=50):
    nChunks = len(array)%nElemPerChunk
    return [array[iChunk*nElemPerChunk:nElemPerChunk * (iChunk +1)] for iChunk in range(nChunks)]


class CharmFemtoTask(luigi.Task):
    data_version = luigi.Parameter()
    cleaner = luigi.Parameter()
    kStarBW = luigi.IntParameter(default=-1)


def get_labels(version, sample):
    return [label for (_, label) in datasets[version][sample]]


class MakeDirTreeTask(luigi.Task):
    data_version = luigi.Parameter()

    def dir_tree(self):
        return [
            os.path.join(analysisBasePath, self.data_version),
            os.path.join(analysisBasePath, self.data_version,  'distr'),
            os.path.join(analysisBasePath, self.data_version,  'cf'),
            os.path.join(analysisBasePath, self.data_version,  'comparisons'),
        ]

    def run(self):
        for folder in self.dir_tree():
            print("creating directory:", folder)
            if not os.path.exists(folder):
                os.mkdir(folder)

    def output(self):
        return [luigi.LocalTarget(dir) for dir in self.dir_tree()]


class MakeDistrTask(CharmFemtoTask):
    sample = luigi.Parameter()
    kStarBW = luigi.FloatParameter(default=-1)
    nJobs = luigi.IntParameter()
    doSyst = luigi.IntParameter(default='false')
    doQA = luigi.Parameter(default='false')
    origin = luigi.Parameter(default='')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.origin = self.origin if self.origin == '' else f'_{self.origin}'

    def run(self):
        print("Merging the distributions")
        
        oFile = f"~/an/{pair}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}"
        if self.kStarBW > 0:
            oFile +=f"_kStarBW{self.kStarBW}MeV"
        oFile += f"{self.origin}.root"

        inFileNames = []
        for label in get_labels(self.data_version, self.sample):
            inFileName = f"/home/daniel/an/{pair}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}"
            if self.kStarBW > 0:
                inFileName += f"_kStarBW{self.kStarBW}MeV"
            inFileName += f"_{label}{self.origin}.root"
            inFileNames.append(inFileName)
        inFiles = " ".join(inFileNames)
        if inFiles != '':
            fempy.info(f"hadd {oFile} {inFiles}")

            os.system(f"cp {inFileNames[0]} {oFile}")

            with alive_bar(len(inFileNames[1:]), force_tty=True) as bar:
                for chunk in Chunkyfy(inFileNames[1:]): # merge one by one to avoid exploding the RAM
                    print('Merging:', chunk)
                    command = f'hadd -v 1 -a -j 16 {oFile} {" ".join(chunk)}'
                    print("piririr", command)
                    os.system(command)
                    bar()
            # sum charge combos
            fempy.info(f"python3 SumPairs.py {oFile}")
            os.system(f"python3 SumPairs.py {oFile}")

    def requires(self):
        # print("Compiling MakeDistr.C macro...")
        # os.system('root -l -b -q -e ".L /home/daniel/phsw/fempy/MakeDistr.C+"')

        # print("Executing projection with the following parameters:\n",
        #       [(f,
        #         f'{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV_{lab}{self.origin}',
        #         self.data_version,
        #         self.sample,
        #         self.kStarBW,
        #         self.nJobs,
        #         self.cleaner,
        #         self.doSyst,
        #         self.origin,
        #         "true",
        #         ) for (f, lab) in datasets[self.data_version][self.sample]])

        if len(datasets[self.data_version][self.sample]) == 0:
            return None

        suffixes = []
        files = []
        for (f, lab) in datasets[self.data_version][self.sample]:
            suffix=f'{self.sample}_{self.cleaner}'
            if self.kStarBW > 0:
                suffix += f'_kStarBW{self.kStarBW}MeV'
            suffix += f'_{lab}{self.origin}'
            suffixes.append(suffix)
            files.append(f)
        return [MakeOneDistrTask(
            inFile=file,
            suffix=suffix,
            data_version=self.data_version,
            sample=self.sample,
            kStarBW=self.kStarBW,
            nJobs=self.nJobs,
            cleaner=self.cleaner,
            doSyst=self.doSyst,
            origin=self.origin,
            # doOnlyFD=self.doOnlyFD,
            doQA=self.doQA,
            # uniq="true" if self.doSyst == 'false' and self.doOnlyFD == 'false' else 'false',
            uniq='false',
            doKStarSlices='true' if 'Dstar' in analysisBasePath else 'false',
            nSelToKeep=20) for (file, suffix) in zip(files, suffixes)]

    def output(self):
        oFiles = []
        for (_, datasetSuffix) in datasets[self.data_version][self.sample]:
            fn = f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}'
            if self.kStarBW > 0:
                fn += f'_kStarBW{self.kStarBW}MeV'
            fn += f'_{datasetSuffix}{self.origin}.root'
            oFiles.append(fn)

        fnMerged = f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}'
        if self.kStarBW > 0:
            fnMerged += f'_kStarBW{self.kStarBW}MeV'
        fnMerged += f'{self.origin}.root'
        oFiles.append(fn)
        # print("ffff ", fnMerged)
        return [luigi.LocalTarget(f) for f in oFiles] + [luigi.LocalTarget(fnMerged)]


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
    doKStarSlices = luigi.Parameter(default='false')

    nSelToKeep = luigi.IntParameter(default=0)
    # doOnlyFD = luigi.Parameter(default='false')
    doQA=luigi.Parameter(default='false')

    macro = '/home/daniel/phsw/fempy/MakeDistr.C'

    def run(self):
        if pair == 'DPi':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DpionMCtruthReco_Trees0"
            elif '/sbl/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBLeft_DplusPion_Trees0"
            elif '/sbr/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBRight_DplusPion_Trees0"
            elif self.sample == 'data':
                treeDir = "HM_CharmFemto_DplusPion_Trees0"
            else:
                treeDir = "HM_CharmFemto_Dpion_Trees0"
        elif pair == 'DPiClosure':
            # only data
            treeDir = "HM_CharmFemto_Dpion_Trees0"
        elif pair == 'DK':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DkaonMCtruthReco_Trees0"
            elif '/sbl/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBLeft_Dkaon_Trees0"
            elif '/sbr/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBRight_Dkaon_Trees0"
            else:
                treeDir = "HM_CharmFemto_Dkaon_Trees0"

        elif pair == 'DstarPi':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DstarPionMCtruthReco_Trees0"
            else:
                treeDir = "HM_CharmFemto_DstarPion_Trees0"
        elif pair == 'DstarK':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DstarKaonMCtruthReco_Trees0"
            if 'data' in str(self.sample) or 'mcgp' in str(self.sample) or 'mchf' in str(self.sample):
                treeDir = "HM_CharmFemto_DstarKaon_Trees0"
            else:
                fempy.error("not implemented.")

        isMC = 'true' if self.sample in ['mcgp', 'mchf', 'mcgptruthreco', 'mcgptruthgen', 'mchftruthreco', 'mchftruthgen'] else 'false'
        if 'truthreco' in str(self.sample):
            cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_mctruth{self.origin}.yml'
        else:
            cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_syst{self.origin}.yml'

        if not os.path.isfile(cfg):
            fempy.error(f"the config file {cfg} doesn't exist.")

        oDir = f"{analysisBasePath}/{self.data_version}/distr"
        params = (
            f'"{self.inFile}", '
            f'"{cfg}", '
            f'{isMC}, '
            f'"{oDir}", '
            f'"{pair}", '
            f'"{self.suffix}", '
            f'"{treeDir}", '
            f'{self.uniq}, '
            f'{self.doKStarSlices}, '
            f'{self.kStarBW}, '
            f'{self.doQA}, '
            f'{self.doSyst}, '
            f'{self.nSelToKeep}, '
            f'{self.nJobs}'
        )
        fempy.info(f'Running one distr root -l -b -q \'{self.macro}+({params})')
        os.system(f'root -l -b -q \'{self.macro}+({params})\'')

    def output(self):
        target = f'{analysisBasePath}/{self.data_version}/distr/Distr_{self.suffix}.root'
        # print("Looking for file: ", target)
        return luigi.LocalTarget(target)

    def requires(self):
        if pair == 'DPi':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DpionMCtruthReco_Trees0"
            elif '/sbl/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBLeft_DplusPion_Trees0"
            elif '/sbr/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBRight_DplusPion_Trees0"
            elif self.sample == 'data':
                treeDir = "HM_CharmFemto_DplusPion_Trees0"
            else:
                treeDir = "HM_CharmFemto_Dpion_Trees0"
        elif pair == 'DPiClosure':
            # only data
            treeDir = "HM_CharmFemto_Dpion_Trees0"
        elif pair == 'DK':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DkaonMCtruthReco_Trees0"
            elif '/sbl/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBLeft_Dkaon_Trees0"
            elif '/sbr/' in str(self.inFile):
                treeDir = "HM_CharmFemto_SBRight_Dkaon_Trees0"
            else:
                treeDir = "HM_CharmFemto_Dkaon_Trees0"

        elif pair == 'DstarPi':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DstarPionMCtruthReco_Trees0"
            else:
                treeDir = "HM_CharmFemto_DstarPion_Trees0"
        elif pair == 'DstarK':
            if 'truthreco' in str(self.sample):
                treeDir = "HM_CharmFemto_DstarKaonMCtruthReco_Trees0"
            if 'data' in str(self.sample) or 'mcgp' in str(self.sample) or 'mchf' in str(self.sample):
                treeDir = "HM_CharmFemto_DstarKaon_Trees0"
            else:
                fempy.error("not implemented.")

        isMC = 'true' if self.sample in ['mcgp', 'mchf', 'mcgptruthreco', 'mcgptruthgen', 'mchftruthreco', 'mchftruthgen'] else 'false'
        if 'truthreco' in str(self.sample):
            cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_mctruth{self.origin}.yml'
        else:
            cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_syst{self.origin}.yml'

        if not os.path.isfile(cfg):
            fempy.error(f"the config file {cfg} doesn't exist.")
        
        oDir = f"{analysisBasePath}/{self.data_version}/distr"
        
        params = (
            f'"{self.inFile}", '
            f'"{cfg}", '
            f'{isMC}, '
            f'"{oDir}", '
            f'"{pair}", '
            f'"{self.suffix}", '
            f'"{treeDir}", '
            f'{self.uniq}, '
            f'{self.doKStarSlices}, '
            f'{self.kStarBW}, '
            f'{self.doQA}, '
            f'{self.doSyst}, '
            f'{self.nSelToKeep}, '
            f'{self.nJobs}'
        )
        print(f'root -l -b -q \'{self.macro}+({params})\'')
        # if pair == 'DPi':
        #     if self.sample == 'data':
        #         treeDir = "HM_CharmFemto_DplusPion_Trees0"
        #     else:
        #         treeDir = "HM_CharmFemto_Dpion_Trees0"
        # elif pair == 'DK':
        #     treeDir = "HM_CharmFemto_Dkaon_Trees0"
        # elif pair == 'DstarPi':
        #     treeDir = "HM_CharmFemto_DstarPion_Trees0"

        # isMC = 'true' if self.sample in ['mcgp', 'mchf', 'mcgptruthreco', 'mcgptruthgen', 'mchftruthreco', 'mchftruthgen'] else 'false'
        # if self.sample == 'mcgptruthreco':
        #     cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_mctruth{self.origin}.yml'
        # else:
        #     cfg = f'/home/daniel/an/{pair}/cfg_selection_{self.cleaner}_syst{self.origin}.yml'

        # oDir = f"{analysisBasePath}/{self.data_version}/distr"
        # params = f'"{self.inFile}", "{cfg}", {isMC}, {self.nSelToKeep}, "{oDir}", "{pair}", "{self.suffix}", {self.kStarBW}, "{treeDir}", {self.uniq}, {self.doOnlyFD}, {self.doSyst}, {self.nJobs}'
        # print(f'lalalalal root -l -b -q \'{self.macro}+({params})')

        return MakeDirTreeTask(data_version=self.data_version)


class ComputeCFTask(luigi.Task):
    data_version = luigi.Parameter()
    sample = luigi.Parameter()
    cleaner = luigi.Parameter()
    uniq = luigi.BoolParameter(default=False)
    fitMass = luigi.BoolParameter(default=False)
    kStarBW = luigi.FloatParameter()
    sgnFitFunc = luigi.OptionalParameter(default=None)
    bkgFitFunc = luigi.OptionalParameter(default=None)
    aliFitter = luigi.BoolParameter(default=False)
    charmMassBW = luigi.OptionalFloatParameter(default=None)  # MeV
    lowKStarCharmMassBW = luigi.FloatParameter(default=-1)
    lowKStarMax = luigi.FloatParameter(default=-1)
    origin = luigi.OptionalParameter(default=None)

    def make_suffix(self):
        suffix = f'{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV'
        if self.sgnFitFunc is not None:
            suffix += f'_{self.sgnFitFunc}'
        if self.charmMassBW is not None:
            suffix += f'_charmMassBW{self.charmMassBW}MeV'
        
        if self.lowKStarCharmMassBW > 0 and self.lowKStarMax > 0:
            suffix += f'_lowKStarCharmMassBW{self.lowKStarCharmMassBW:.1f}MeV_until{self.lowKStarMax:.0f}MeV'
        if self.aliFitter:
            suffix += '_aliFitter'
        if self.uniq:
            suffix += '_uniq'

        return suffix

    def run(self):
        inFile = f"~/an/{pair}/{self.data_version}/distr/Distr_{self.sample}_{self.cleaner}_kStarBW{self.kStarBW}MeV"
        if self.origin == 'true':
            inFile += '_true'
        inFile += '.root'

        parameters = f'--pair {pair} '
        parameters += f'{inFile} '
        parameters += f'{analysisBasePath}/{self.data_version}/cf/ '
        if self.sgnFitFunc is not None and self.bkgFitFunc is not None: 
            parameters += f'--sgnFitFunc {self.sgnFitFunc} '
            parameters += f'--bkgFitFunc {self.bkgFitFunc} '
        if self.charmMassBW is not None:
            parameters += f'--charmMassBW {self.charmMassBW} '
        parameters += f'--suffix {self.make_suffix()} '
        parameters += f'--kStarBW {self.kStarBW} '

        if self.aliFitter:
            parameters += '--aliFitter '
        if self.lowKStarCharmMassBW > 0 and self.lowKStarMax > 0:
            parameters += f'--lowKStarCharmMassBW {self.lowKStarCharmMassBW:.1f} {self.lowKStarMax:.0f} '

        if (self.uniq):
            parameters += '--uniq '
        if self.fitMass:
            parameters += '--fitMass  '

        command = f"python3 ComputeDstarCFTrick.py {parameters}"
        fempy.info(f"Running: {command}")
        os.system(command)

    def output(self):
        return luigi.LocalTarget(f'{analysisBasePath}/{self.data_version}/cf/RawCF_{self.make_suffix()}.root')

    def requires(self):
        return MakeDirTreeTask(data_version=self.data_version)

class ComputeGenCFTask(CharmFemtoTask):
    charmMassBW = luigi.OptionalParameter(default=None)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.suffix = MakeSuffix(cleaner=self.cleaner, kStarBW=self.kStarBW, charmMassBW=self.charmMassBW)


    def run(self):
        inFile = os.path.join(analysisBasePath, self.data_version, f'cf/RawCF_data{self.suffix}.root')
        cfg = f'/home/daniel/phsw/fempy/cfg_gencf_{pair}_{self.kStarBW}MeV.yml'
        oFile = os.path.join(analysisBasePath, self.data_version, f'cf/GenCF{self.suffix}.root')
        command = f'python3 ComputeGenCF.py {inFile} {cfg} {oFile}'

        os.system(command)

    def output(self):
        return luigi.LocalTarget(os.path.join(analysisBasePath, self.data_version, f'cf/GenCF{self.suffix}.root'))


class FitGenCFTask(CharmFemtoTask):
    charmMassBW = luigi.OptionalParameter(default=None)
    bl = luigi.OptionalStrParameter(default=None)
    fitRange = luigi.OptionalTupleParameter(default=(10, 500))
    bootstrap = luigi.OptionalIntParameter(default=0)
    syst = luigi.OptionalBoolParameter(default=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.oSuffix = MakeSuffix(cleaner=self.cleaner, kStarBW=self.kStarBW, charmMassBW=self.charmMassBW, bl=self.bl, fitRange=self.fitRange, bootstrap=self.bootstrap, syst=self.syst)
        self.inSuffix = MakeSuffix(cleaner=self.cleaner, kStarBW=self.kStarBW, charmMassBW=self.charmMassBW)


    def run(self):
        inFile = os.path.join(analysisBasePath, self.data_version, f'cf/GenCF{self.inSuffix}.root')
        oFile = os.path.join(analysisBasePath, self.data_version, f'cf/GenCFFit{self.oSuffix}.root')
        command = f'python3 FitCF.py --pair {pair} {inFile} {oFile} --fitRange {int(self.fitRange[0])} {int(self.fitRange[1])}'
        if self.bl is not None:
            command += ' --bl ' + self.bl
        if self.bootstrap is not None:
            command += ' --bootstrap ' + str(self.bootstrap)
        if self.syst is not None:
            command += ' --syst '

        print(command)
        os.system(command)

    def output(self):
        oFile = os.path.join(analysisBasePath, self.data_version, f'cf/GenCFFit{self.oSuffix}.root')
        print(oFile)
        return luigi.LocalTarget(oFile)


class DPiTask(CharmFemtoTask):
    def run(self):
        for mc in ['mcgp', 'mcgptruthreco']:
            fromc = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromc.root"
            fromb = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromb.root"
            fromu = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromu.root"
            fromq = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromq.root"
            true = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_true.root"

            if os.path.isfile(fromc) and os.path.isfile(fromb):
                print(f'---> hadd -n 10 {fromq} {fromc} {fromb}')
                os.system(f'hadd -n 10 {fromq} {fromc} {fromb}')
            if os.path.isfile(fromc) and os.path.isfile(fromb) and os.path.isfile(fromu):
                print(f'---> hadd -n 10 {true} {fromc} {fromb} {fromu}')
                os.system(f'hadd -n 10 {true} {fromc} {fromb} {fromu}')
        
    def requires(self):
        # yield MakeDirTreeTask(data_version=self.data_version)

        if False and 'data' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="data", nJobs=1, cleaner=self.cleaner, doSyst='true', doQA='false')

        if False and 'mcgp' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromu')
        if False and 'mchf' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')
        if 'mcgptruthreco' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromu')
        if False and 'mchftruthreco' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mchftruthreco", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mchftruthreco", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')
            yield MakeDistrTask(data_version=self.data_version, sample="mchftruthreco", nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromu')



class DKTask(CharmFemtoTask):
    def run(self):
        # Merge MCGP
        for mc in ['mcgp', 'mcgptruthreco']:
            fromc = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromc.root"
            fromb = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromb.root"
            fromu = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromu.root"
            fromq = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_fromq.root"
            true = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_true.root"

            if os.path.isfile(fromc) and os.path.isfile(fromb):
                print(f'---> hadd -n 10 {fromq} {fromc} {fromb}')
                os.system(f'hadd -n 10 {fromq} {fromc} {fromb}')
            if os.path.isfile(fromc) and os.path.isfile(fromb) and os.path.isfile(fromu):
                print(f'---> hadd -n 10 {true} {fromc} {fromb} {fromu}')
                os.system(f'hadd -n 10 {true} {fromc} {fromb} {fromu}')

    def requires(self):
        yield MakeDirTreeTask(data_version=self.data_version)
        yield MakeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='true', doQA='true')

        if 'mcgp' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromu')
        if 'mchf' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')

        if 'mcgptruthreco' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromu')
        if 'mchftruthreco' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mchftruthreco", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mchftruthreco", kStarBW=self.kStarBW, nJobs=1, cleaner=self.cleaner, doSyst='false', doQA='true', origin='fromb')


class DstarKTask(CharmFemtoTask):
    def run(self):
        print(f'{analysisBasePath}/{self.data_version}/comparisons/rawcf/cfg_compare_stdSBCorr.yml')

        # yield CompareGraphs(f'~/an/DstarPi/comparisons/rawcf/cfg_compare_stdSBCorr.yml')

    def requires(self):
        
        if 'data' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="data", nJobs=12, kStarBW=self.kStarBW, cleaner=self.cleaner, doSyst='true', doQA='false')
            yield ComputeCFTask(
                data_version=self.data_version,
                sample="data",
                cleaner=self.cleaner,
                kStarBW=self.kStarBW,
                aliFitter=False),
        if 'mchf' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", nJobs=12, kStarBW=self.kStarBW, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", nJobs=12, kStarBW=self.kStarBW, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromb')

        if 'mcgp' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromb')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromu')
            # yield ComputeCFTask(
            #     data_version=self.data_version,
            #     sample="mcgp",
            #     cleaner=self.cleaner,
            #     kStarBW=self.kStarBW,
            #     charmMassBW=0.2,
            #     origin='true'),
        if 'data' in datasets[self.data_version] and 'mcgp' in datasets[self.data_version]:
            yield ComputeGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW)
            yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 400))
            yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 500))
            yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 1000))

            yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 500), bootstrap=30, syst=True)

class DstarPiTask(CharmFemtoTask):
    def run(self):
        pass
        # for mc in ['mcgp', 'mcgptruthreco']:
        #     fromc = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_kStarBW{self.kStarBW}MeV_fromc.root"
        #     fromb = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_kStarBW{self.kStarBW}MeV_fromb.root"
        #     fromu = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_kStarBW{self.kStarBW}MeV_fromu.root"
        #     fromq = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_kStarBW{self.kStarBW}MeV_fromq.root"
        #     true = f"{analysisBasePath}/{self.data_version}/distr/Distr_{mc}_{self.cleaner}_kStarBW{self.kStarBW}MeV_true.root"

        #     if os.path.isfile(fromc) and os.path.isfile(fromb):
        #         print(f'---> hadd -n 10 {fromq} {fromc} {fromb}')
        #         os.system(f'hadd -n 10 {fromq} {fromc} {fromb}')
        #     if os.path.isfile(fromc) and os.path.isfile(fromb) and os.path.isfile(fromu):
        #         print(f'---> hadd -n 10 {true} {fromc} {fromb} {fromu}')
        #         os.system(f'hadd -n 10 {true} {fromc} {fromb} {fromu}')
    #     yield CompareGraphs(f'{analysisBasePath}/{self.data_version}/comparisons/uniq/cfg_compare_uniq.yml')

    def requires(self):
        if 'data' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="data", kStarBW=self.kStarBW, nJobs=10, cleaner=self.cleaner, doSyst='true', doQA='false')
            yield ComputeCFTask(
                data_version=self.data_version,
                sample="data",
                cleaner=self.cleaner,
                kStarBW=self.kStarBW,
                sgnFitFunc="gaus",
                bkgFitFunc="powex",
                charmMassBW=0.2,
                lowKStarCharmMassBW=0.4,
                lowKStarMax=50,
                fitMass=True,
                aliFitter=False)
            # yield ComputeCFTask(
            #     data_version=self.data_version,
            #     sample="data",
            #     cleaner=self.cleaner,
            #     kStarBW=self.kStarBW,
            #     sgnFitFunc="gaus",
            #     bkgFitFunc="powex",
            #     charmMassBW=0.2,
            #     lowKStarCharmMassBW=0.4,
            #     lowKStarMax=50,
            #     fitMass=True,
            #     uniq=True,
            #     aliFitter=False)
        if 'mchf' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mchf", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromb')

        if 'mcgp' in datasets[self.data_version]:
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromc')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromb')
            yield MakeDistrTask(data_version=self.data_version, sample="mcgp", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromu')
            print('-------------> ', self.kStarBW)
            yield ComputeCFTask(
                data_version=self.data_version,
                sample="mcgp",
                cleaner=self.cleaner,
                kStarBW=self.kStarBW,
                origin='true')
        if 'mcgptruthreco' in datasets[self.data_version]:
            # yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromc')
            # yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromb')
            # yield MakeDistrTask(data_version=self.data_version, sample="mcgptruthreco", kStarBW=self.kStarBW, nJobs=12, cleaner=self.cleaner, doSyst='false', doQA='false', origin='fromu')
            yield ComputeCFTask(
                data_version=self.data_version,
                sample="mcgptruthreco",
                cleaner=self.cleaner,
                kStarBW=self.kStarBW,
                origin='true')


        # if 'data' in datasets[self.data_version] and 'mcgp' in datasets[self.data_version]:
        #     yield ComputeGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW)
        #     # yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 400))
        #     # yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 500))
        #     # yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 1000))

        #     yield FitGenCFTask(data_version=self.data_version, cleaner=self.cleaner, kStarBW=self.kStarBW, bl='const', fitRange=(10, 500), bootstrap=30, syst=True)


        # yield ComputeCFTask(
        #     data_version=self.data_version,
        #     sample="data",
        #     cleaner=self.cleaner,
        #     kStarBW=self.kStarBW,
        #     sgnFitFunc="gaus",
        #     bkgFitFunc="powex",
        #     charmMassBW=0.2,
        #     lowKStarCharmMassBW=0.4,
        #     lowKStarMax=50,
        #     fitMass=True,
        #     aliFitter=True)
        # yield ComputeCFTask(
        #     data_version=self.data_version,
        #     sample="data",
        #     cleaner=self.cleaner,
        #     kStarBW=self.kStarBW,
        #     sgnFitFunc="gaus",
        #     bkgFitFunc="powex",
        #     charmMassBW=0.2,
        #     lowKStarCharmMassBW=0.4,
        #     lowKStarMax=50,
        #     fitMass=True,
        #     aliFitter=False)
        # yield ComputeCFTask(
        #     data_version=self.data_version,
        #     sample="data",
        #     cleaner=self.cleaner,
        #     kStarBW=self.kStarBW,
        #     sgnFitFunc="hat",
        #     bkgFitFunc="powex",
        #     charmMassBW=0.2,
        #     lowKStarCharmMassBW=0.4,
        #     lowKStarMax=50,
        #     fitMass=True,
        #     aliFitter=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("pair", choices=('DstarPi', 'DstarK', 'DPi', 'DK', 'DPiClosure'))
    args = parser.parse_args()

    pair = args.pair

    analysisBasePath = f'/home/daniel/an/{pair}'
    if pair == 'DstarPi':

        datasets = {
            "16_mixpions": {
                "mcgp": [
                    ("/data/DstarPi/16_mixpions/mcgp/AnalysisResults_4044.root", "4044"),
                    ("/data/DstarPi/16_mixpions/mcgp/AnalysisResults_4046.root", "4046"),
                    ("/data/DstarPi/16_mixpions/mcgp/AnalysisResults_4045.root", "4045"),
                ]
            },
            "18_fixmix": {
                "data": [
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4654.root", "4654"),
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4655.root", "4655"),
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4656.root", "4656"),
                    ("/data/DstarPi/18_fixmix/data/AnalysisResults_4657.root", "4657"),
                ],
            },
            "19_mix150": {
                "mcgptruthreco": [
                    ("/data/DstarPi/19_mix150/mcgptruthreco/AnalysisResults_4217.root", "4217"),
                    ("/data/DstarPi/19_mix150/mcgptruthreco/AnalysisResults_4218.root", "4218"),
                    ("/data/DstarPi/19_mix150/mcgptruthreco/AnalysisResults_4219.root", "4219"),
                ]
            },
            "20_luuksel": {
                "data": [
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/001/AnalysisResults.root", "4726_001"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/002/AnalysisResults.root", "4726_002"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/003/AnalysisResults.root", "4726_003"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/004/AnalysisResults.root", "4726_004"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/005/AnalysisResults.root", "4726_005"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/006/AnalysisResults.root", "4726_006"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/007/AnalysisResults.root", "4726_007"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/008/AnalysisResults.root", "4726_008"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/009/AnalysisResults.root", "4726_009"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/010/AnalysisResults.root", "4726_010"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/011/AnalysisResults.root", "4726_011"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/012/AnalysisResults.root", "4726_012"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/013/AnalysisResults.root", "4726_013"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/014/AnalysisResults.root", "4726_014"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/015/AnalysisResults.root", "4726_015"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/016/AnalysisResults.root", "4726_016"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/017/AnalysisResults.root", "4726_017"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/018/AnalysisResults.root", "4726_018"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/019/AnalysisResults.root", "4726_019"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/020/AnalysisResults.root", "4726_020"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/021/AnalysisResults.root", "4726_021"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/022/AnalysisResults.root", "4726_022"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/023/AnalysisResults.root", "4726_023"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/024/AnalysisResults.root", "4726_024"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/025/AnalysisResults.root", "4726_025"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/026/AnalysisResults.root", "4726_026"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/027/AnalysisResults.root", "4726_027"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/028/AnalysisResults.root", "4726_028"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/029/AnalysisResults.root", "4726_029"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/030/AnalysisResults.root", "4726_030"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/031/AnalysisResults.root", "4726_031"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/032/AnalysisResults.root", "4726_032"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/033/AnalysisResults.root", "4726_033"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/034/AnalysisResults.root", "4726_034"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/035/AnalysisResults.root", "4726_035"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/036/AnalysisResults.root", "4726_036"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/037/AnalysisResults.root", "4726_037"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/038/AnalysisResults.root", "4726_038"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/039/AnalysisResults.root", "4726_039"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/040/AnalysisResults.root", "4726_040"),
                    ("/data/DstarPi/20_luuksel/data/merge_4726/Stage_1/041/AnalysisResults.root", "4726_041"),
                    ("/data/DstarPi/20_luuksel/data/merge_4727/Stage_3/001/AnalysisResults.root", "4727_001"),
                    ("/data/DstarPi/20_luuksel/data/merge_4727/Stage_3/002/AnalysisResults.root", "4727_002"),
                    ("/data/DstarPi/20_luuksel/data/merge_4727/Stage_3/003/AnalysisResults.root", "4727_003"),
                    ("/data/DstarPi/20_luuksel/data/merge_4727/Stage_3/004/AnalysisResults.root", "4727_004"),
                    ("/data/DstarPi/20_luuksel/data/merge_4727/Stage_3/006/AnalysisResults.root", "4727_006"),
                    ("/data/DstarPi/20_luuksel/data/merge_4727/Stage_3/008/AnalysisResults.root", "4727_008"),
                    ("/data/DstarPi/20_luuksel/data/merge_4727/Stage_3/009/AnalysisResults.root", "4727_009"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/001/AnalysisResults.root", "4729_001"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/002/AnalysisResults.root", "4729_002"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/003/AnalysisResults.root", "4729_003"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/004/AnalysisResults.root", "4729_004"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/005/AnalysisResults.root", "4729_005"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/006/AnalysisResults.root", "4729_006"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/007/AnalysisResults.root", "4729_007"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/008/AnalysisResults.root", "4729_008"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/009/AnalysisResults.root", "4729_009"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/010/AnalysisResults.root", "4729_010"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/011/AnalysisResults.root", "4729_011"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/012/AnalysisResults.root", "4729_012"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/013/AnalysisResults.root", "4729_013"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/014/AnalysisResults.root", "4729_014"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/015/AnalysisResults.root", "4729_015"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/016/AnalysisResults.root", "4729_016"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/017/AnalysisResults.root", "4729_017"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/018/AnalysisResults.root", "4729_018"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/019/AnalysisResults.root", "4729_019"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/020/AnalysisResults.root", "4729_020"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/021/AnalysisResults.root", "4729_021"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/022/AnalysisResults.root", "4729_022"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/023/AnalysisResults.root", "4729_023"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/024/AnalysisResults.root", "4729_024"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/025/AnalysisResults.root", "4729_025"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/026/AnalysisResults.root", "4729_026"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/027/AnalysisResults.root", "4729_027"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/028/AnalysisResults.root", "4729_028"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/029/AnalysisResults.root", "4729_029"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/030/AnalysisResults.root", "4729_030"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/031/AnalysisResults.root", "4729_031"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/032/AnalysisResults.root", "4729_032"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/033/AnalysisResults.root", "4729_033"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/034/AnalysisResults.root", "4729_034"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/035/AnalysisResults.root", "4729_035"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/036/AnalysisResults.root", "4729_036"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/037/AnalysisResults.root", "4729_037"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/038/AnalysisResults.root", "4729_038"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/039/AnalysisResults.root", "4729_039"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/040/AnalysisResults.root", "4729_040"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/041/AnalysisResults.root", "4729_041"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/042/AnalysisResults.root", "4729_042"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/043/AnalysisResults.root", "4729_043"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/044/AnalysisResults.root", "4729_044"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/045/AnalysisResults.root", "4729_045"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/046/AnalysisResults.root", "4729_046"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/047/AnalysisResults.root", "4729_047"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/048/AnalysisResults.root", "4729_048"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/049/AnalysisResults.root", "4729_049"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/050/AnalysisResults.root", "4729_050"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/051/AnalysisResults.root", "4729_051"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/052/AnalysisResults.root", "4729_052"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/053/AnalysisResults.root", "4729_053"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/057/AnalysisResults.root", "4729_057"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/058/AnalysisResults.root", "4729_058"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/059/AnalysisResults.root", "4729_059"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/060/AnalysisResults.root", "4729_060"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/061/AnalysisResults.root", "4729_061"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/062/AnalysisResults.root", "4729_062"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/063/AnalysisResults.root", "4729_063"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/064/AnalysisResults.root", "4729_064"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/065/AnalysisResults.root", "4729_065"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/066/AnalysisResults.root", "4729_066"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/067/AnalysisResults.root", "4729_067"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/068/AnalysisResults.root", "4729_068"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/069/AnalysisResults.root", "4729_069"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/070/AnalysisResults.root", "4729_070"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/071/AnalysisResults.root", "4729_071"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/072/AnalysisResults.root", "4729_072"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/073/AnalysisResults.root", "4729_073"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/074/AnalysisResults.root", "4729_074"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/075/AnalysisResults.root", "4729_075"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/076/AnalysisResults.root", "4729_076"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/077/AnalysisResults.root", "4729_077"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/078/AnalysisResults.root", "4729_078"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/079/AnalysisResults.root", "4729_079"),
                    ("/data/DstarPi/20_luuksel/data/merge_4729/Stage_1/080/AnalysisResults.root", "4729_080"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/001/AnalysisResults.root", "4730_001"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/002/AnalysisResults.root", "4730_002"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/003/AnalysisResults.root", "4730_003"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/004/AnalysisResults.root", "4730_004"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/005/AnalysisResults.root", "4730_005"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/006/AnalysisResults.root", "4730_006"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/007/AnalysisResults.root", "4730_007"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/008/AnalysisResults.root", "4730_008"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/009/AnalysisResults.root", "4730_009"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/010/AnalysisResults.root", "4730_010"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/011/AnalysisResults.root", "4730_011"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/012/AnalysisResults.root", "4730_012"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/013/AnalysisResults.root", "4730_013"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/014/AnalysisResults.root", "4730_014"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/015/AnalysisResults.root", "4730_015"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/016/AnalysisResults.root", "4730_016"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/017/AnalysisResults.root", "4730_017"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/018/AnalysisResults.root", "4730_018"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/019/AnalysisResults.root", "4730_019"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/020/AnalysisResults.root", "4730_020"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/021/AnalysisResults.root", "4730_021"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/022/AnalysisResults.root", "4730_022"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/023/AnalysisResults.root", "4730_023"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/024/AnalysisResults.root", "4730_024"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/025/AnalysisResults.root", "4730_025"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/026/AnalysisResults.root", "4730_026"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/027/AnalysisResults.root", "4730_027"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/028/AnalysisResults.root", "4730_028"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/029/AnalysisResults.root", "4730_029"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/030/AnalysisResults.root", "4730_030"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/031/AnalysisResults.root", "4730_031"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/032/AnalysisResults.root", "4730_032"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/033/AnalysisResults.root", "4730_033"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/034/AnalysisResults.root", "4730_034"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/035/AnalysisResults.root", "4730_035"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/036/AnalysisResults.root", "4730_036"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/037/AnalysisResults.root", "4730_037"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/038/AnalysisResults.root", "4730_038"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/039/AnalysisResults.root", "4730_039"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/040/AnalysisResults.root", "4730_040"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/041/AnalysisResults.root", "4730_041"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/042/AnalysisResults.root", "4730_042"),
                    ("/data/DstarPi/20_luuksel/data/merge_4730/Stage_1/043/AnalysisResults.root", "4730_043"),
                ],
                'mcgp': [
                    ("/data/DstarPi/20_luuksel/mcgp/AnalysisResults_4251.root", "4251"),
                    ("/data/DstarPi/20_luuksel/mcgp/AnalysisResults_4252.root", "4252"),
                    ("/data/DstarPi/20_luuksel/mcgp/AnalysisResults_4253.root", "4253"),
                ],
                'mchf': [
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/001/AnalysisResults.root", "4306_001"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/002/AnalysisResults.root", "4306_002"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/003/AnalysisResults.root", "4306_003"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/004/AnalysisResults.root", "4306_004"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/005/AnalysisResults.root", "4306_005"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/006/AnalysisResults.root", "4306_006"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/007/AnalysisResults.root", "4306_007"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/008/AnalysisResults.root", "4306_008"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/009/AnalysisResults.root", "4306_009"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/010/AnalysisResults.root", "4306_010"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/011/AnalysisResults.root", "4306_011"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/012/AnalysisResults.root", "4306_012"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/013/AnalysisResults.root", "4306_013"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/014/AnalysisResults.root", "4306_014"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/015/AnalysisResults.root", "4306_015"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/016/AnalysisResults.root", "4306_016"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/017/AnalysisResults.root", "4306_017"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/018/AnalysisResults.root", "4306_018"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/019/AnalysisResults.root", "4306_019"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/020/AnalysisResults.root", "4306_020"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/021/AnalysisResults.root", "4306_021"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/023/AnalysisResults.root", "4306_023"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/024/AnalysisResults.root", "4306_024"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/025/AnalysisResults.root", "4306_025"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4306/Stage_1/026/AnalysisResults.root", "4306_026"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/001/AnalysisResults.root", "4307_001"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/002/AnalysisResults.root", "4307_002"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/003/AnalysisResults.root", "4307_003"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/004/AnalysisResults.root", "4307_004"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/005/AnalysisResults.root", "4307_005"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/006/AnalysisResults.root", "4307_006"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/007/AnalysisResults.root", "4307_007"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/008/AnalysisResults.root", "4307_008"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/009/AnalysisResults.root", "4307_009"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/010/AnalysisResults.root", "4307_010"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/011/AnalysisResults.root", "4307_011"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/012/AnalysisResults.root", "4307_012"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/013/AnalysisResults.root", "4307_013"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/014/AnalysisResults.root", "4307_014"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/015/AnalysisResults.root", "4307_015"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/016/AnalysisResults.root", "4307_016"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/017/AnalysisResults.root", "4307_017"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/018/AnalysisResults.root", "4307_018"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/019/AnalysisResults.root", "4307_019"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/020/AnalysisResults.root", "4307_020"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/021/AnalysisResults.root", "4307_021"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/022/AnalysisResults.root", "4307_022"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/023/AnalysisResults.root", "4307_023"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/024/AnalysisResults.root", "4307_024"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/025/AnalysisResults.root", "4307_025"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/026/AnalysisResults.root", "4307_026"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/027/AnalysisResults.root", "4307_027"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/028/AnalysisResults.root", "4307_028"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/029/AnalysisResults.root", "4307_029"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/030/AnalysisResults.root", "4307_030"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/031/AnalysisResults.root", "4307_031"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/032/AnalysisResults.root", "4307_032"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/033/AnalysisResults.root", "4307_033"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/034/AnalysisResults.root", "4307_034"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/035/AnalysisResults.root", "4307_035"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/036/AnalysisResults.root", "4307_036"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/037/AnalysisResults.root", "4307_037"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/038/AnalysisResults.root", "4307_038"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/039/AnalysisResults.root", "4307_039"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/040/AnalysisResults.root", "4307_040"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/041/AnalysisResults.root", "4307_041"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4307/Stage_1/042/AnalysisResults.root", "4307_042"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/001/AnalysisResults.root", "4308_001"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/002/AnalysisResults.root", "4308_002"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/003/AnalysisResults.root", "4308_003"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/004/AnalysisResults.root", "4308_004"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/005/AnalysisResults.root", "4308_005"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/006/AnalysisResults.root", "4308_006"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/007/AnalysisResults.root", "4308_007"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/008/AnalysisResults.root", "4308_008"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/009/AnalysisResults.root", "4308_009"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/010/AnalysisResults.root", "4308_010"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/011/AnalysisResults.root", "4308_011"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/012/AnalysisResults.root", "4308_012"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/013/AnalysisResults.root", "4308_013"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/014/AnalysisResults.root", "4308_014"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/015/AnalysisResults.root", "4308_015"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/016/AnalysisResults.root", "4308_016"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/017/AnalysisResults.root", "4308_017"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/018/AnalysisResults.root", "4308_018"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/019/AnalysisResults.root", "4308_019"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/020/AnalysisResults.root", "4308_020"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/021/AnalysisResults.root", "4308_021"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/022/AnalysisResults.root", "4308_022"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/023/AnalysisResults.root", "4308_023"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/024/AnalysisResults.root", "4308_024"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/025/AnalysisResults.root", "4308_025"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/026/AnalysisResults.root", "4308_026"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/027/AnalysisResults.root", "4308_027"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/028/AnalysisResults.root", "4308_028"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/029/AnalysisResults.root", "4308_029"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/030/AnalysisResults.root", "4308_030"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/031/AnalysisResults.root", "4308_031"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/032/AnalysisResults.root", "4308_032"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/033/AnalysisResults.root", "4308_033"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/034/AnalysisResults.root", "4308_034"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/035/AnalysisResults.root", "4308_035"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/036/AnalysisResults.root", "4308_036"),
                    ("/data/DstarPi/20_luuksel/mchf/merge_4308/Stage_1/037/AnalysisResults.root", "4308_037"),
                ]
            },
        }

        luigi.build([
            # DstarPiTask(data_version='16_mixpions', kStarBW=10, cleaner='nopc'),
            # DstarPiTask(data_version='18_fixmix', kStarBW=10, cleaner='nopc'),
            # DstarPiTask(data_version='19_mix150', kStarBW=10, cleaner='nopc'),
            # DstarPiTask(data_version='16_mixpions', kStarBW=15, cleaner='nopc'),
            # DstarPiTask(data_version='18_fixmix', kStarBW=15, cleaner='nopc'),
            # DstarPiTask(data_version='19_mix150', kStarBW=15, cleaner='nopc'),
            # DstarPiTask(data_version='16_mixpions', kStarBW=20, cleaner='nopc'),
            # DstarPiTask(data_version='18_fixmix', kStarBW=20, cleaner='nopc'),
            # DstarPiTask(data_version='19_mix150', kStarBW=20, cleaner='nopc'),
            # DstarPiTask(data_version='16_mixpions', kStarBW=25, cleaner='nopc'),
            # DstarPiTask(data_version='18_fixmix', kStarBW=25, cleaner='nopc'),
            # DstarPiTask(data_version='19_mix150', kStarBW=25, cleaner='nopc'),
            # DstarPiTask(data_version='16_mixpions', kStarBW=50, cleaner='nopc'),
            # DstarPiTask(data_version='18_fixmix', kStarBW=50, cleaner='nopc'),
            # DstarPiTask(data_version='19_mix150', kStarBW=50, cleaner='nopc'),
            # DstarPiTask(data_version='20_luuksel', kStarBW=15, cleaner='nopc'),
            # DstarPiTask(data_version='20_luuksel', kStarBW=20, cleaner='nopc'),
            # DstarPiTask(data_version='20_luuksel', kStarBW=25, cleaner='nopc'),
            DstarPiTask(data_version='20_luuksel', kStarBW=50, cleaner='nopc'),
            # DstarPiTask(data_version='20_luuksel', kStarBW=15, cleaner='nopc'),

            # closure test on the fancy fit trick procedure
            # DPiClosureTask(data_version='20_luuksel', kStarBW=50, cleaner='nopc'),
            
            

        ], workers=12)

    elif pair == 'DstarK':
        datasets = {
            "1_trees": {
                "data": [
                    ("/data/DstarK/1_trees/data/AnalysisResults_4708.root", "4708"),
                    ("/data/DstarK/1_trees/data/AnalysisResults_4709.root", "4709"),
                    ("/data/DstarK/1_trees/data/AnalysisResults_4710.root", "4710"),
                    ("/data/DstarK/1_trees/data/AnalysisResults_4711.root", "4711"),
                ],
            },
            "2_luuksel": {
                'mcgp': [
                    ("/data/DstarK/2_luuksel/mcgp/AnalysisResults_4257.root", "4257"),
                    ("/data/DstarK/2_luuksel/mcgp/AnalysisResults_4258.root", "4258"),
                    ("/data/DstarK/2_luuksel/mcgp/AnalysisResults_4259.root", "4259"),
                ],
                'data': [
                    ("/data/DstarK/2_luuksel/data/AnalysisResults_4745.root", "4745"),
                    ("/data/DstarK/2_luuksel/data/AnalysisResults_4746.root", "4746"),
                    ("/data/DstarK/2_luuksel/data/AnalysisResults_4747.root", "4747"),
                    ("/data/DstarK/2_luuksel/data/AnalysisResults_4748.root", "4748"),
                ],
                'mchf': [
                    ("/data/DstarK/2_luuksel/mchf/AnalysisResults_4305.root", "4305"),
                    ("/data/DstarK/2_luuksel/mchf/AnalysisResults_4303.root", "4303"),
                    ("/data/DstarK/2_luuksel/mchf/AnalysisResults_4304.root", "4304"),
                ]
            }
        }

        luigi.build([
            # DstarKTask(data_version='1_trees', kStarBW=25, cleaner='nopc'),
            DstarKTask(data_version='1_trees', kStarBW=50, cleaner='nopc'),
            DstarKTask(data_version='2_luuksel', kStarBW=50, cleaner='nopc'),
        ], workers=12)

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
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/001/AnalysisResults.root", "sbl_4702_001"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/002/AnalysisResults.root", "sbl_4702_002"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/003/AnalysisResults.root", "sbl_4702_003"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/004/AnalysisResults.root", "sbl_4702_004"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/005/AnalysisResults.root", "sbl_4702_005"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/006/AnalysisResults.root", "sbl_4702_006"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/007/AnalysisResults.root", "sbl_4702_007"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/008/AnalysisResults.root", "sbl_4702_008"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/009/AnalysisResults.root", "sbl_4702_009"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/010/AnalysisResults.root", "sbl_4702_010"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/011/AnalysisResults.root", "sbl_4702_011"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/012/AnalysisResults.root", "sbl_4702_012"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/013/AnalysisResults.root", "sbl_4702_013"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/014/AnalysisResults.root", "sbl_4702_014"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/015/AnalysisResults.root", "sbl_4702_015"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/016/AnalysisResults.root", "sbl_4702_016"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/017/AnalysisResults.root", "sbl_4702_017"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/018/AnalysisResults.root", "sbl_4702_018"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/019/AnalysisResults.root", "sbl_4702_019"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/020/AnalysisResults.root", "sbl_4702_020"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/021/AnalysisResults.root", "sbl_4702_021"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/022/AnalysisResults.root", "sbl_4702_022"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/023/AnalysisResults.root", "sbl_4702_023"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/024/AnalysisResults.root", "sbl_4702_024"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/025/AnalysisResults.root", "sbl_4702_025"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/026/AnalysisResults.root", "sbl_4702_026"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/027/AnalysisResults.root", "sbl_4702_027"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/028/AnalysisResults.root", "sbl_4702_028"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/029/AnalysisResults.root", "sbl_4702_029"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/030/AnalysisResults.root", "sbl_4702_030"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/031/AnalysisResults.root", "sbl_4702_031"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/032/AnalysisResults.root", "sbl_4702_032"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/033/AnalysisResults.root", "sbl_4702_033"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/034/AnalysisResults.root", "sbl_4702_034"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/035/AnalysisResults.root", "sbl_4702_035"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4702/Stage_1/036/AnalysisResults.root", "sbl_4702_036"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/001/AnalysisResults.root", "sbl_4703_001"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/002/AnalysisResults.root", "sbl_4703_002"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/003/AnalysisResults.root", "sbl_4703_003"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/004/AnalysisResults.root", "sbl_4703_004"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/005/AnalysisResults.root", "sbl_4703_005"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/006/AnalysisResults.root", "sbl_4703_006"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/007/AnalysisResults.root", "sbl_4703_007"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/008/AnalysisResults.root", "sbl_4703_008"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/009/AnalysisResults.root", "sbl_4703_009"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/010/AnalysisResults.root", "sbl_4703_010"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/011/AnalysisResults.root", "sbl_4703_011"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/012/AnalysisResults.root", "sbl_4703_012"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/013/AnalysisResults.root", "sbl_4703_013"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/014/AnalysisResults.root", "sbl_4703_014"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/015/AnalysisResults.root", "sbl_4703_015"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/016/AnalysisResults.root", "sbl_4703_016"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/017/AnalysisResults.root", "sbl_4703_017"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/018/AnalysisResults.root", "sbl_4703_018"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/019/AnalysisResults.root", "sbl_4703_019"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/020/AnalysisResults.root", "sbl_4703_020"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/021/AnalysisResults.root", "sbl_4703_021"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/022/AnalysisResults.root", "sbl_4703_022"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/023/AnalysisResults.root", "sbl_4703_023"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/024/AnalysisResults.root", "sbl_4703_024"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/025/AnalysisResults.root", "sbl_4703_025"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/026/AnalysisResults.root", "sbl_4703_026"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/027/AnalysisResults.root", "sbl_4703_027"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/028/AnalysisResults.root", "sbl_4703_028"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4703/Stage_1/029/AnalysisResults.root", "sbl_4703_029"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/001/AnalysisResults.root", "sbl_4704_001"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/002/AnalysisResults.root", "sbl_4704_002"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/003/AnalysisResults.root", "sbl_4704_003"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/004/AnalysisResults.root", "sbl_4704_004"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/005/AnalysisResults.root", "sbl_4704_005"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/006/AnalysisResults.root", "sbl_4704_006"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/007/AnalysisResults.root", "sbl_4704_007"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/008/AnalysisResults.root", "sbl_4704_008"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/009/AnalysisResults.root", "sbl_4704_009"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/010/AnalysisResults.root", "sbl_4704_010"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/011/AnalysisResults.root", "sbl_4704_011"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/012/AnalysisResults.root", "sbl_4704_012"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/013/AnalysisResults.root", "sbl_4704_013"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/014/AnalysisResults.root", "sbl_4704_014"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/015/AnalysisResults.root", "sbl_4704_015"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/016/AnalysisResults.root", "sbl_4704_016"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/017/AnalysisResults.root", "sbl_4704_017"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/018/AnalysisResults.root", "sbl_4704_018"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/019/AnalysisResults.root", "sbl_4704_019"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/020/AnalysisResults.root", "sbl_4704_020"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/021/AnalysisResults.root", "sbl_4704_021"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/022/AnalysisResults.root", "sbl_4704_022"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/023/AnalysisResults.root", "sbl_4704_023"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/024/AnalysisResults.root", "sbl_4704_024"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/025/AnalysisResults.root", "sbl_4704_025"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/026/AnalysisResults.root", "sbl_4704_026"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/027/AnalysisResults.root", "sbl_4704_027"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/028/AnalysisResults.root", "sbl_4704_028"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/029/AnalysisResults.root", "sbl_4704_029"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/030/AnalysisResults.root", "sbl_4704_030"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4704/Stage_2/031/AnalysisResults.root", "sbl_4704_031"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/001/AnalysisResults.root", "sbl_4705_001"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/002/AnalysisResults.root", "sbl_4705_002"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/003/AnalysisResults.root", "sbl_4705_003"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/004/AnalysisResults.root", "sbl_4705_004"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/005/AnalysisResults.root", "sbl_4705_005"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/006/AnalysisResults.root", "sbl_4705_006"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/007/AnalysisResults.root", "sbl_4705_007"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/008/AnalysisResults.root", "sbl_4705_008"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/009/AnalysisResults.root", "sbl_4705_009"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/010/AnalysisResults.root", "sbl_4705_010"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/011/AnalysisResults.root", "sbl_4705_011"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/012/AnalysisResults.root", "sbl_4705_012"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/013/AnalysisResults.root", "sbl_4705_013"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/014/AnalysisResults.root", "sbl_4705_014"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/015/AnalysisResults.root", "sbl_4705_015"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/016/AnalysisResults.root", "sbl_4705_016"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/017/AnalysisResults.root", "sbl_4705_017"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/018/AnalysisResults.root", "sbl_4705_018"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/019/AnalysisResults.root", "sbl_4705_019"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/020/AnalysisResults.root", "sbl_4705_020"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/021/AnalysisResults.root", "sbl_4705_021"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/022/AnalysisResults.root", "sbl_4705_022"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/023/AnalysisResults.root", "sbl_4705_023"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/024/AnalysisResults.root", "sbl_4705_024"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/025/AnalysisResults.root", "sbl_4705_025"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/026/AnalysisResults.root", "sbl_4705_026"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/027/AnalysisResults.root", "sbl_4705_027"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/028/AnalysisResults.root", "sbl_4705_028"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/029/AnalysisResults.root", "sbl_4705_029"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/030/AnalysisResults.root", "sbl_4705_030"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/031/AnalysisResults.root", "sbl_4705_031"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/032/AnalysisResults.root", "sbl_4705_032"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/033/AnalysisResults.root", "sbl_4705_033"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/034/AnalysisResults.root", "sbl_4705_034"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/035/AnalysisResults.root", "sbl_4705_035"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/036/AnalysisResults.root", "sbl_4705_036"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/037/AnalysisResults.root", "sbl_4705_037"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/038/AnalysisResults.root", "sbl_4705_038"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/039/AnalysisResults.root", "sbl_4705_039"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/040/AnalysisResults.root", "sbl_4705_040"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/041/AnalysisResults.root", "sbl_4705_041"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/042/AnalysisResults.root", "sbl_4705_042"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/043/AnalysisResults.root", "sbl_4705_043"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/044/AnalysisResults.root", "sbl_4705_044"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/045/AnalysisResults.root", "sbl_4705_045"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/046/AnalysisResults.root", "sbl_4705_046"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/047/AnalysisResults.root", "sbl_4705_047"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/048/AnalysisResults.root", "sbl_4705_048"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/049/AnalysisResults.root", "sbl_4705_049"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/050/AnalysisResults.root", "sbl_4705_050"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/051/AnalysisResults.root", "sbl_4705_051"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/052/AnalysisResults.root", "sbl_4705_052"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/053/AnalysisResults.root", "sbl_4705_053"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/054/AnalysisResults.root", "sbl_4705_054"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/055/AnalysisResults.root", "sbl_4705_055"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/056/AnalysisResults.root", "sbl_4705_056"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/057/AnalysisResults.root", "sbl_4705_057"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/058/AnalysisResults.root", "sbl_4705_058"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/059/AnalysisResults.root", "sbl_4705_059"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/060/AnalysisResults.root", "sbl_4705_060"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/061/AnalysisResults.root", "sbl_4705_061"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/062/AnalysisResults.root", "sbl_4705_062"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/063/AnalysisResults.root", "sbl_4705_063"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/064/AnalysisResults.root", "sbl_4705_064"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/065/AnalysisResults.root", "sbl_4705_065"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/066/AnalysisResults.root", "sbl_4705_066"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/067/AnalysisResults.root", "sbl_4705_067"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/068/AnalysisResults.root", "sbl_4705_068"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/069/AnalysisResults.root", "sbl_4705_069"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/070/AnalysisResults.root", "sbl_4705_070"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/071/AnalysisResults.root", "sbl_4705_071"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/072/AnalysisResults.root", "sbl_4705_072"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/073/AnalysisResults.root", "sbl_4705_073"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/074/AnalysisResults.root", "sbl_4705_074"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/075/AnalysisResults.root", "sbl_4705_075"),
                    ("/data/DPi/11_mix150/data/sbl/merge_4705/Stage_1/076/AnalysisResults.root", "sbl_4705_076"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/001/AnalysisResults.root", "sbr_4702_001"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/002/AnalysisResults.root", "sbr_4702_002"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/003/AnalysisResults.root", "sbr_4702_003"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/004/AnalysisResults.root", "sbr_4702_004"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/005/AnalysisResults.root", "sbr_4702_005"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/006/AnalysisResults.root", "sbr_4702_006"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/007/AnalysisResults.root", "sbr_4702_007"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/008/AnalysisResults.root", "sbr_4702_008"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/009/AnalysisResults.root", "sbr_4702_009"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/010/AnalysisResults.root", "sbr_4702_010"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/011/AnalysisResults.root", "sbr_4702_011"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/012/AnalysisResults.root", "sbr_4702_012"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/013/AnalysisResults.root", "sbr_4702_013"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/014/AnalysisResults.root", "sbr_4702_014"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/015/AnalysisResults.root", "sbr_4702_015"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/016/AnalysisResults.root", "sbr_4702_016"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/017/AnalysisResults.root", "sbr_4702_017"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/018/AnalysisResults.root", "sbr_4702_018"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/019/AnalysisResults.root", "sbr_4702_019"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/020/AnalysisResults.root", "sbr_4702_020"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/021/AnalysisResults.root", "sbr_4702_021"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/022/AnalysisResults.root", "sbr_4702_022"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/023/AnalysisResults.root", "sbr_4702_023"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/024/AnalysisResults.root", "sbr_4702_024"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/025/AnalysisResults.root", "sbr_4702_025"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/026/AnalysisResults.root", "sbr_4702_026"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/027/AnalysisResults.root", "sbr_4702_027"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/028/AnalysisResults.root", "sbr_4702_028"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/029/AnalysisResults.root", "sbr_4702_029"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/030/AnalysisResults.root", "sbr_4702_030"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/031/AnalysisResults.root", "sbr_4702_031"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/032/AnalysisResults.root", "sbr_4702_032"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/033/AnalysisResults.root", "sbr_4702_033"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/034/AnalysisResults.root", "sbr_4702_034"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/035/AnalysisResults.root", "sbr_4702_035"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4702/Stage_1/036/AnalysisResults.root", "sbr_4702_036"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/001/AnalysisResults.root", "sbr_4703_001"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/002/AnalysisResults.root", "sbr_4703_002"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/003/AnalysisResults.root", "sbr_4703_003"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/004/AnalysisResults.root", "sbr_4703_004"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/005/AnalysisResults.root", "sbr_4703_005"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/006/AnalysisResults.root", "sbr_4703_006"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/007/AnalysisResults.root", "sbr_4703_007"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/008/AnalysisResults.root", "sbr_4703_008"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/009/AnalysisResults.root", "sbr_4703_009"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/010/AnalysisResults.root", "sbr_4703_010"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/011/AnalysisResults.root", "sbr_4703_011"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/012/AnalysisResults.root", "sbr_4703_012"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/013/AnalysisResults.root", "sbr_4703_013"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/014/AnalysisResults.root", "sbr_4703_014"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/015/AnalysisResults.root", "sbr_4703_015"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/016/AnalysisResults.root", "sbr_4703_016"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/017/AnalysisResults.root", "sbr_4703_017"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/018/AnalysisResults.root", "sbr_4703_018"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/019/AnalysisResults.root", "sbr_4703_019"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/020/AnalysisResults.root", "sbr_4703_020"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/021/AnalysisResults.root", "sbr_4703_021"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/022/AnalysisResults.root", "sbr_4703_022"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/023/AnalysisResults.root", "sbr_4703_023"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/024/AnalysisResults.root", "sbr_4703_024"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/025/AnalysisResults.root", "sbr_4703_025"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/026/AnalysisResults.root", "sbr_4703_026"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/027/AnalysisResults.root", "sbr_4703_027"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/028/AnalysisResults.root", "sbr_4703_028"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4703/Stage_1/029/AnalysisResults.root", "sbr_4703_029"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/001/AnalysisResults.root", "sbr_4704_001"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/002/AnalysisResults.root", "sbr_4704_002"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/003/AnalysisResults.root", "sbr_4704_003"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/004/AnalysisResults.root", "sbr_4704_004"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/005/AnalysisResults.root", "sbr_4704_005"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/006/AnalysisResults.root", "sbr_4704_006"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/007/AnalysisResults.root", "sbr_4704_007"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/008/AnalysisResults.root", "sbr_4704_008"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/009/AnalysisResults.root", "sbr_4704_009"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/010/AnalysisResults.root", "sbr_4704_010"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/011/AnalysisResults.root", "sbr_4704_011"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/012/AnalysisResults.root", "sbr_4704_012"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/013/AnalysisResults.root", "sbr_4704_013"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/014/AnalysisResults.root", "sbr_4704_014"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/015/AnalysisResults.root", "sbr_4704_015"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/016/AnalysisResults.root", "sbr_4704_016"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/017/AnalysisResults.root", "sbr_4704_017"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/018/AnalysisResults.root", "sbr_4704_018"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/019/AnalysisResults.root", "sbr_4704_019"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/020/AnalysisResults.root", "sbr_4704_020"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/021/AnalysisResults.root", "sbr_4704_021"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/022/AnalysisResults.root", "sbr_4704_022"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/023/AnalysisResults.root", "sbr_4704_023"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/024/AnalysisResults.root", "sbr_4704_024"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/025/AnalysisResults.root", "sbr_4704_025"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/026/AnalysisResults.root", "sbr_4704_026"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/027/AnalysisResults.root", "sbr_4704_027"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/028/AnalysisResults.root", "sbr_4704_028"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/029/AnalysisResults.root", "sbr_4704_029"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/030/AnalysisResults.root", "sbr_4704_030"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4704/Stage_2/031/AnalysisResults.root", "sbr_4704_031"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/001/AnalysisResults.root", "sbr_4705_001"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/002/AnalysisResults.root", "sbr_4705_002"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/003/AnalysisResults.root", "sbr_4705_003"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/004/AnalysisResults.root", "sbr_4705_004"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/005/AnalysisResults.root", "sbr_4705_005"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/006/AnalysisResults.root", "sbr_4705_006"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/007/AnalysisResults.root", "sbr_4705_007"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/008/AnalysisResults.root", "sbr_4705_008"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/009/AnalysisResults.root", "sbr_4705_009"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/010/AnalysisResults.root", "sbr_4705_010"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/011/AnalysisResults.root", "sbr_4705_011"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/012/AnalysisResults.root", "sbr_4705_012"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/013/AnalysisResults.root", "sbr_4705_013"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/014/AnalysisResults.root", "sbr_4705_014"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/015/AnalysisResults.root", "sbr_4705_015"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/016/AnalysisResults.root", "sbr_4705_016"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/017/AnalysisResults.root", "sbr_4705_017"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/018/AnalysisResults.root", "sbr_4705_018"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/019/AnalysisResults.root", "sbr_4705_019"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/020/AnalysisResults.root", "sbr_4705_020"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/021/AnalysisResults.root", "sbr_4705_021"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/022/AnalysisResults.root", "sbr_4705_022"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/023/AnalysisResults.root", "sbr_4705_023"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/024/AnalysisResults.root", "sbr_4705_024"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/025/AnalysisResults.root", "sbr_4705_025"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/026/AnalysisResults.root", "sbr_4705_026"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/027/AnalysisResults.root", "sbr_4705_027"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/028/AnalysisResults.root", "sbr_4705_028"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/029/AnalysisResults.root", "sbr_4705_029"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/030/AnalysisResults.root", "sbr_4705_030"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/031/AnalysisResults.root", "sbr_4705_031"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/032/AnalysisResults.root", "sbr_4705_032"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/033/AnalysisResults.root", "sbr_4705_033"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/034/AnalysisResults.root", "sbr_4705_034"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/035/AnalysisResults.root", "sbr_4705_035"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/036/AnalysisResults.root", "sbr_4705_036"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/037/AnalysisResults.root", "sbr_4705_037"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/038/AnalysisResults.root", "sbr_4705_038"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/039/AnalysisResults.root", "sbr_4705_039"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/040/AnalysisResults.root", "sbr_4705_040"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/041/AnalysisResults.root", "sbr_4705_041"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/042/AnalysisResults.root", "sbr_4705_042"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/043/AnalysisResults.root", "sbr_4705_043"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/044/AnalysisResults.root", "sbr_4705_044"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/045/AnalysisResults.root", "sbr_4705_045"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/046/AnalysisResults.root", "sbr_4705_046"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/047/AnalysisResults.root", "sbr_4705_047"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/048/AnalysisResults.root", "sbr_4705_048"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/049/AnalysisResults.root", "sbr_4705_049"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/050/AnalysisResults.root", "sbr_4705_050"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/051/AnalysisResults.root", "sbr_4705_051"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/052/AnalysisResults.root", "sbr_4705_052"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/053/AnalysisResults.root", "sbr_4705_053"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/054/AnalysisResults.root", "sbr_4705_054"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/055/AnalysisResults.root", "sbr_4705_055"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/056/AnalysisResults.root", "sbr_4705_056"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/057/AnalysisResults.root", "sbr_4705_057"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/058/AnalysisResults.root", "sbr_4705_058"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/059/AnalysisResults.root", "sbr_4705_059"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/060/AnalysisResults.root", "sbr_4705_060"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/061/AnalysisResults.root", "sbr_4705_061"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/062/AnalysisResults.root", "sbr_4705_062"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/063/AnalysisResults.root", "sbr_4705_063"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/064/AnalysisResults.root", "sbr_4705_064"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/065/AnalysisResults.root", "sbr_4705_065"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/066/AnalysisResults.root", "sbr_4705_066"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/067/AnalysisResults.root", "sbr_4705_067"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/068/AnalysisResults.root", "sbr_4705_068"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/069/AnalysisResults.root", "sbr_4705_069"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/070/AnalysisResults.root", "sbr_4705_070"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/071/AnalysisResults.root", "sbr_4705_071"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/072/AnalysisResults.root", "sbr_4705_072"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/073/AnalysisResults.root", "sbr_4705_073"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/074/AnalysisResults.root", "sbr_4705_074"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/075/AnalysisResults.root", "sbr_4705_075"),
                    ("/data/DPi/11_mix150/data/sbr/merge_4705/Stage_1/076/AnalysisResults.root", "sbr_4705_076"),
                ],
                "mcgptruthreco": [
                    ("/data/DPi/11_mix150/mcgptruthreco/AnalysisResults_4199.root", "4199"),
                    ("/data/DPi/11_mix150/mcgptruthreco/AnalysisResults_4200.root", "4200"),
                    ("/data/DPi/11_mix150/mcgptruthreco/AnalysisResults_4201.root", "4201"),
                ],
                "mchftruthreco": [
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/001/AnalysisResults.root", "4202_001"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/002/AnalysisResults.root", "4202_002"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/003/AnalysisResults.root", "4202_003"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/005/AnalysisResults.root", "4202_005"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/006/AnalysisResults.root", "4202_006"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/007/AnalysisResults.root", "4202_007"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/008/AnalysisResults.root", "4202_008"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/009/AnalysisResults.root", "4202_009"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/010/AnalysisResults.root", "4202_010"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/011/AnalysisResults.root", "4202_011"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/012/AnalysisResults.root", "4202_012"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/013/AnalysisResults.root", "4202_013"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/014/AnalysisResults.root", "4202_014"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/015/AnalysisResults.root", "4202_015"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/016/AnalysisResults.root", "4202_016"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/017/AnalysisResults.root", "4202_017"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/018/AnalysisResults.root", "4202_018"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/019/AnalysisResults.root", "4202_019"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/020/AnalysisResults.root", "4202_020"),
                    ("/data/==DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/021/AnalysisResults.root", "4202_021"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/022/AnalysisResults.root", "4202_022"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/023/AnalysisResults.root", "4202_023"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/024/AnalysisResults.root", "4202_024"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/025/AnalysisResults.root", "4202_025"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/026/AnalysisResults.root", "4202_026"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/027/AnalysisResults.root", "4202_027"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/028/AnalysisResults.root", "4202_028"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/029/AnalysisResults.root", "4202_029"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/030/AnalysisResults.root", "4202_030"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/031/AnalysisResults.root", "4202_031"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/032/AnalysisResults.root", "4202_032"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/033/AnalysisResults.root", "4202_033"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/034/AnalysisResults.root", "4202_034"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/035/AnalysisResults.root", "4202_035"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/036/AnalysisResults.root", "4202_036"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/037/AnalysisResults.root", "4202_037"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/038/AnalysisResults.root", "4202_038"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/039/AnalysisResults.root", "4202_039"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4202/Stage_1/040/AnalysisResults.root", "4202_040"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/001/AnalysisResults.root", "4203_001"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/002/AnalysisResults.root", "4203_002"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/003/AnalysisResults.root", "4203_003"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/004/AnalysisResults.root", "4203_004"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/005/AnalysisResults.root", "4203_005"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/006/AnalysisResults.root", "4203_006"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/007/AnalysisResults.root", "4203_007"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/008/AnalysisResults.root", "4203_008"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/009/AnalysisResults.root", "4203_009"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/010/AnalysisResults.root", "4203_010"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/011/AnalysisResults.root", "4203_011"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/012/AnalysisResults.root", "4203_012"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/013/AnalysisResults.root", "4203_013"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/014/AnalysisResults.root", "4203_014"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/015/AnalysisResults.root", "4203_015"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/017/AnalysisResults.root", "4203_017"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/019/AnalysisResults.root", "4203_019"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/020/AnalysisResults.root", "4203_020"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/021/AnalysisResults.root", "4203_021"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/022/AnalysisResults.root", "4203_022"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/023/AnalysisResults.root", "4203_023"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/024/AnalysisResults.root", "4203_024"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/025/AnalysisResults.root", "4203_025"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/026/AnalysisResults.root", "4203_026"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/027/AnalysisResults.root", "4203_027"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/028/AnalysisResults.root", "4203_028"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/029/AnalysisResults.root", "4203_029"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/030/AnalysisResults.root", "4203_030"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/031/AnalysisResults.root", "4203_031"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/032/AnalysisResults.root", "4203_032"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/033/AnalysisResults.root", "4203_033"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/034/AnalysisResults.root", "4203_034"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/035/AnalysisResults.root", "4203_035"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/036/AnalysisResults.root", "4203_036"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/037/AnalysisResults.root", "4203_037"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/038/AnalysisResults.root", "4203_038"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/039/AnalysisResults.root", "4203_039"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/040/AnalysisResults.root", "4203_040"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/041/AnalysisResults.root", "4203_041"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/042/AnalysisResults.root", "4203_042"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4203/Stage_1/043/AnalysisResults.root", "4203_043"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/001/AnalysisResults.root", "4204_001"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/002/AnalysisResults.root", "4204_002"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/003/AnalysisResults.root", "4204_003"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/004/AnalysisResults.root", "4204_004"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/005/AnalysisResults.root", "4204_005"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/006/AnalysisResults.root", "4204_006"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/007/AnalysisResults.root", "4204_007"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/008/AnalysisResults.root", "4204_008"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/009/AnalysisResults.root", "4204_009"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/010/AnalysisResults.root", "4204_010"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/011/AnalysisResults.root", "4204_011"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/012/AnalysisResults.root", "4204_012"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/014/AnalysisResults.root", "4204_014"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/016/AnalysisResults.root", "4204_016"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/017/AnalysisResults.root", "4204_017"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/018/AnalysisResults.root", "4204_018"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/019/AnalysisResults.root", "4204_019"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/020/AnalysisResults.root", "4204_020"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/021/AnalysisResults.root", "4204_021"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/022/AnalysisResults.root", "4204_022"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/023/AnalysisResults.root", "4204_023"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/024/AnalysisResults.root", "4204_024"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/025/AnalysisResults.root", "4204_025"),
                    ("/data/DPi/11_mix150/mchftruthreco/merge_4204/Stage_1/026/AnalysisResults.root", "4204_026"),
                ],
            },
            "12_mix150_nosepsb": {
                "data": [
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/059/AnalysisResults.root", "4720_059"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/064/AnalysisResults.root", "4720_064"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/043/AnalysisResults.root", "4720_043"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/007/AnalysisResults.root", "4720_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/017/AnalysisResults.root", "4720_017"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/062/AnalysisResults.root", "4720_062"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/069/AnalysisResults.root", "4720_069"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/004/AnalysisResults.root", "4720_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/073/AnalysisResults.root", "4720_073"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/002/AnalysisResults.root", "4720_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/041/AnalysisResults.root", "4720_041"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/067/AnalysisResults.root", "4720_067"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/045/AnalysisResults.root", "4720_045"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/047/AnalysisResults.root", "4720_047"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/022/AnalysisResults.root", "4720_022"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/050/AnalysisResults.root", "4720_050"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/016/AnalysisResults.root", "4720_016"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/066/AnalysisResults.root", "4720_066"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/026/AnalysisResults.root", "4720_026"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/079/AnalysisResults.root", "4720_079"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/006/AnalysisResults.root", "4720_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/072/AnalysisResults.root", "4720_072"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/010/AnalysisResults.root", "4720_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/055/AnalysisResults.root", "4720_055"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/030/AnalysisResults.root", "4720_030"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/028/AnalysisResults.root", "4720_028"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/001/AnalysisResults.root", "4720_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/076/AnalysisResults.root", "4720_076"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/027/AnalysisResults.root", "4720_027"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/056/AnalysisResults.root", "4720_056"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/032/AnalysisResults.root", "4720_032"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/048/AnalysisResults.root", "4720_048"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/070/AnalysisResults.root", "4720_070"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/078/AnalysisResults.root", "4720_078"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/019/AnalysisResults.root", "4720_019"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/034/AnalysisResults.root", "4720_034"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/003/AnalysisResults.root", "4720_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/049/AnalysisResults.root", "4720_049"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/068/AnalysisResults.root", "4720_068"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/036/AnalysisResults.root", "4720_036"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/053/AnalysisResults.root", "4720_053"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/074/AnalysisResults.root", "4720_074"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/031/AnalysisResults.root", "4720_031"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/065/AnalysisResults.root", "4720_065"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/040/AnalysisResults.root", "4720_040"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/023/AnalysisResults.root", "4720_023"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/008/AnalysisResults.root", "4720_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/075/AnalysisResults.root", "4720_075"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/080/AnalysisResults.root", "4720_080"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/005/AnalysisResults.root", "4720_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/052/AnalysisResults.root", "4720_052"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/018/AnalysisResults.root", "4720_018"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/044/AnalysisResults.root", "4720_044"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/033/AnalysisResults.root", "4720_033"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/011/AnalysisResults.root", "4720_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/058/AnalysisResults.root", "4720_058"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/021/AnalysisResults.root", "4720_021"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/054/AnalysisResults.root", "4720_054"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/046/AnalysisResults.root", "4720_046"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/042/AnalysisResults.root", "4720_042"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/015/AnalysisResults.root", "4720_015"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/012/AnalysisResults.root", "4720_012"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/024/AnalysisResults.root", "4720_024"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/029/AnalysisResults.root", "4720_029"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/063/AnalysisResults.root", "4720_063"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/009/AnalysisResults.root", "4720_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/071/AnalysisResults.root", "4720_071"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/037/AnalysisResults.root", "4720_037"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/035/AnalysisResults.root", "4720_035"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/057/AnalysisResults.root", "4720_057"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/039/AnalysisResults.root", "4720_039"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/077/AnalysisResults.root", "4720_077"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/060/AnalysisResults.root", "4720_060"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/038/AnalysisResults.root", "4720_038"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/020/AnalysisResults.root", "4720_020"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/061/AnalysisResults.root", "4720_061"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/013/AnalysisResults.root", "4720_013"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/051/AnalysisResults.root", "4720_051"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/025/AnalysisResults.root", "4720_025"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/081/AnalysisResults.root", "4720_081"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/014/AnalysisResults.root", "4720_014"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl2/AnalysisResults.root", "4719_rl2"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/007/AnalysisResults.root", "4719_rl1_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/004/AnalysisResults.root", "4719_rl1_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/002/AnalysisResults.root", "4719_rl1_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/006/AnalysisResults.root", "4719_rl1_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/010/AnalysisResults.root", "4719_rl1_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/001/AnalysisResults.root", "4719_rl1_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/003/AnalysisResults.root", "4719_rl1_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/008/AnalysisResults.root", "4719_rl1_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/005/AnalysisResults.root", "4719_rl1_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/011/AnalysisResults.root", "4719_rl1_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/009/AnalysisResults.root", "4719_rl1_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/007/AnalysisResults.root", "4718_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/017/AnalysisResults.root", "4718_017"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/004/AnalysisResults.root", "4718_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/002/AnalysisResults.root", "4718_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/022/AnalysisResults.root", "4718_022"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/016/AnalysisResults.root", "4718_016"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/026/AnalysisResults.root", "4718_026"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/006/AnalysisResults.root", "4718_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/010/AnalysisResults.root", "4718_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/030/AnalysisResults.root", "4718_030"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/028/AnalysisResults.root", "4718_028"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/001/AnalysisResults.root", "4718_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/027/AnalysisResults.root", "4718_027"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/019/AnalysisResults.root", "4718_019"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/003/AnalysisResults.root", "4718_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/031/AnalysisResults.root", "4718_031"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/023/AnalysisResults.root", "4718_023"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/008/AnalysisResults.root", "4718_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/005/AnalysisResults.root", "4718_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/018/AnalysisResults.root", "4718_018"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/011/AnalysisResults.root", "4718_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/021/AnalysisResults.root", "4718_021"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/015/AnalysisResults.root", "4718_015"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/012/AnalysisResults.root", "4718_012"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/024/AnalysisResults.root", "4718_024"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/029/AnalysisResults.root", "4718_029"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/009/AnalysisResults.root", "4718_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/020/AnalysisResults.root", "4718_020"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/013/AnalysisResults.root", "4718_013"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/025/AnalysisResults.root", "4718_025"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/014/AnalysisResults.root", "4718_014"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/007/AnalysisResults.root", "4717_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/017/AnalysisResults.root", "4717_017"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/004/AnalysisResults.root", "4717_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/002/AnalysisResults.root", "4717_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/022/AnalysisResults.root", "4717_022"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/016/AnalysisResults.root", "4717_016"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/026/AnalysisResults.root", "4717_026"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/006/AnalysisResults.root", "4717_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/010/AnalysisResults.root", "4717_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/030/AnalysisResults.root", "4717_030"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/028/AnalysisResults.root", "4717_028"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/001/AnalysisResults.root", "4717_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/027/AnalysisResults.root", "4717_027"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/032/AnalysisResults.root", "4717_032"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/019/AnalysisResults.root", "4717_019"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/034/AnalysisResults.root", "4717_034"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/003/AnalysisResults.root", "4717_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/031/AnalysisResults.root", "4717_031"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/023/AnalysisResults.root", "4717_023"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/008/AnalysisResults.root", "4717_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/005/AnalysisResults.root", "4717_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/018/AnalysisResults.root", "4717_018"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/033/AnalysisResults.root", "4717_033"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/011/AnalysisResults.root", "4717_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/021/AnalysisResults.root", "4717_021"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/015/AnalysisResults.root", "4717_015"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/012/AnalysisResults.root", "4717_012"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/024/AnalysisResults.root", "4717_024"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/029/AnalysisResults.root", "4717_029"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/009/AnalysisResults.root", "4717_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/035/AnalysisResults.root", "4717_035"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/020/AnalysisResults.root", "4717_020"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/013/AnalysisResults.root", "4717_013"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/025/AnalysisResults.root", "4717_025"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/014/AnalysisResults.root", "4717_014"),
                ],
            },
            "12_mix150_nosepsb_sli": {
                "data": [
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/059/AnalysisResults.root", "4720_059"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/064/AnalysisResults.root", "4720_064"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/043/AnalysisResults.root", "4720_043"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/007/AnalysisResults.root", "4720_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/017/AnalysisResults.root", "4720_017"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/062/AnalysisResults.root", "4720_062"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/069/AnalysisResults.root", "4720_069"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/004/AnalysisResults.root", "4720_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/073/AnalysisResults.root", "4720_073"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/002/AnalysisResults.root", "4720_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/041/AnalysisResults.root", "4720_041"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/067/AnalysisResults.root", "4720_067"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/045/AnalysisResults.root", "4720_045"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/047/AnalysisResults.root", "4720_047"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/022/AnalysisResults.root", "4720_022"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/050/AnalysisResults.root", "4720_050"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/016/AnalysisResults.root", "4720_016"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/066/AnalysisResults.root", "4720_066"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/026/AnalysisResults.root", "4720_026"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/079/AnalysisResults.root", "4720_079"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/006/AnalysisResults.root", "4720_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/072/AnalysisResults.root", "4720_072"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/010/AnalysisResults.root", "4720_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/055/AnalysisResults.root", "4720_055"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/030/AnalysisResults.root", "4720_030"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/028/AnalysisResults.root", "4720_028"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/001/AnalysisResults.root", "4720_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/076/AnalysisResults.root", "4720_076"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/027/AnalysisResults.root", "4720_027"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/056/AnalysisResults.root", "4720_056"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/032/AnalysisResults.root", "4720_032"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/048/AnalysisResults.root", "4720_048"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/070/AnalysisResults.root", "4720_070"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/078/AnalysisResults.root", "4720_078"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/019/AnalysisResults.root", "4720_019"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/034/AnalysisResults.root", "4720_034"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/003/AnalysisResults.root", "4720_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/049/AnalysisResults.root", "4720_049"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/068/AnalysisResults.root", "4720_068"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/036/AnalysisResults.root", "4720_036"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/053/AnalysisResults.root", "4720_053"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/074/AnalysisResults.root", "4720_074"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/031/AnalysisResults.root", "4720_031"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/065/AnalysisResults.root", "4720_065"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/040/AnalysisResults.root", "4720_040"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/023/AnalysisResults.root", "4720_023"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/008/AnalysisResults.root", "4720_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/075/AnalysisResults.root", "4720_075"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/080/AnalysisResults.root", "4720_080"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/005/AnalysisResults.root", "4720_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/052/AnalysisResults.root", "4720_052"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/018/AnalysisResults.root", "4720_018"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/044/AnalysisResults.root", "4720_044"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/033/AnalysisResults.root", "4720_033"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/011/AnalysisResults.root", "4720_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/058/AnalysisResults.root", "4720_058"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/021/AnalysisResults.root", "4720_021"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/054/AnalysisResults.root", "4720_054"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/046/AnalysisResults.root", "4720_046"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/042/AnalysisResults.root", "4720_042"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/015/AnalysisResults.root", "4720_015"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/012/AnalysisResults.root", "4720_012"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/024/AnalysisResults.root", "4720_024"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/029/AnalysisResults.root", "4720_029"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/063/AnalysisResults.root", "4720_063"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/009/AnalysisResults.root", "4720_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/071/AnalysisResults.root", "4720_071"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/037/AnalysisResults.root", "4720_037"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/035/AnalysisResults.root", "4720_035"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/057/AnalysisResults.root", "4720_057"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/039/AnalysisResults.root", "4720_039"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/077/AnalysisResults.root", "4720_077"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/060/AnalysisResults.root", "4720_060"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/038/AnalysisResults.root", "4720_038"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/020/AnalysisResults.root", "4720_020"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/061/AnalysisResults.root", "4720_061"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/013/AnalysisResults.root", "4720_013"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/051/AnalysisResults.root", "4720_051"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/025/AnalysisResults.root", "4720_025"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/081/AnalysisResults.root", "4720_081"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4720/Stage_1/014/AnalysisResults.root", "4720_014"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl2/AnalysisResults.root", "4719_rl2"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/007/AnalysisResults.root", "4719_rl1_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/004/AnalysisResults.root", "4719_rl1_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/002/AnalysisResults.root", "4719_rl1_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/006/AnalysisResults.root", "4719_rl1_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/010/AnalysisResults.root", "4719_rl1_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/001/AnalysisResults.root", "4719_rl1_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/003/AnalysisResults.root", "4719_rl1_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/008/AnalysisResults.root", "4719_rl1_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/005/AnalysisResults.root", "4719_rl1_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/011/AnalysisResults.root", "4719_rl1_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4719/rl1/Stage_2/009/AnalysisResults.root", "4719_rl1_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/007/AnalysisResults.root", "4718_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/017/AnalysisResults.root", "4718_017"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/004/AnalysisResults.root", "4718_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/002/AnalysisResults.root", "4718_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/022/AnalysisResults.root", "4718_022"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/016/AnalysisResults.root", "4718_016"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/026/AnalysisResults.root", "4718_026"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/006/AnalysisResults.root", "4718_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/010/AnalysisResults.root", "4718_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/030/AnalysisResults.root", "4718_030"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/028/AnalysisResults.root", "4718_028"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/001/AnalysisResults.root", "4718_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/027/AnalysisResults.root", "4718_027"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/019/AnalysisResults.root", "4718_019"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/003/AnalysisResults.root", "4718_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/031/AnalysisResults.root", "4718_031"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/023/AnalysisResults.root", "4718_023"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/008/AnalysisResults.root", "4718_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/005/AnalysisResults.root", "4718_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/018/AnalysisResults.root", "4718_018"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/011/AnalysisResults.root", "4718_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/021/AnalysisResults.root", "4718_021"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/015/AnalysisResults.root", "4718_015"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/012/AnalysisResults.root", "4718_012"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/024/AnalysisResults.root", "4718_024"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/029/AnalysisResults.root", "4718_029"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/009/AnalysisResults.root", "4718_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/020/AnalysisResults.root", "4718_020"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/013/AnalysisResults.root", "4718_013"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/025/AnalysisResults.root", "4718_025"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4718/Stage_1/014/AnalysisResults.root", "4718_014"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/007/AnalysisResults.root", "4717_007"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/017/AnalysisResults.root", "4717_017"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/004/AnalysisResults.root", "4717_004"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/002/AnalysisResults.root", "4717_002"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/022/AnalysisResults.root", "4717_022"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/016/AnalysisResults.root", "4717_016"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/026/AnalysisResults.root", "4717_026"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/006/AnalysisResults.root", "4717_006"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/010/AnalysisResults.root", "4717_010"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/030/AnalysisResults.root", "4717_030"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/028/AnalysisResults.root", "4717_028"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/001/AnalysisResults.root", "4717_001"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/027/AnalysisResults.root", "4717_027"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/032/AnalysisResults.root", "4717_032"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/019/AnalysisResults.root", "4717_019"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/034/AnalysisResults.root", "4717_034"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/003/AnalysisResults.root", "4717_003"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/031/AnalysisResults.root", "4717_031"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/023/AnalysisResults.root", "4717_023"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/008/AnalysisResults.root", "4717_008"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/005/AnalysisResults.root", "4717_005"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/018/AnalysisResults.root", "4717_018"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/033/AnalysisResults.root", "4717_033"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/011/AnalysisResults.root", "4717_011"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/021/AnalysisResults.root", "4717_021"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/015/AnalysisResults.root", "4717_015"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/012/AnalysisResults.root", "4717_012"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/024/AnalysisResults.root", "4717_024"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/029/AnalysisResults.root", "4717_029"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/009/AnalysisResults.root", "4717_009"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/035/AnalysisResults.root", "4717_035"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/020/AnalysisResults.root", "4717_020"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/013/AnalysisResults.root", "4717_013"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/025/AnalysisResults.root", "4717_025"),
                    ("/data/DPi/12_mix150_nosepsb/data/merge_4717/Stage_1/014/AnalysisResults.root", "4717_014"),
                ]
            }
        }

        luigi.build([
            DPiTask(data_version='10_tree_syst', cleaner='nopc'),
            DPiTask(data_version='10_tree_syst', cleaner='newpc'),
            DPiTask(data_version='10_tree_syst', cleaner='oldpc'),
            DPiTask(data_version='11_mix150', cleaner='nopc'),
            DPiTask(data_version='11_mix150', cleaner='newpc'),
            DPiTask(data_version='11_mix150', cleaner='oldpc'),
            DPiTask(data_version='12_mix150_nosepsb', cleaner='nopc'),
        ], workers=12)

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
            },
            "11_mix150": {
                "data": [
                    ("/data/DK/11_mix150/data/sgn/AnalysisResults_4694.root", "4694"),
                    ("/data/DK/11_mix150/data/sgn/AnalysisResults_4695.root", "4695"),
                    ("/data/DK/11_mix150/data/sgn/AnalysisResults_4696.root", "4696"),
                    ("/data/DK/11_mix150/data/sgn/AnalysisResults_4697.root", "4697"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/001/AnalysisResults.root", "sbl_4702_001"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/002/AnalysisResults.root", "sbl_4702_002"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/003/AnalysisResults.root", "sbl_4702_003"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/004/AnalysisResults.root", "sbl_4702_004"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/005/AnalysisResults.root", "sbl_4702_005"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/006/AnalysisResults.root", "sbl_4702_006"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/007/AnalysisResults.root", "sbl_4702_007"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/008/AnalysisResults.root", "sbl_4702_008"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/009/AnalysisResults.root", "sbl_4702_009"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/010/AnalysisResults.root", "sbl_4702_010"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/011/AnalysisResults.root", "sbl_4702_011"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/012/AnalysisResults.root", "sbl_4702_012"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/013/AnalysisResults.root", "sbl_4702_013"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/014/AnalysisResults.root", "sbl_4702_014"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/015/AnalysisResults.root", "sbl_4702_015"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/016/AnalysisResults.root", "sbl_4702_016"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/017/AnalysisResults.root", "sbl_4702_017"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/018/AnalysisResults.root", "sbl_4702_018"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/019/AnalysisResults.root", "sbl_4702_019"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/020/AnalysisResults.root", "sbl_4702_020"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/021/AnalysisResults.root", "sbl_4702_021"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/022/AnalysisResults.root", "sbl_4702_022"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/023/AnalysisResults.root", "sbl_4702_023"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/024/AnalysisResults.root", "sbl_4702_024"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/025/AnalysisResults.root", "sbl_4702_025"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/026/AnalysisResults.root", "sbl_4702_026"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/027/AnalysisResults.root", "sbl_4702_027"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/028/AnalysisResults.root", "sbl_4702_028"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/029/AnalysisResults.root", "sbl_4702_029"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/030/AnalysisResults.root", "sbl_4702_030"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/031/AnalysisResults.root", "sbl_4702_031"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/032/AnalysisResults.root", "sbl_4702_032"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/033/AnalysisResults.root", "sbl_4702_033"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/034/AnalysisResults.root", "sbl_4702_034"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/035/AnalysisResults.root", "sbl_4702_035"),
                    ("/data/DK/11_mix150/data/sbl/merge_4702/Stage_1/036/AnalysisResults.root", "sbl_4702_036"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/001/AnalysisResults.root", "sbl_4703_001"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/002/AnalysisResults.root", "sbl_4703_002"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/003/AnalysisResults.root", "sbl_4703_003"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/004/AnalysisResults.root", "sbl_4703_004"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/005/AnalysisResults.root", "sbl_4703_005"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/006/AnalysisResults.root", "sbl_4703_006"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/007/AnalysisResults.root", "sbl_4703_007"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/008/AnalysisResults.root", "sbl_4703_008"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/009/AnalysisResults.root", "sbl_4703_009"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/010/AnalysisResults.root", "sbl_4703_010"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/011/AnalysisResults.root", "sbl_4703_011"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/012/AnalysisResults.root", "sbl_4703_012"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/013/AnalysisResults.root", "sbl_4703_013"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/014/AnalysisResults.root", "sbl_4703_014"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/015/AnalysisResults.root", "sbl_4703_015"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/016/AnalysisResults.root", "sbl_4703_016"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/017/AnalysisResults.root", "sbl_4703_017"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/018/AnalysisResults.root", "sbl_4703_018"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/019/AnalysisResults.root", "sbl_4703_019"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/020/AnalysisResults.root", "sbl_4703_020"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/021/AnalysisResults.root", "sbl_4703_021"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/022/AnalysisResults.root", "sbl_4703_022"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/023/AnalysisResults.root", "sbl_4703_023"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/024/AnalysisResults.root", "sbl_4703_024"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/025/AnalysisResults.root", "sbl_4703_025"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/026/AnalysisResults.root", "sbl_4703_026"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/027/AnalysisResults.root", "sbl_4703_027"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/028/AnalysisResults.root", "sbl_4703_028"),
                    ("/data/DK/11_mix150/data/sbl/merge_4703/Stage_1/029/AnalysisResults.root", "sbl_4703_029"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/001/AnalysisResults.root", "sbl_4705_001"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/002/AnalysisResults.root", "sbl_4705_002"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/003/AnalysisResults.root", "sbl_4705_003"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/004/AnalysisResults.root", "sbl_4705_004"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/005/AnalysisResults.root", "sbl_4705_005"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/006/AnalysisResults.root", "sbl_4705_006"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/007/AnalysisResults.root", "sbl_4705_007"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/008/AnalysisResults.root", "sbl_4705_008"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/009/AnalysisResults.root", "sbl_4705_009"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/010/AnalysisResults.root", "sbl_4705_010"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/011/AnalysisResults.root", "sbl_4705_011"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/012/AnalysisResults.root", "sbl_4705_012"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/013/AnalysisResults.root", "sbl_4705_013"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/014/AnalysisResults.root", "sbl_4705_014"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/015/AnalysisResults.root", "sbl_4705_015"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/016/AnalysisResults.root", "sbl_4705_016"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/017/AnalysisResults.root", "sbl_4705_017"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/018/AnalysisResults.root", "sbl_4705_018"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/019/AnalysisResults.root", "sbl_4705_019"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/020/AnalysisResults.root", "sbl_4705_020"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/021/AnalysisResults.root", "sbl_4705_021"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/022/AnalysisResults.root", "sbl_4705_022"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/023/AnalysisResults.root", "sbl_4705_023"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/024/AnalysisResults.root", "sbl_4705_024"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/025/AnalysisResults.root", "sbl_4705_025"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/026/AnalysisResults.root", "sbl_4705_026"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/027/AnalysisResults.root", "sbl_4705_027"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/028/AnalysisResults.root", "sbl_4705_028"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/029/AnalysisResults.root", "sbl_4705_029"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/030/AnalysisResults.root", "sbl_4705_030"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/031/AnalysisResults.root", "sbl_4705_031"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/032/AnalysisResults.root", "sbl_4705_032"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/033/AnalysisResults.root", "sbl_4705_033"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/034/AnalysisResults.root", "sbl_4705_034"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/035/AnalysisResults.root", "sbl_4705_035"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/036/AnalysisResults.root", "sbl_4705_036"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/037/AnalysisResults.root", "sbl_4705_037"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/038/AnalysisResults.root", "sbl_4705_038"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/039/AnalysisResults.root", "sbl_4705_039"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/040/AnalysisResults.root", "sbl_4705_040"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/041/AnalysisResults.root", "sbl_4705_041"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/042/AnalysisResults.root", "sbl_4705_042"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/043/AnalysisResults.root", "sbl_4705_043"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/044/AnalysisResults.root", "sbl_4705_044"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/045/AnalysisResults.root", "sbl_4705_045"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/046/AnalysisResults.root", "sbl_4705_046"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/047/AnalysisResults.root", "sbl_4705_047"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/048/AnalysisResults.root", "sbl_4705_048"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/049/AnalysisResults.root", "sbl_4705_049"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/050/AnalysisResults.root", "sbl_4705_050"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/051/AnalysisResults.root", "sbl_4705_051"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/052/AnalysisResults.root", "sbl_4705_052"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/053/AnalysisResults.root", "sbl_4705_053"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/054/AnalysisResults.root", "sbl_4705_054"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/055/AnalysisResults.root", "sbl_4705_055"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/056/AnalysisResults.root", "sbl_4705_056"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/057/AnalysisResults.root", "sbl_4705_057"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/058/AnalysisResults.root", "sbl_4705_058"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/059/AnalysisResults.root", "sbl_4705_059"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/060/AnalysisResults.root", "sbl_4705_060"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/061/AnalysisResults.root", "sbl_4705_061"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/062/AnalysisResults.root", "sbl_4705_062"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/063/AnalysisResults.root", "sbl_4705_063"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/064/AnalysisResults.root", "sbl_4705_064"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/065/AnalysisResults.root", "sbl_4705_065"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/066/AnalysisResults.root", "sbl_4705_066"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/067/AnalysisResults.root", "sbl_4705_067"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/068/AnalysisResults.root", "sbl_4705_068"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/069/AnalysisResults.root", "sbl_4705_069"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/070/AnalysisResults.root", "sbl_4705_070"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/071/AnalysisResults.root", "sbl_4705_071"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/072/AnalysisResults.root", "sbl_4705_072"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/073/AnalysisResults.root", "sbl_4705_073"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/074/AnalysisResults.root", "sbl_4705_074"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/075/AnalysisResults.root", "sbl_4705_075"),
                    ("/data/DK/11_mix150/data/sbl/merge_4705/Stage_1/076/AnalysisResults.root", "sbl_4705_076"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/001/AnalysisResults.root", "sbr_4702_001"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/002/AnalysisResults.root", "sbr_4702_002"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/003/AnalysisResults.root", "sbr_4702_003"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/004/AnalysisResults.root", "sbr_4702_004"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/005/AnalysisResults.root", "sbr_4702_005"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/006/AnalysisResults.root", "sbr_4702_006"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/007/AnalysisResults.root", "sbr_4702_007"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/008/AnalysisResults.root", "sbr_4702_008"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/009/AnalysisResults.root", "sbr_4702_009"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/010/AnalysisResults.root", "sbr_4702_010"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/011/AnalysisResults.root", "sbr_4702_011"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/012/AnalysisResults.root", "sbr_4702_012"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/013/AnalysisResults.root", "sbr_4702_013"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/014/AnalysisResults.root", "sbr_4702_014"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/015/AnalysisResults.root", "sbr_4702_015"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/016/AnalysisResults.root", "sbr_4702_016"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/017/AnalysisResults.root", "sbr_4702_017"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/018/AnalysisResults.root", "sbr_4702_018"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/019/AnalysisResults.root", "sbr_4702_019"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/020/AnalysisResults.root", "sbr_4702_020"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/021/AnalysisResults.root", "sbr_4702_021"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/022/AnalysisResults.root", "sbr_4702_022"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/023/AnalysisResults.root", "sbr_4702_023"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/024/AnalysisResults.root", "sbr_4702_024"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/025/AnalysisResults.root", "sbr_4702_025"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/026/AnalysisResults.root", "sbr_4702_026"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/027/AnalysisResults.root", "sbr_4702_027"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/028/AnalysisResults.root", "sbr_4702_028"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/029/AnalysisResults.root", "sbr_4702_029"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/030/AnalysisResults.root", "sbr_4702_030"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/031/AnalysisResults.root", "sbr_4702_031"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/032/AnalysisResults.root", "sbr_4702_032"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/033/AnalysisResults.root", "sbr_4702_033"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/034/AnalysisResults.root", "sbr_4702_034"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/035/AnalysisResults.root", "sbr_4702_035"),
                    ("/data/DK/11_mix150/data/sbr/merge_4702/Stage_1/036/AnalysisResults.root", "sbr_4702_036"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/001/AnalysisResults.root", "sbr_4703_001"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/002/AnalysisResults.root", "sbr_4703_002"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/003/AnalysisResults.root", "sbr_4703_003"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/004/AnalysisResults.root", "sbr_4703_004"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/005/AnalysisResults.root", "sbr_4703_005"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/006/AnalysisResults.root", "sbr_4703_006"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/007/AnalysisResults.root", "sbr_4703_007"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/008/AnalysisResults.root", "sbr_4703_008"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/009/AnalysisResults.root", "sbr_4703_009"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/010/AnalysisResults.root", "sbr_4703_010"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/011/AnalysisResults.root", "sbr_4703_011"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/012/AnalysisResults.root", "sbr_4703_012"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/013/AnalysisResults.root", "sbr_4703_013"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/014/AnalysisResults.root", "sbr_4703_014"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/015/AnalysisResults.root", "sbr_4703_015"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/016/AnalysisResults.root", "sbr_4703_016"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/017/AnalysisResults.root", "sbr_4703_017"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/018/AnalysisResults.root", "sbr_4703_018"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/019/AnalysisResults.root", "sbr_4703_019"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/020/AnalysisResults.root", "sbr_4703_020"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/021/AnalysisResults.root", "sbr_4703_021"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/022/AnalysisResults.root", "sbr_4703_022"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/023/AnalysisResults.root", "sbr_4703_023"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/024/AnalysisResults.root", "sbr_4703_024"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/025/AnalysisResults.root", "sbr_4703_025"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/026/AnalysisResults.root", "sbr_4703_026"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/027/AnalysisResults.root", "sbr_4703_027"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/028/AnalysisResults.root", "sbr_4703_028"),
                    ("/data/DK/11_mix150/data/sbr/merge_4703/Stage_1/029/AnalysisResults.root", "sbr_4703_029"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/001/AnalysisResults.root", "sbr_4705_001"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/002/AnalysisResults.root", "sbr_4705_002"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/003/AnalysisResults.root", "sbr_4705_003"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/004/AnalysisResults.root", "sbr_4705_004"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/005/AnalysisResults.root", "sbr_4705_005"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/006/AnalysisResults.root", "sbr_4705_006"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/007/AnalysisResults.root", "sbr_4705_007"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/008/AnalysisResults.root", "sbr_4705_008"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/009/AnalysisResults.root", "sbr_4705_009"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/010/AnalysisResults.root", "sbr_4705_010"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/011/AnalysisResults.root", "sbr_4705_011"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/012/AnalysisResults.root", "sbr_4705_012"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/013/AnalysisResults.root", "sbr_4705_013"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/014/AnalysisResults.root", "sbr_4705_014"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/015/AnalysisResults.root", "sbr_4705_015"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/016/AnalysisResults.root", "sbr_4705_016"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/017/AnalysisResults.root", "sbr_4705_017"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/018/AnalysisResults.root", "sbr_4705_018"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/019/AnalysisResults.root", "sbr_4705_019"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/020/AnalysisResults.root", "sbr_4705_020"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/021/AnalysisResults.root", "sbr_4705_021"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/022/AnalysisResults.root", "sbr_4705_022"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/023/AnalysisResults.root", "sbr_4705_023"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/024/AnalysisResults.root", "sbr_4705_024"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/025/AnalysisResults.root", "sbr_4705_025"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/026/AnalysisResults.root", "sbr_4705_026"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/027/AnalysisResults.root", "sbr_4705_027"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/028/AnalysisResults.root", "sbr_4705_028"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/029/AnalysisResults.root", "sbr_4705_029"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/030/AnalysisResults.root", "sbr_4705_030"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/031/AnalysisResults.root", "sbr_4705_031"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/032/AnalysisResults.root", "sbr_4705_032"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/033/AnalysisResults.root", "sbr_4705_033"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/034/AnalysisResults.root", "sbr_4705_034"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/035/AnalysisResults.root", "sbr_4705_035"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/036/AnalysisResults.root", "sbr_4705_036"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/037/AnalysisResults.root", "sbr_4705_037"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/038/AnalysisResults.root", "sbr_4705_038"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/039/AnalysisResults.root", "sbr_4705_039"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/040/AnalysisResults.root", "sbr_4705_040"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/041/AnalysisResults.root", "sbr_4705_041"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/042/AnalysisResults.root", "sbr_4705_042"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/043/AnalysisResults.root", "sbr_4705_043"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/044/AnalysisResults.root", "sbr_4705_044"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/045/AnalysisResults.root", "sbr_4705_045"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/046/AnalysisResults.root", "sbr_4705_046"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/047/AnalysisResults.root", "sbr_4705_047"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/048/AnalysisResults.root", "sbr_4705_048"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/049/AnalysisResults.root", "sbr_4705_049"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/050/AnalysisResults.root", "sbr_4705_050"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/051/AnalysisResults.root", "sbr_4705_051"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/052/AnalysisResults.root", "sbr_4705_052"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/053/AnalysisResults.root", "sbr_4705_053"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/054/AnalysisResults.root", "sbr_4705_054"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/055/AnalysisResults.root", "sbr_4705_055"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/056/AnalysisResults.root", "sbr_4705_056"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/057/AnalysisResults.root", "sbr_4705_057"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/058/AnalysisResults.root", "sbr_4705_058"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/059/AnalysisResults.root", "sbr_4705_059"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/060/AnalysisResults.root", "sbr_4705_060"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/061/AnalysisResults.root", "sbr_4705_061"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/062/AnalysisResults.root", "sbr_4705_062"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/063/AnalysisResults.root", "sbr_4705_063"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/064/AnalysisResults.root", "sbr_4705_064"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/065/AnalysisResults.root", "sbr_4705_065"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/066/AnalysisResults.root", "sbr_4705_066"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/067/AnalysisResults.root", "sbr_4705_067"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/068/AnalysisResults.root", "sbr_4705_068"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/069/AnalysisResults.root", "sbr_4705_069"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/070/AnalysisResults.root", "sbr_4705_070"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/071/AnalysisResults.root", "sbr_4705_071"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/072/AnalysisResults.root", "sbr_4705_072"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/073/AnalysisResults.root", "sbr_4705_073"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/074/AnalysisResults.root", "sbr_4705_074"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/075/AnalysisResults.root", "sbr_4705_075"),
                    ("/data/DK/11_mix150/data/sbr/merge_4705/Stage_1/076/AnalysisResults.root", "sbr_4705_076"),
                ],
                "mcgptruthreco": [
                    ("/data/DK/11_mix150/mcgptruthreco/AnalysisResults_4199.root", "4199"),
                    ("/data/DK/11_mix150/mcgptruthreco/AnalysisResults_4200.root", "4200"),
                    ("/data/DK/11_mix150/mcgptruthreco/AnalysisResults_4201.root", "4201"),
                ],
                "mchftruthreco": [
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/001/AnalysisResults.root", "4205_001"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/002/AnalysisResults.root", "4205_002"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/003/AnalysisResults.root", "4205_003"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/004/AnalysisResults.root", "4205_004"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/005/AnalysisResults.root", "4205_005"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/006/AnalysisResults.root", "4205_006"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/007/AnalysisResults.root", "4205_007"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/008/AnalysisResults.root", "4205_008"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/009/AnalysisResults.root", "4205_009"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/010/AnalysisResults.root", "4205_010"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/011/AnalysisResults.root", "4205_011"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/012/AnalysisResults.root", "4205_012"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/013/AnalysisResults.root", "4205_013"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/014/AnalysisResults.root", "4205_014"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/015/AnalysisResults.root", "4205_015"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/016/AnalysisResults.root", "4205_016"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/017/AnalysisResults.root", "4205_017"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/018/AnalysisResults.root", "4205_018"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/019/AnalysisResults.root", "4205_019"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/020/AnalysisResults.root", "4205_020"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/021/AnalysisResults.root", "4205_021"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/022/AnalysisResults.root", "4205_022"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/023/AnalysisResults.root", "4205_023"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/024/AnalysisResults.root", "4205_024"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/025/AnalysisResults.root", "4205_025"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/026/AnalysisResults.root", "4205_026"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/027/AnalysisResults.root", "4205_027"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/028/AnalysisResults.root", "4205_028"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/029/AnalysisResults.root", "4205_029"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/030/AnalysisResults.root", "4205_030"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/031/AnalysisResults.root", "4205_031"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/032/AnalysisResults.root", "4205_032"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/033/AnalysisResults.root", "4205_033"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/034/AnalysisResults.root", "4205_034"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/035/AnalysisResults.root", "4205_035"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/036/AnalysisResults.root", "4205_036"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/037/AnalysisResults.root", "4205_037"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/038/AnalysisResults.root", "4205_038"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/039/AnalysisResults.root", "4205_039"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/040/AnalysisResults.root", "4205_040"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/041/AnalysisResults.root", "4205_041"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/042/AnalysisResults.root", "4205_042"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/043/AnalysisResults.root", "4205_043"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/044/AnalysisResults.root", "4205_044"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4205/Stage_1/045/AnalysisResults.root", "4205_045"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/001/AnalysisResults.root", "4206_001"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/002/AnalysisResults.root", "4206_002"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/003/AnalysisResults.root", "4206_003"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/004/AnalysisResults.root", "4206_004"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/005/AnalysisResults.root", "4206_005"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/006/AnalysisResults.root", "4206_006"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/007/AnalysisResults.root", "4206_007"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/008/AnalysisResults.root", "4206_008"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/009/AnalysisResults.root", "4206_009"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/010/AnalysisResults.root", "4206_010"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/011/AnalysisResults.root", "4206_011"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/012/AnalysisResults.root", "4206_012"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/013/AnalysisResults.root", "4206_013"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/014/AnalysisResults.root", "4206_014"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/015/AnalysisResults.root", "4206_015"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/016/AnalysisResults.root", "4206_016"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/017/AnalysisResults.root", "4206_017"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/018/AnalysisResults.root", "4206_018"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/019/AnalysisResults.root", "4206_019"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/020/AnalysisResults.root", "4206_020"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/021/AnalysisResults.root", "4206_021"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/022/AnalysisResults.root", "4206_022"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/023/AnalysisResults.root", "4206_023"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/024/AnalysisResults.root", "4206_024"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/025/AnalysisResults.root", "4206_025"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/026/AnalysisResults.root", "4206_026"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/027/AnalysisResults.root", "4206_027"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/028/AnalysisResults.root", "4206_028"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/029/AnalysisResults.root", "4206_029"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/030/AnalysisResults.root", "4206_030"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/031/AnalysisResults.root", "4206_031"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/032/AnalysisResults.root", "4206_032"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/033/AnalysisResults.root", "4206_033"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/034/AnalysisResults.root", "4206_034"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/035/AnalysisResults.root", "4206_035"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/036/AnalysisResults.root", "4206_036"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/037/AnalysisResults.root", "4206_037"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/038/AnalysisResults.root", "4206_038"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/039/AnalysisResults.root", "4206_039"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/040/AnalysisResults.root", "4206_040"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/041/AnalysisResults.root", "4206_041"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/042/AnalysisResults.root", "4206_042"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/043/AnalysisResults.root", "4206_043"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/044/AnalysisResults.root", "4206_044"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4206/Stage_1/045/AnalysisResults.root", "4206_045"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/001/AnalysisResults.root", "4207_001"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/002/AnalysisResults.root", "4207_002"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/003/AnalysisResults.root", "4207_003"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/004/AnalysisResults.root", "4207_004"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/005/AnalysisResults.root", "4207_005"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/006/AnalysisResults.root", "4207_006"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/007/AnalysisResults.root", "4207_007"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/008/AnalysisResults.root", "4207_008"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/009/AnalysisResults.root", "4207_009"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/010/AnalysisResults.root", "4207_010"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/011/AnalysisResults.root", "4207_011"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/012/AnalysisResults.root", "4207_012"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/013/AnalysisResults.root", "4207_013"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/014/AnalysisResults.root", "4207_014"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/015/AnalysisResults.root", "4207_015"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/016/AnalysisResults.root", "4207_016"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/017/AnalysisResults.root", "4207_017"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/018/AnalysisResults.root", "4207_018"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/019/AnalysisResults.root", "4207_019"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/020/AnalysisResults.root", "4207_020"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/021/AnalysisResults.root", "4207_021"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/022/AnalysisResults.root", "4207_022"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/023/AnalysisResults.root", "4207_023"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/024/AnalysisResults.root", "4207_024"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/025/AnalysisResults.root", "4207_025"),
                    ("/data/DK/11_mix150/mchftruthreco/merge_4207/Stage_1/026/AnalysisResults.root", "4207_026"),
                ],
            },
            "12_mix150_nosepsb": {
                "data": [
                    ("/data/DK/12_mix150_nosepsb/data/AnalysisResults_4712.root", "4712"),
                    ("/data/DK/12_mix150_nosepsb/data/AnalysisResults_4713.root", "4713"),
                    ("/data/DK/12_mix150_nosepsb/data/AnalysisResults_4714.root", "4714"),
                    # ("/data/DK/12_mix150_nosepsb/data/AnalysisResults_4721.root", "4721"),
                ]
            }
        }

        luigi.build([
            # DKTask(data_version='10_tree_syst', cleaner='nopc'),
            # DKTask(data_version='10_tree_syst', cleaner='newpc'),
            # DKTask(data_version='10_tree_syst', cleaner='oldpc'),
            # DKTask(data_version='11_mix150', cleaner='nopc'),
            # DKTask(data_version='11_mix150', cleaner='newpc'),
            # DKTask(data_version='11_mix150', cleaner='oldpc'),
            DKTask(data_version='12_mix150_nosepsb', cleaner='nopc'),
        ], workers=12)
