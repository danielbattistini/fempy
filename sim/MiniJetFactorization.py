import os
import sys
import luigi

from ROOT import EColor

from MJFactorizationPlots import MJFactorizatonPlots


luigi.interface.core.log_level="INFO"

baseDir = "/home/daniel/an/mjfactorization/sim"

def MakeSuffix(
    pdg1=None,
    pdg2=None,
    alignment=None,
    computeKStarMode=None,
    minMult=None,
    maxMult=None,
    minMajorSpheri=None,
    maxMajorSpheri=None,
    minMinorSpheri=None,
    maxMinorSpheri=None,
    minLambdaSpheriMax=None,
    maxLambdaSpheriMax=None,
    minLambdaSpheriMid=None,
    maxLambdaSpheriMid=None,
    minLambdaSpheriMin=None,
    maxLambdaSpheriMin=None,
    minPrincipalAxisTheta=None,
    maxPrincipalAxisTheta=None,
    seed=None,

):
    suffix = ''
    if (pdg1 == 2212):
        suffix += 'p'
    elif (pdg1 == 22):
        suffix += 'gamma'
    elif (pdg1 == 211):
        suffix += 'pi'
    elif (pdg1 == 411):
        suffix += 'D'
    else:
        print(f'{pdg1} not implemented. Exit')
        sys.exit()
        
    if (pdg2 == 2212):
        suffix += 'p'
    elif (pdg2 == -2212):
        suffix += 'pbar'
    elif (pdg2 == 3122):
        suffix += 'L'
    elif (pdg2 == 22):
        suffix += 'gamma'
    elif (pdg2 == 211):
        suffix += 'pi'
    elif (pdg2 == -211):
        suffix += 'pibar'
    else:
        print(f'{pdg2} not implemented. Exit')
        sys.exit()

    suffix += f'_align{alignment}'
    suffix += f'_kStar{computeKStarMode}'

    # multiplicity
    if minMult is not None:
        suffix += f'_mult{minMult:.0f}'
    if maxMult is not None:
        suffix += f'-{maxMult:.0f}'

    # major and minor sphericity
    if minMajorSpheri is not None:
        suffix += f'_majorSphi{minMajorSpheri:.1f}'
    if maxMajorSpheri is not None:
        suffix += f'-{maxMajorSpheri:.1f}'
    if minMinorSpheri is not None:
        suffix += f'_minorSphi{minMinorSpheri:.1f}'
    if maxMinorSpheri is not None:
        suffix += f'-{maxMinorSpheri:.1f}'

    # sphericity eigenvalues
    if minLambdaSpheriMax is not None:
        suffix += f'_LambdaSpheriMax{minLambdaSpheriMax:.1f}'
    if maxLambdaSpheriMax is not None:
        suffix += f'-{maxLambdaSpheriMax:.1f}'
    if minLambdaSpheriMid is not None:
        suffix += f'_LambdaSpheriMid{minLambdaSpheriMid:.1f}'
    if maxLambdaSpheriMid is not None:
        suffix += f'-{maxLambdaSpheriMid:.1f}'
    if minLambdaSpheriMin is not None:
        suffix += f'_LambdaSpheriMin{minLambdaSpheriMin:.1f}'
    if maxLambdaSpheriMin is not None:
        suffix += f'-{maxLambdaSpheriMin:.1f}'

    # theta direction of jet
    if minPrincipalAxisTheta is not None:
        suffix += f'_pat{minPrincipalAxisTheta:.1f}'
    if maxPrincipalAxisTheta is not None:
        suffix += f'-{maxPrincipalAxisTheta:.1f}'
    
    if seed != None:
        suffix += f'_{seed}'

    return suffix


class SubSimulationTask(luigi.Task):
    nEvents=luigi.IntParameter()
    pdg1=luigi.IntParameter()
    pdg2=luigi.IntParameter()
    seed=luigi.IntParameter()
    alignment=luigi.Parameter()
    computeKStarMode=luigi.Parameter()
    minMult=luigi.FloatParameter()
    maxMult=luigi.FloatParameter()
    minMajorSpheri=luigi.FloatParameter()
    maxMajorSpheri=luigi.FloatParameter()
    minMinorSpheri=luigi.FloatParameter()
    maxMinorSpheri=luigi.FloatParameter()
    minLambdaSpheriMax=luigi.FloatParameter()
    maxLambdaSpheriMax=luigi.FloatParameter()
    minLambdaSpheriMid=luigi.FloatParameter()
    maxLambdaSpheriMid=luigi.FloatParameter()
    minLambdaSpheriMin=luigi.FloatParameter()
    maxLambdaSpheriMin=luigi.FloatParameter()
    minPrincipalAxisTheta=luigi.FloatParameter()
    maxPrincipalAxisTheta=luigi.FloatParameter()

    def __init__(self, *args, **kwargs):
        super(SubSimulationTask, self).__init__(*args, **kwargs)

        self.suffix = MakeSuffix(
            pdg1=self.pdg1,
            pdg2=self.pdg2,
            seed=self.seed,
            alignment=self.alignment,
            computeKStarMode=self.computeKStarMode,
            minMult=self.minMult,
            maxMult=self.maxMult,
            minMajorSpheri=self.minMajorSpheri,
            maxMajorSpheri=self.maxMajorSpheri,
            minMinorSpheri=self.minMinorSpheri,
            maxMinorSpheri=self.maxMinorSpheri,
            minLambdaSpheriMax=self.minLambdaSpheriMax,
            maxLambdaSpheriMax=self.maxLambdaSpheriMax,
            minLambdaSpheriMid=self.minLambdaSpheriMid,
            maxLambdaSpheriMid=self.maxLambdaSpheriMid,
            minLambdaSpheriMin=self.minLambdaSpheriMin,
            maxLambdaSpheriMin=self.maxLambdaSpheriMin,
            minPrincipalAxisTheta=self.minPrincipalAxisTheta,
            maxPrincipalAxisTheta=self.maxPrincipalAxisTheta,
        )
        
        self.outFileNameRoot = f'{baseDir}/part/AnalysisResults'
        if self.suffix != '':
            self.outFileNameRoot += f'_{self.suffix}.root'
        else:
            self.outFileNameRoot += '.root'

    
    def run(self):

        params = f'{self.nEvents}, '\
                f'{self.pdg1}, '\
                f'{self.pdg2}, '\
                f'{self.seed}, '\
                f'k{self.alignment if self.alignment is not None else "No"}, '\
                f'k{self.computeKStarMode if self.computeKStarMode is not None else "Custom"}, '\
                f'{self.minMult if self.minMult is not None else "0"}, '\
                f'{self.maxMult if self.maxMult is not None else "1.e10"}, '\
                f'{self.minMajorSpheri if self.minMajorSpheri is not None else "0"}, '\
                f'{self.maxMajorSpheri if self.maxMajorSpheri is not None else "1"}, '\
                f'{self.minMinorSpheri if self.minMinorSpheri is not None else "0"}, '\
                f'{self.maxMinorSpheri if self.maxMinorSpheri is not None else "1"}, '\
                f'{self.minLambdaSpheriMax if self.minLambdaSpheriMax is not None else "0"}, '\
                f'{self.maxLambdaSpheriMax if self.maxLambdaSpheriMax is not None else "1.e10"}, '\
                f'{self.minLambdaSpheriMid if self.minLambdaSpheriMid is not None else "0"}, '\
                f'{self.maxLambdaSpheriMid if self.maxLambdaSpheriMid is not None else "1.e10"}, '\
                f'{self.minLambdaSpheriMin if self.minLambdaSpheriMin is not None else "0"}, '\
                f'{self.maxLambdaSpheriMin if self.maxLambdaSpheriMin is not None else "1.e10"}, '\
                f'{self.minPrincipalAxisTheta if self.minPrincipalAxisTheta is not None else "-3.14"}, '\
                f'{self.maxPrincipalAxisTheta if self.maxPrincipalAxisTheta is not None else "+3.14"}, '\
                f'\\\"{self.outFileNameRoot}\\\"'
        
        command = f'root -l -b -q -e ".x MiniJetFactorization.cc+({params})" > {self.outFileNameRoot[:-5]}.log'
        print(f'Running the macro: MiniJetFactorization.cc+({params})')
        os.system(command)

    def output(self):
        return luigi.LocalTarget(self.outFileNameRoot)


class SimulationTask(luigi.Task):
    nEvents=luigi.IntParameter(default=1001)
    pdg1=luigi.IntParameter(default=None)
    pdg2=luigi.IntParameter(default=None)
    seedFrom=luigi.IntParameter(default=1)
    seedTo=luigi.IntParameter(default=16)
    alignment=luigi.Parameter(default='No')
    computeKStarMode=luigi.Parameter(default='Custom')
    minMult=luigi.FloatParameter(default=None)
    maxMult=luigi.FloatParameter(default=None)
    minMajorSpheri=luigi.FloatParameter(default=None)
    maxMajorSpheri=luigi.FloatParameter(default=None)
    minMinorSpheri=luigi.FloatParameter(default=None)
    maxMinorSpheri=luigi.FloatParameter(default=None)
    minLambdaSpheriMax=luigi.FloatParameter(default=None)
    maxLambdaSpheriMax=luigi.FloatParameter(default=None)
    minLambdaSpheriMid=luigi.FloatParameter(default=None)
    maxLambdaSpheriMid=luigi.FloatParameter(default=None)
    minLambdaSpheriMin=luigi.FloatParameter(default=None)
    maxLambdaSpheriMin=luigi.FloatParameter(default=None)
    minPrincipalAxisTheta=luigi.FloatParameter(default=None)
    maxPrincipalAxisTheta=luigi.FloatParameter(default=None)

    def __init__(self, *args, **kwargs):
        super(SimulationTask, self).__init__(*args, **kwargs)

        suffixMerged = MakeSuffix(
                pdg1=self.pdg1,
                pdg2=self.pdg2,
                alignment=self.alignment,
                computeKStarMode=self.computeKStarMode,
                minMult=self.minMult,
                maxMult=self.maxMult,
                minMajorSpheri=self.minMajorSpheri,
                maxMajorSpheri=self.maxMajorSpheri,
                minMinorSpheri=self.minMinorSpheri,
                maxMinorSpheri=self.maxMinorSpheri,
                minLambdaSpheriMax=self.minLambdaSpheriMax,
                maxLambdaSpheriMax=self.maxLambdaSpheriMax,
                minLambdaSpheriMid=self.minLambdaSpheriMid,
                maxLambdaSpheriMid=self.maxLambdaSpheriMid,
                minLambdaSpheriMin=self.minLambdaSpheriMin,
                maxLambdaSpheriMin=self.maxLambdaSpheriMin,
                minPrincipalAxisTheta=self.minPrincipalAxisTheta,
                maxPrincipalAxisTheta=self.maxPrincipalAxisTheta,
            ) 
        
        self.mergedFile = f'{baseDir}/AnalysisResults_{suffixMerged}.root'

        suffixes = [
            MakeSuffix(
                pdg1=self.pdg1,
                pdg2=self.pdg2,
                seed=seed,
                alignment=self.alignment,
                computeKStarMode=self.computeKStarMode,
                minMult=self.minMult,
                maxMult=self.maxMult,
                minMajorSpheri=self.minMajorSpheri,
                maxMajorSpheri=self.maxMajorSpheri,
                minMinorSpheri=self.minMinorSpheri,
                maxMinorSpheri=self.maxMinorSpheri,
                minLambdaSpheriMax=self.minLambdaSpheriMax,
                maxLambdaSpheriMax=self.maxLambdaSpheriMax,
                minLambdaSpheriMid=self.minLambdaSpheriMid,
                maxLambdaSpheriMid=self.maxLambdaSpheriMid,
                minLambdaSpheriMin=self.minLambdaSpheriMin,
                maxLambdaSpheriMin=self.maxLambdaSpheriMin,
                minPrincipalAxisTheta=self.minPrincipalAxisTheta,
                maxPrincipalAxisTheta=self.maxPrincipalAxisTheta,
            ) for seed in range(self.seedFrom, self.seedTo+1)]
        
        
        self.filesToMerge = [f'{baseDir}/part/AnalysisResults_{suffix}.root' for suffix in suffixes]

    def run(self):
        command = 'hadd'
        command = f'hadd -f {self.mergedFile} {" ".join(self.filesToMerge)}'

        print(command)
        os.system(command)

    def output(self):
        return [luigi.LocalTarget(self.mergedFile)] + [luigi.LocalTarget(f) for f in self.filesToMerge]



    def requires(self):
        print("Compiling the macro...")
        os.system('root -l -b -q -e ".L MiniJetFactorization.cc+"')

        yield [
            SubSimulationTask(
                nEvents=self.nEvents,
                pdg1=self.pdg1,
                pdg2=self.pdg2,
                seed=seed,
                alignment=self.alignment,
                computeKStarMode=self.computeKStarMode,
                minMult=self.minMult,
                maxMult=self.maxMult,
                minMajorSpheri=self.minMajorSpheri,
                maxMajorSpheri=self.maxMajorSpheri,
                minMinorSpheri=self.minMinorSpheri,
                maxMinorSpheri=self.maxMinorSpheri,
                minLambdaSpheriMax=self.minLambdaSpheriMax,
                maxLambdaSpheriMax=self.maxLambdaSpheriMax,
                minLambdaSpheriMid=self.minLambdaSpheriMid,
                maxLambdaSpheriMid=self.maxLambdaSpheriMid,
                minLambdaSpheriMin=self.minLambdaSpheriMin,
                maxLambdaSpheriMin=self.maxLambdaSpheriMin,
                minPrincipalAxisTheta=self.minPrincipalAxisTheta,
                maxPrincipalAxisTheta=self.maxPrincipalAxisTheta,
            ) for seed in range(self.seedFrom, self.seedTo+1)]


class MiniJetFactorizationTask(luigi.Task):
    def run(self):
        # Make plot
        # MJFactorizatonPlots({
        #         MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', alignment='Spheri3DFull') : EColor.kBlack,
        #         MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, alignment='Spheri3DFull') : EColor.kRed,
        #         MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull') : EColor.kBlue,
        #     })
        MJFactorizatonPlots({
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, alignment='No') : EColor.kBlue,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, alignment='Spheri3DFull') : EColor.kCyan,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No') : EColor.kRed,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull') : EColor.kRed+2,
            })
        MJFactorizatonPlots({
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', alignment='No') : EColor.kBlue,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', alignment='Spheri3DFull') : EColor.kCyan,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No') : EColor.kRed,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull') : EColor.kRed+2,
            })
        MJFactorizatonPlots({
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=4, maxMult=10, alignment='No') : EColor.kBlue,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=4, maxMult=10, alignment='Spheri3DFull') : EColor.kCyan,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=4, maxMult=10, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No') : EColor.kRed,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=4, maxMult=10, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull') : EColor.kRed+2,
            })
        
        
        
        
        MJFactorizatonPlots({
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, alignment='No') : EColor.kBlue,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, alignment='Spheri3DFull') : EColor.kCyan,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.7, maxMajorSpheri=1, alignment='No') : EColor.kRed,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.7, maxMajorSpheri=1, alignment='Spheri3DFull') : EColor.kRed+2,
            })
        
        MJFactorizatonPlots({
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', alignment='No') : EColor.kRed,
                MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', alignment='Spheri3DFull') : EColor.kRed+2,
            })
        # MJFactorizatonPlots({
        #         MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No') : EColor.kRed,
        #         MakeSuffix(pdg1=2212, pdg2=-2212, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull') : EColor.kBlue,
        #     })
        # MJFactorizatonPlots({
        #         MakeSuffix(pdg1=2212, pdg2=3122, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No') : EColor.kRed,
        #         MakeSuffix(pdg1=2212, pdg2=3122, computeKStarMode='Custom', minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull') : EColor.kBlue,
        #     })
        # MJFactorizatonPlots({
        #         MakeSuffix(pdg1=2212, pdg2=3122, computeKStarMode='Custom', minMult=25, maxMult=30, alignment='No') : EColor.kRed,
        #         MakeSuffix(pdg1=2212, pdg2=3122, computeKStarMode='Custom', minMult=25, maxMult=30, alignment='Spheri3DFull') : EColor.kBlue,
        #     })
        # MJFactorizatonPlots({
        #         MakeSuffix(pdg1=22, pdg2=22, computeKStarMode='Custom', alignment='No') : EColor.kBlack,
        #         MakeSuffix(pdg1=22, pdg2=22, minMult=5, maxMult=8, computeKStarMode='Custom', alignment='No') : EColor.kRed,
        #         MakeSuffix(pdg1=22, pdg2=22, minMult=5, maxMult=8, computeKStarMode='Custom', alignment='Spheri3DFull') : EColor.kBlue,
        #     })
        # MJFactorizatonPlots({
        #         MakeSuffix(pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, computeKStarMode='Custom', minMajorSpheri=0.8, maxMajorSpheri=0.9, alignment='No',) : EColor.kRed,
        #         MakeSuffix(pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, computeKStarMode='Custom', minMajorSpheri=0.8, maxMajorSpheri=0.9, alignment='Spheri3DFull',) : EColor.kBlue,
        #     })
        # MJFactorizatonPlots({
        #         MakeSuffix(pdg1=411, pdg2=-211, computeKStarMode='Custom', alignment='No',) : EColor.kBlack,
        #     })




    def requires(self):
        ncpu = 16
        return [
            SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, alignment='No', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=500000, pdg1=2212, pdg2=-2212, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=500000, pdg1=2212, pdg2=-2212, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            
            SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, minMult=4, maxMult=10, alignment='No', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, minMult=4, maxMult=10, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=500000, pdg1=2212, pdg2=-2212, minMult=4, maxMult=10, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=500000, pdg1=2212, pdg2=-2212, minMult=4, maxMult=10, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            
            SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, alignment='No', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=500000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=500000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=300000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, minMajorSpheri=0.7, maxMajorSpheri=1, alignment='No', seedFrom=1, seedTo=ncpu),
            SimulationTask(nEvents=300000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, minMajorSpheri=0.7, maxMajorSpheri=1, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),

            # SimulationTask(nEvents=10000, pdg1=2212, pdg2=-2212, computeKStarMode='FD', alignment='No', seedFrom=1, seedTo=32),
            
            # SimulationTask(nEvents=10000, pdg1=2212, pdg2=-2212, alignment='No', seedFrom=1, seedTo=32),

            # SimulationTask(nEvents=10000, pdg1=2212, pdg2=-2212, alignment='Spheri3DFull', seedFrom=1, seedTo=32),

            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=3122, alignment='No', seedFrom=1, seedTo=ncpu),
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=3122, minMult=25, maxMult=30, alignment='No', seedFrom=1, seedTo=ncpu*3),
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=3122, minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No', seedFrom=1, seedTo=ncpu*6),
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=3122, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=3122, minMult=25, maxMult=30, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu*3),
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=3122, minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu*6),

            
            # # spherical
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, minMajorSpheri=0.8, maxMajorSpheri=0.9, alignment='No', seedFrom=1, seedTo=ncpu*6),
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=-2212, minMult=25, maxMult=30, minMajorSpheri=0.8, maxMajorSpheri=0.9, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu*6),
            
            # SimulationTask(nEvents=100000, pdg1=2212, pdg2=2212, minMult=25, maxMult=30, alignment='No', seedFrom=1, seedTo=22),

            # SimulationTask(nEvents=20000, pdg1=22, pdg2=22, alignment='No', seedFrom=1, seedTo=8),
            # SimulationTask(nEvents=20000, pdg1=22, pdg2=22, minMult=5, maxMult=8, alignment='No', seedFrom=1, seedTo=32),
            # SimulationTask(nEvents=20000, pdg1=22, pdg2=22, minMult=5, maxMult=8, alignment='Spheri3DFull', seedFrom=1, seedTo=32),

            # # check kstar computation
            # SimulationTask(nEvents=1000, pdg1=211, pdg2=-211, computeKStarMode='FD', alignment='No', seedFrom=1, seedTo=1),
            # SimulationTask(nEvents=1000, pdg1=211, pdg2=-211, computeKStarMode='Custom', alignment='No', seedFrom=1, seedTo=1),

            # # Dpi
            # SimulationTask(nEvents=100000, pdg1=411, pdg2=-211, alignment='No', seedFrom=1, seedTo=ncpu),
            # SimulationTask(nEvents=100000, pdg1=411, pdg2=211, minMult=25, maxMult=30, alignment='No', seedFrom=1, seedTo=ncpu*3),
            # SimulationTask(nEvents=100000, pdg1=411, pdg2=211, minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='No', seedFrom=1, seedTo=ncpu*6),
            # SimulationTask(nEvents=100000, pdg1=411, pdg2=211, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu),
            # SimulationTask(nEvents=100000, pdg1=411, pdg2=211, minMult=25, maxMult=30, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu*3),
            # SimulationTask(nEvents=100000, pdg1=411, pdg2=211, minMult=25, maxMult=30, minMajorSpheri=0.1, maxMajorSpheri=0.3, alignment='Spheri3DFull', seedFrom=1, seedTo=ncpu*6),

        ]

if __name__ == '__main__':
    luigi.build([
        MiniJetFactorizationTask(),
    ], workers=16)