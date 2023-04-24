import os
import luigi


class ComputeDistrTask(luigi.Task):
    nEvents = luigi.IntParameter(default=10000)
    
    def run(self):
        # compile the macro
        os.system('root -l -b -q -e ".L MiniJetFactorization.cc+"')

        command = 'root -l -b -q -e ".L MiniJetFactorization_cc.so && MiniJetFactorization('\
            f'{self.nEvents}'\
            '")'
        

        print(command)

if __name__ == '__main__':
    luigi.build([
        ComputeDistrTask(),
    ], workers=1)