'''
Script to control the submission of jobs in the batch farm.

Run with:
python3 JobManager.py [-sleepTime time] [--alwaysFree ncpu]
'''

import os
import subprocess
import time
import argparse

from tabulate import tabulate

parser = argparse.ArgumentParser()
parser.add_argument('--sleep', type=int, default=300, metavar='sleep', help='Sleep time between each update (seconds)')
parser.add_argument('--free', type=int, default=30, metavar='ncpu', help='Number of CPUs to be left always free')
args = parser.parse_args()

sleepTime = args.sleep
nAlwaysFree = args.free

while True:
    # My jobs
    cmd = 'squeue -h --me -t R  | wc -l'
    myRunningJobs = int(subprocess.check_output(cmd, shell=True, universal_newlines='\n'))
    
    cmd = 'squeue -h --me -t PD | grep -v JobHeldUser | wc -l'
    myPendingJobs = int(subprocess.check_output(cmd, shell=True, universal_newlines='\n'))
    
    cmd = 'squeue -h --me -t PD | grep JobHeldUser | wc -l'
    myHeldJobs = int(subprocess.check_output(cmd, shell=True, universal_newlines='\n'))

    # Their jobs
    cmd = 'squeue -h -t R  | grep -v $(whoami) | wc -l'
    theirRunningJobs = int(subprocess.check_output(cmd, shell=True, universal_newlines='\n'))
    
    cmd = 'squeue -h -t PD | grep -v $(whoami) | grep -v JobHeldUser | wc -l'
    theirPendingJobs = int(subprocess.check_output(cmd, shell=True, universal_newlines='\n'))
    
    cmd = 'squeue -h -t PD | grep -v $(whoami) | grep JobHeldUser | wc -l'
    theirHeldJobs = int(subprocess.check_output(cmd, shell=True, universal_newlines='\n'))
    
    theirTotalJobs = theirRunningJobs + theirPendingJobs + theirHeldJobs

    cmd = 'sinfo -aeN -o "%N %c %P" -h | grep -v test | grep -v slimfast | awk \'{print $1 " " $2}\' | uniq | ' \
          'awk \'{s+=$2} END {print s}\''
    nCpu = int(subprocess.check_output(cmd,  shell=True, universal_newlines='\n'))
    nCpu = 324 #! The ones with at least 2 GB of memory
    
    nRunningJobs = int(subprocess.check_output('squeue -h -t R | wc -l',  shell=True, universal_newlines='\n'))
    nFreeCPU = nCpu - nRunningJobs
    
    print(f'{"Running jobs: ":<20}{nRunningJobs:>7}')
    print(f'{"Total n. of CPUS: ":<20}{nCpu:>7}')
    print(f'{"Available CPUs: ":<20}{nFreeCPU:>7}')
    print(tabulate([['Your jobs', myRunningJobs, myPendingJobs, myHeldJobs],
                    ['Their jobs', theirRunningJobs, theirPendingJobs, theirHeldJobs]],
                    headers=['', 'running', 'pending', 'held']))
    if myPendingJobs == 0 and nFreeCPU > nAlwaysFree - theirRunningJobs:
        nSubmit = nCpu - nAlwaysFree - nRunningJobs
        cmd = f'squeue -h --me -t PD | grep JobHeldUser | awk \'{{print $1}}\' | head -n {nSubmit} | tr \'\n\' \',\''
        myHeldJobIDs = str(subprocess.check_output(cmd, shell=True, universal_newlines='\n'))[:-1]

        print('Releasing the following jobs jobs: ', myHeldJobIDs)
        os.system(f'scontrol release {myHeldJobIDs}')

    print('\n\n')
    time.sleep(sleepTime)
