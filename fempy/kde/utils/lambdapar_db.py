'''Database with Lambda parameters'''

import sys

def sum_lambda_pars(params):
    '''Computes the sum of all the lambda parameters for a given analysis.'''

    return sum([params[h][l][0] for h in params for l in params[h]])


def get_lambda_par(analysis=None):
    '''Obtain lambda parameter for a certain analysis.'''

    if analysis == 'Dpi':
        lambda_pars = {}

        # origins of D: true, non-prompt, D*, comb background
        # origins of pi: true, decays and material, strong decays, misidentified

        lambda_pars['true'] = {}
        lambda_pars['true']['true'] = (0.403, 'gen')
        lambda_pars['true']['secondaries'] = (0.0023, 'flat')
        lambda_pars['true']['strong'] = (0.0550, 'flat')
        lambda_pars['true']['misidentified'] = (0.0047, 'flat')
        lambda_pars['nonprompt'] = {}
        lambda_pars['nonprompt']['true'] = (0.0443, 'flat')
        lambda_pars['nonprompt']['secondaries'] = (0.0003, 'flat')
        lambda_pars['nonprompt']['strong'] = (0.0060, 'flat')
        lambda_pars['nonprompt']['misidentified'] = (0.0005, 'flat')
        lambda_pars['Dstar'] = {}
        lambda_pars['Dstar']['true'] = (0.1680, 'coulomb')
        lambda_pars['Dstar']['secondaries'] = (0.0010, 'flat')
        lambda_pars['Dstar']['strong'] = (0.0229, 'flat')
        lambda_pars['Dstar']['misidentified'] = (0.0019, 'flat')
        lambda_pars['combbkg'] = {}
        lambda_pars['combbkg']['true'] = (0.2514, 'sideband')
        lambda_pars['combbkg']['secondaries'] = (0.0014, 'sideband')
        lambda_pars['combbkg']['strong'] = (0.0343, 'sideband')
        lambda_pars['combbkg']['misidentified'] = (0.0029, 'sideband')

        return lambda_pars, sum_lambda_pars(lambda_pars)

    print('Error: analysis not implemented. Exit!')
    sys.exit()


lp, s = get_lambda_par('Dpi')
print(s)
