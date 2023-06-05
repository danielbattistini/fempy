def GetFemtoDreamPairId(chargeComb):
    '''
    Get ParticleX_ParticleY
    particleX is the heaviest of the two

    id: sc, oc, all

    return string
    '''

    if not isinstance(chargeComb, list):
        chargeComb = [chargeComb]
    
    pairStr ={}
    for comb in chargeComb:
        if comb=='sc':
            pairStr['pp'] = 'Particle0_Particle2'
            pairStr['mm'] = 'Particle1_Particle3'
        elif comb=='oc':
            pairStr['pm'] = 'Particle1_Particle2'
            pairStr['mp'] = 'Particle0_Particle3'
        elif comb=='all':
            pairStr['pp'] = 'Particle0_Particle2'
            pairStr['mm'] = 'Particle1_Particle3'
            pairStr['pm'] = 'Particle1_Particle2'
            pairStr['mp'] = 'Particle0_Particle3'
        else:
            print('\033[31mError: particle pair not implemented\0330m')
        

    return pairStr