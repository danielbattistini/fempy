#!/usr/bin/env python3
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    ParticleSmearingSigmas,
    addCKFTracks,
#    CKFPerformanceConfig,
TrackSelectorConfig,
    addKalmanTracks,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,

)
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
import pathlib
import acts
import acts.examples
import alice3

u = acts.UnitConstants

geo_dir = pathlib.Path.cwd()


# ### IA: logging details
myVerboseGen = acts.logging.INFO
myVerboseDigi = acts.logging.INFO
myVerboseFatras = acts.logging.INFO
myVerboseSeeding = acts.logging.INFO
myVerboseCKFTracks = acts.logging.INFO
myVerboseAmbiguityResolution = acts.logging.WARNING #INFO

import sys

heavyion = False

# ### IA: params by hand
IA_MF = 1.0 #float(sys.argv[1]) #T
IA_nEvents = 2000 #int(sys.argv[2]) #200
IA_mult = 7 #8000#7#12 #7#1 #7 #12 #5
#    IA_ptMin = 0.05 #0.7
#    IA_ptMax = 5.0
#    IA_ptMin = 1.0 #0.7
#    IA_ptMax = 1.01
IA_ptMin = 0.1#5#0.1 #0.1#0.1 #0.1#05 #5.0 #1.0
IA_ptMax = 100#5 #100 #2#100#100#2#2 #5 #0.1#0501 #5.001 #001
IA_etaMin = -0.08#5 #2.5005 #1.0
IA_etaMax = 0.08#5 #2.500501 #001

IA_nPileupPYTHIA = 5 #int(sys.argv[3]) #5



IA_nMeasurementsMin = 7#9#11#9#10#11 #7#5 #4#7 #7
IA_maximumSharedHits = 1#3 #3



#    IA_outputDirName = "output/IA_ImprovedHistBinning_changed_initialVarInflation_ckf_"+str(IA_nEvents)+"Ev_B"+str(IA_MF) \
#IA_outputDirName = "output/IA_ImprovedHistBinning_2xWorseResolutionOuter_ckf_"+str(IA_nEvents)+"Ev_B"+str(IA_MF) \

#IA_outputDirName = "output/IA_ToyFemto_protons_Barrel_Gauss2kStar_0.025_ckf_"+str(IA_nEvents)+"Ev_B"+str(IA_MF) \
IA_outputDirName = "output/pions_ckf_nMeas"+str(IA_nMeasurementsMin)+"_nSharedHits"+str(IA_maximumSharedHits)+"_"+str(IA_nEvents)+"Ev_B"+str(IA_MF) \
    +"_gun"+"_mult"+str(IA_mult)+"_eta"+str(IA_etaMin)+"-"+str(IA_etaMax) \
    +"_pt"+str(IA_ptMin)+"-"+str(IA_ptMax)+"GeV"


outputDir = pathlib.Path.cwd() / IA_outputDirName


if not outputDir.exists():
    outputDir.mkdir(mode = 0o777, parents= True, exist_ok= True)

detector, trackingGeometry, decorators = alice3.buildALICE3Geometry(
#    geo_dir, False, False, acts.logging.VERBOSE)
    geo_dir, True, False, acts.logging.INFO)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, IA_MF * u.T))
rnd = acts.examples.RandomNumbers(seed=42)


s = acts.examples.Sequencer(events=IA_nEvents, numThreads=-1)

if not heavyion:
    addParticleGun(
        s,
        MomentumConfig( IA_ptMin * u.GeV, IA_ptMax * u.GeV, transverse=True ),
#        MomentumConfig(1.0 * u.GeV, 1.001 * u.GeV, transverse=True),
#        MomentumConfig(0.5 * u.GeV, 0.501 * u.GeV, transverse=True),
#        MomentumConfig(0.2 * u.GeV, 0.201 * u.GeV, transverse=True),
#        MomentumConfig(0.5 * u.GeV, 2.0 * u.GeV, transverse=True),
#         EtaConfig(-1.0, 1.0, uniform=True),
         EtaConfig(IA_etaMin, IA_etaMax, uniform=True),
#        EtaConfig(0., 1., uniform=True),
#        EtaConfig(2, 4, uniform=True),
#        EtaConfig(0.02, 0.1, uniform=True),
#        EtaConfig(3.2, 3.21, uniform=True),
#        EtaConfig(1.4, 3.8, uniform=True),
#        EtaConfig(0, 0.1, uniform=True),
#        PhiConfig(0, 0.1), #, uniform=True),

#        ParticleConfig(IA_mult, acts.PdgParticle.eMuon, randomizeCharge=True),
#        ParticleConfig(IA_mult, acts.PdgParticle.ePionPlus, randomizeCharge=True),
#        ParticleConfig(IA_mult, acts.PdgParticle.ePionPlus, randomizeCharge=False), #True), # !!!!!!!!!
#        ParticleConfig(IA_mult, acts.PdgParticle.eMuon, randomizeCharge=True),
        ParticleConfig(IA_mult, acts.PdgParticle.ePionPlus, randomizeCharge=True),
#        ParticleConfig(IA_mult, acts.PdgParticle.ePionPlus, randomizeCharge=True),
#        ParticleConfig(IA_mult, acts.PdgParticle.eProton, randomizeCharge=False), #True), # !!!!!!!!!
        rnd=rnd,
        logLevel = myVerboseGen, #acts.logging.VERBOSE,
        outputDirRoot=outputDir,
    )

else:
    s = addPythia8(
        s,
        npileup=IA_nPileupPYTHIA,
        beam=acts.PdgParticle.eProton,  #eLead,
        cmsEnergy=5 * acts.UnitConstants.TeV,
        # hardProcess=["Top:qqbar2ttbar=on"],
        # npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
#            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 0.1 * u.mm, 5.0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
        outputDirRoot=outputDir,
    )




s = addFatras(
    s,
    trackingGeometry,
    field,
    enableInteractions=True,
    rnd=rnd,
    preSelectParticles= ParticleSelectorConfig(
#        eta=(0.0, 4.0), pt=(500 * u.MeV, None), removeNeutral=False),
        eta=(-4.0, 4.0), pt=(10 * u.MeV, None), removeNeutral=False),

#        postSelectParticles=ParticleSelectorConfig(
#    #        eta=(0.0, 4.0), pt=(500 * u.MeV, None), removeNeutral=False),
#            eta=(0.0, 4.0), pt=(50 * u.MeV, None), removeNeutral=False),


    outputDirRoot=outputDir,
    logLevel = myVerboseFatras, #acts.logging.VERBOSE,
)

#s.run()

#s=addGeant4(
#    s,
#    detector,
#    trackingGeometry,
#    field,
#    preSelectParticles=ParticleSelectorConfig(
#        eta=(-3.0, 3.0),
#        absZ=(0, 1e4),
#        rho=(0, 1e3),
#        pt=(50 * u.MeV, None),
#        removeNeutral=True,
#    ),
#    outputDirRoot=outputDir,
#    # outputDirCsv=outputDir,
#    rnd=rnd,
#)

#s.run()

s = addDigitization(
    s,
    trackingGeometry,
    field,
#    digiConfigFile=geo_dir / "alice3-smearing-config_InnerBarrelRes_5mkm.json", #"alice3-smearing-config_DEFAULT.json",
    digiConfigFile=geo_dir / "alice3-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
    logLevel = myVerboseDigi, #acts.logging.VERBOSE,
)







s = addSeeding(
    s,
    trackingGeometry,
    field,
#    TruthSeedRanges(pt=(0.5 * u.GeV, None), eta=(0, 4.0), nHits=(7, None)),
#    TruthSeedRanges(pt=(0.045 * u.GeV, None), eta=(-4.0, 4.0), nHits=(4, None)),
    SeedFinderConfigArg(
        r=(None, 200 * u.mm),
        deltaR=(1 * u.mm, 60 * u.mm),  # deltaR=(1. * u.mm, 60 * u.mm),
#        collisionRegion=(-250 * u.mm, 250 * u.mm),
        collisionRegion=(-100 * u.mm, 100 * u.mm),
        z=(-2000 * u.mm, 2000 * u.mm),
        maxSeedsPerSpM=1,
        sigmaScattering=3,#3, #1.,#5.,
        radLengthPerSeed=0.01,#0.008,  #0.005,  # default: float radLengthPerSeed = 0.05; (https://github.com/acts-project/acts/blob/main/Core/include/Acts/Seeding/SeedFinderConfig.hpp)
        minPt=50* u.MeV, #60 * u.MeV, # GOOD IA! was: 500 * u.MeV,
        impactMax=1. * u.mm, #1. * u.mm,
        cotThetaMax=27.2899,
        seedConfirmation=True,
        centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-620 * u.mm,
            zMaxSeedConf=620 * u.mm,
            rMaxSeedConf=4.9 * u.mm, #36 * u.mm,  # IA: dramatically affects acceptance at eta ~4! <5 * u.mm  gives best results
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
        forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
            zMinSeedConf=-1220 * u.mm,
            zMaxSeedConf=1220 * u.mm,
            rMaxSeedConf=15 * u.mm,  #36 * u.mm,
            nTopForLargeR=1,
            nTopForSmallR=2,
        ),
#        skipPreviousTopSP=True,
        useVariableMiddleSPRange=True,
        # deltaRMiddleMinSPRange=10 * u.mm,
        # deltaRMiddleMaxSPRange=10 * u.mm,
        deltaRMiddleSPRange=(1 * u.mm, 10 * u.mm),
    ),
    SeedFinderOptionsArg(bFieldInZ= IA_MF * u.T, beamPos=(0 * u.mm, 0 * u.mm)),
    SeedFilterConfigArg(
        seedConfirmation=True,
        maxSeedsPerSpMConf=1,
        maxQualitySeedsPerSpMConf=1,
    ),
    SpacePointGridConfigArg(
        # zBinEdges=[
        # -4000.0,
        # -2500.0,
        # -2000.0,
        # -1320.0,
        # -625.0,
        # -350.0,
        # -250.0,
        # 250.0,
        # 350.0,
        # 625.0,
        # 1320.0,
        # 2000.0,
        # 2500.0,
        # 4000.0,
        # ],
        impactMax=1. * u.mm,
        phiBinDeflectionCoverage=3,
    ),
    SeedingAlgorithmConfigArg(
        # zBinNeighborsTop=[
        # [0, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 0],
        # ],
        # zBinNeighborsBottom=[
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 1],
        # [0, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # [-1, 0],
        # ],
        # numPhiNeighbors=1,
    ),
    geoSelectionConfigFile=geo_dir /
    "geoSelection-alice3-cfg10.json",
#    seedingAlgorithm=SeedingAlgorithm.TruthSmeared, # ADDED BY IA
    outputDirRoot=outputDir,
#    initialVarInflation = (50,50,50,50,50,50)  # ADDED BY IA
#    initialVarInflation = (0.2,0.2,0.2,0.2,0.2,0.2)  # ADDED BY IA
)






s = addCKFTracks(
    s,
    trackingGeometry,
    field,
#    CKFPerformanceConfig(ptMin=0.06* u.GeV, nMeasurementsMin=IA_nMeasurementsMin),   # IA, was: 500.0 * u.MeV
    TrackSelectorConfig(pt=(0.05*u.GeV, 100*u.GeV), nMeasurementsMin=IA_nMeasurementsMin),   # IA, was: 500.0 * u.MeV
    outputDirRoot=outputDir,
    writeTrajectories=False, #True,
    logLevel = myVerboseCKFTracks,
)


#s.run()


s = addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(maximumSharedHits=IA_maximumSharedHits, nMeasurementsMin=IA_nMeasurementsMin),
#    CKFPerformanceConfig(
#        ptMin=0.06* u.GeV, #1.0 * u.GeV if ttbar_pu200 else 0.0,   # IA, was: 0.5* u.GeV
#        nMeasurementsMin=IA_nMeasurementsMin,
#    ),
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
#    writeTrajectories=False,
    logLevel = myVerboseAmbiguityResolution,
)



#trackFinder = acts.examples.TrackFindingAlgorithm(
#       level=acts.logging.INFO,
#       measurementSelectorCfg=acts.MeasurementSelector.Config(
#           [(acts.GeometryIdentifier(), ([], [15.0], [10]))]
#       ),
#       inputMeasurements="measurements",
#       inputSourceLinks="sourcelinks",
#       inputInitialTrackParameters="estimatedparameters",
#       outputTracks="ckfTracks",
#       findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
#           trackingGeometry, field
#       ),
#   )
#s.addAlgorithm(trackFinder)

#trackConverter = acts.examples.TracksToTrajectories(
#    level=acts.logging.INFO,
#    inputTracks=trackFinder.config.outputTracks,
#    outputTrajectories="trajectories-from-tracks",
#)
#s.addAlgorithm(trackConverter)

#s = addKalmanTracks(
#    s,
#    trackingGeometry,
#    field,
#    False, #    directNavigation,
#    0 * u.GeV, #    reverseFilteringMomThreshold,
#    inputProtoTracks = "trajectories-from-tracks" # "truth_particle_tracks" # "seed-prototracks",
#    #outputDirRoot=outputDir
##    writeTrajectories=True,
#)





s.run()



import shutil
shutil.copyfile( 'full_chain_acts_27.py', IA_outputDirName+'/full_chain_IA_check_steps.py' )
shutil.copyfile( 'bashRunChainWithVariousParams.sh', IA_outputDirName+'/bashRunChainWithVariousParams.sh' )





