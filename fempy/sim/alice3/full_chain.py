#!/usr/bin/env python3
from acts.examples.reconstruction import (
    addSeeding,
    # TruthSeedRanges,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    # SeedingAlgorithm,
    # ParticleSmearingSigmas,
    addCKFTracks,
#    CKFPerformanceConfig,
TrackSelectorConfig,
    # addKalmanTracks,
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

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('oDir')
parser.add_argument('--evt', type=int, required=True)
parser.add_argument('-B', default=1, type=float)
parser.add_argument('--dipole', default=False, action='store_true')
parser.add_argument('--pileup', type=int, default=1)
parser.add_argument('--hf', action='store_true', default=False)
parser.add_argument('--forced-decays', action='store_true', default=False)
parser.add_argument('--geom', default='20230731_fixed_endcaps')
parser.add_argument('--seed', default=42, type=int)
args = parser.parse_args()

oDir = args.oDir
geo_dir = pathlib.Path(f'/home/ktas/ge86rim/phsw/fempy/fempy/sim/alice3/geom/{args.geom}')

# IA: logging details
myVerboseGen = acts.logging.INFO
myVerboseDigi = acts.logging.INFO
myVerboseFatras = acts.logging.INFO
myVerboseSeeding = acts.logging.INFO
myVerboseCKFTracks = acts.logging.INFO
myVerboseAmbiguityResolution = acts.logging.WARNING #INFO

IA_MF = args.B # magnetic field

IA_nMeasurementsMin = 5
IA_maximumSharedHits = 1#3 #3

IA_outputDirName = "output/"
IA_outputDirName += f"evt-{args.evt}"
IA_outputDirName += f"_B-{args.B:.1f}T"
IA_outputDirName += f"_dipole" if args.dipole else ''
IA_outputDirName += f"_geom-{args.geom[:8]}"
IA_outputDirName += f"_pileup-{args.pileup}"
IA_outputDirName += f"_hf" if args.hf else ''
IA_outputDirName += f"_forced-decays" if args.forced_decays else ''
IA_outputDirName += f"_seed-{args.seed}"

outputDir = pathlib.Path(oDir) / IA_outputDirName

if not outputDir.exists():
    outputDir.mkdir(mode = 0o777, parents= True, exist_ok= True)

detector, trackingGeometry, decorators = alice3.buildALICE3Geometry(geo_dir, True, False, acts.logging.INFO)
if args.dipole:
    field = acts.examples.MagneticFieldMapXyz('geom/fieldMap_solenoid_dipoles_4-10-21_xyz_in_mm.txt')
else:
    field = acts.ConstantBField(acts.Vector3(0.0, 0.0, IA_MF * u.T))
rnd = acts.examples.RandomNumbers(seed=args.seed)

s = acts.examples.Sequencer(events=args.evt, numThreads=1)
s.trackFpes=False


# Monash tune
settings = [
    'Tune:pp = 14',
    'SoftQCD:nonDiffractive = on'
]

# override the monash tune and activate only the processes that lead to charm production
if args.hf:
    settings = [
        "HardQCD:gg2ccbar = on",
        "HardQCD:qqbar2ccbar = on",
        "HardQCD:gg2bbbar = on",
        "HardQCD:qqbar2bbbar = on",
    ]

if args.forced_decays:
    settings += [
        "421:onMode = off",
        "411:onMode = off",
        "413:onMode = off",
        "423:onMode = off",
        "421:onIfMatch = 211 321",
        "411:onIfMatch = 211 211 321",
        "413:onIfMatch = 211 421",
        "423:onIfMatch = 22 421",
    ]

s = addPythia8(
    s,
    npileup=args.pileup,
    beam=acts.PdgParticle.eProton,
    cmsEnergy= 13 * acts.UnitConstants.TeV,
    hardProcess=settings,
    vtxGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 0.1 * u.mm, 5.0 * u.ns),
        mean=acts.Vector4(0, 0, 0, 0),
    ),
    rnd=rnd,
    outputDirRoot=outputDir,
    # printPythiaEventListing='long'
)

# interaction w detector material like geant but simplified 
s = addFatras(
    s,
    trackingGeometry,
    field,
    enableInteractions=True,
    rnd=rnd,
    preSelectParticles= ParticleSelectorConfig(
        eta=(-4.0, 4.0), pt=(10 * u.MeV, None), removeNeutral=False),

#        postSelectParticles=ParticleSelectorConfig(
#    #        eta=(0.0, 4.0), pt=(500 * u.MeV, None), removeNeutral=False),
#            eta=(0.0, 4.0), pt=(50 * u.MeV, None), removeNeutral=False),
    outputDirRoot=outputDir,
    logLevel = myVerboseFatras, #acts.logging.VERBOSE,
)

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
    SeedFinderConfigArg(
        r=(None, 200 * u.mm),
        deltaR=(1 * u.mm, 60 * u.mm),
        collisionRegion=(-100 * u.mm, 100 * u.mm),
        z=(-2000 * u.mm, 2000 * u.mm),
        maxSeedsPerSpM=1,
        sigmaScattering=3,#3, #1.,#5.,
        radLengthPerSeed=0.01,#0.008,  #0.005,  # default: float radLengthPerSeed = 0.05; (https://github.com/acts-project/acts/blob/main/Core/include/Acts/Seeding/SeedFinderConfig.hpp)
        minPt=100* u.MeV, #60 * u.MeV, # GOOD IA! was: 500 * u.MeV,
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
    # CKFPerformanceConfig(ptMin=0.06* u.GeV, nMeasurementsMin=IA_nMeasurementsMin),   # IA, was: 500.0 * u.MeV
    TrackSelectorConfig(pt=(0.05*u.GeV, 100*u.GeV), nMeasurementsMin=IA_nMeasurementsMin),   # IA, was: 500.0 * u.MeV
    outputDirRoot=outputDir,
    writeTrajectories=True,
    logLevel = myVerboseCKFTracks,
)

s = addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(maximumSharedHits=IA_maximumSharedHits, nMeasurementsMin=IA_nMeasurementsMin),
#    CKFPerformanceConfig(
#        ptMin=0.06* u.GeV, #1.0 * u.GeV if ttbar_pu200 else 0.0,   # IA, was: 0.5* u.GeV
#        nMeasurementsMin=IA_nMeasurementsMin,
#    ),
    outputDirRoot=outputDir,
    # outputDirCsv=outputDir,
    writeTrajectories=True,
    logLevel = myVerboseAmbiguityResolution,
)

# trackFinder = acts.examples.TrackFindingAlgorithm(
#       level=acts.logging.INFO,
#       measurementSelectorCfg=acts.MeasurementSelector.Config(
#           [(acts.GeometryIdentifier(), ([], [15.0], [10]))]
#       ),
#       inputMeasurements="measurements",
#       inputSourceLinks="sourcelinks",
#       inputInitialTrackParameters="estimatedparameters",
#       outputTracks="ckfTracks",
#     #   outputTracks="ckfTracks",
#       findTracks=acts.examples.TrackFindingAlgorithm.makeTrackFinderFunction(
#           trackingGeometry, field, acts.ActsPythonBindings.logging.Level.INFO
#       ),
#   )
# s.addAlgorithm(trackFinder)

# trackConverter = acts.examples.TracksToTrajectories(
#    level=acts.logging.INFO,
#    inputTracks=trackFinder.config.outputTracks,
#    outputTrajectories="trajectories-from-tracks",
# )
# s.addAlgorithm(trackConverter)

# s = addKalmanTracks(
#    s,
#    trackingGeometry,
#    field,
#    False, #    directNavigation,
#    0 * u.GeV, #    reverseFilteringMomThreshold,
#    inputProtoTracks = "trajectories-from-tracks" # "truth_particle_tracks" # "seed-prototracks",
#    #outputDirRoot=outputDir
# #    writeTrajectories=True,
# )

s.run()

# import shutil
# shutil.copyfile( 'full_chain_acts_27.py', IA_outputDirName+'/full_chain_IA_check_steps.py' )
# shutil.copyfile( 'bashRunChainWithVariousParams.sh', IA_outputDirName+'/bashRunChainWithVariousParams.sh' )
