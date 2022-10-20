#include <string>
#include <vector>

#include <TClonesArray.h>
#include <TFile.h>
#include <TStopwatch.h>

#include "AliPythia8.h"

using namespace Pythia8;

namespace {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
enum processes { kSoftQCD = 0, kHardQCD };
enum triggers { kMB = 0, kHighMultV0, kHF};

}  // namespace

bool IsDetectable(int absPdg) {
    if (absPdg == 11 ||   // electrons
        absPdg == 13 ||   // muons
        absPdg == 211 ||  // pions
        absPdg == 321 ||  // kaons
        absPdg == 2212    // protons
    )
        return true;

    return false;
}

void CreateHMDataset(int nEvents, double maxRunTime, double energy, triggers trigger, tunes tune, processes process, int seed, std::string outDir) {
    TStopwatch timer;
    timer.Start();

    //__________________________________________________________
    // create and configure pythia generator
    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

    AliPythia8 pythia;
    if (process == kSoftQCD) {
        pythia.ReadString("SoftQCD:all = on");
    } else if (process == kHardQCD) {
        pythia.ReadString("HardQCD:hardccbar = on");
        pythia.ReadString("HardQCD:hardbbbar = on");
    }

    // set tune
    if (tune == kMonash) {
        pythia.ReadString(Form("Tune:pp = 14"));
    } else if (tune == kCRMode0) {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 2.9");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.43");
        pythia.ReadString("ColourReconnection:timeDilationMode = 0");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.12");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode2) {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 0.3");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.20");
        pythia.ReadString("ColourReconnection:timeDilationMode = 2");
        pythia.ReadString("ColourReconnection:timeDilationPar = 0.18");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.15");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    } else if (tune == kCRMode3) {
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("ColourReconnection:mode = 1");
        pythia.ReadString("ColourReconnection:allowDoubleJunRem = off");
        pythia.ReadString("ColourReconnection:m0 = 0.3");
        pythia.ReadString("ColourReconnection:allowJunctions = on");
        pythia.ReadString("ColourReconnection:junctionCorrection = 1.15");
        pythia.ReadString("ColourReconnection:timeDilationMode = 3");
        pythia.ReadString("ColourReconnection:timeDilationPar = 0.073");
        pythia.ReadString("StringPT:sigma = 0.335");
        pythia.ReadString("StringZ:aLund = 0.36");
        pythia.ReadString("StringZ:bLund = 0.56");
        pythia.ReadString("StringFlav:probQQtoQ = 0.078");
        pythia.ReadString("StringFlav:ProbStoUD = 0.2");
        pythia.ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        pythia.ReadString("MultiPartonInteractions:pT0Ref = 2.05");
        pythia.ReadString("BeamRemnants:remnantMode = 1");
        pythia.ReadString("BeamRemnants:saturation = 5");
    }

    // init
    pythia.ReadString("Random:setSeed = on");
    pythia.ReadString(Form("Random:seed %d", seed));
    pythia.Initialize(2212, 2212, energy);

    std::string oFileName = "pythia_HMpp";
    if (energy == 13000) {
        oFileName += "13TeV";
    } else {
        printf("Error: energy not implemented. Exit!");
        return;
    }

    if (tune == kMonash) {
        oFileName += "_tune-Monash";
    } else if (tune == kCRMode0) {
        oFileName += "_tune-CRMode0";
    } else if (tune == kCRMode2) {
        oFileName += "_tune-CRMode2";
    } else if (tune == kCRMode3) {
        oFileName += "_tune-CRMode3";
    } else {
        printf("Error: tune not implemented. Exit!");
        return;
    }

    if (process == kSoftQCD) {
        oFileName += "_process-SoftQCD";
    } else if (process == kHardQCD) {
        oFileName += "_process-HardQCD";
    } else {
        printf("Error: process not implemented. Exit!");
        return;
    }
    oFileName += ".root";

    //__________________________________________________________
    // define outputs
    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    TTree* tEvents = new TTree("tEvents", "tEvents");
    tEvents->Branch("particles", "TClonesArray", &particles);

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        if (timer.RealTime() > maxRunTime * 3600) {
            int time = int(timer.RealTime());  // in seconds
            int nHours = int(time / 3600);
            int nMins = int((time % 3600) / 60);

            printf("Reached max run time: %d:%d\n", nHours, nMins);
            break;
        } else {
            timer.Continue();
        }

        pythia.GenerateEvent();
        pythia.ImportParticles(particles, "All");


        if (trigger == kMB) {
            int nCh = 0;
            for (auto iPart = 2; iPart < particles->GetEntries(); ++iPart) {
                TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
                int status = std::abs(particle->GetStatusCode());

                if (status == 1) { // count particles in the final state
                    nCh++;
                }
            }
            if (nCh > 2) { // reject diffractive events
                tEvents->Fill();
            }
        } else if (trigger == kHighMultV0) {
            // evaluate multiplicity at forward rapidity
            int nChForward = 0;
            for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
                TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
                int pdg = std::abs(particle->GetPdgCode());
                int status = std::abs(particle->GetStatusCode());
                float eta = particle->Eta();

                if (IsDetectable(pdg) && ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1)) &&
                    status == 1) {  // V0A and V0C acceptance
                    nChForward++;
                }
            }
            if (nChForward > 130) {
                tEvents->Fill();
            }
        } else if (trigger == kHF) { // look for c or b quarks
            for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
                TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
                int pdg = std::abs(particle->GetPdgCode());

                if (std::abs(pdg) == 4 || std::abs(pdg) == 5) {
                    tEvents->Fill();
                    break;
                }
            }
        }
    }

    // save root output file
    TFile outFile(Form("%s/%s", outDir.data(), oFileName.data()), "recreate");
    tEvents->Write();
    outFile.Close();
}