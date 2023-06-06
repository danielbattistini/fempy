#include <TClonesArray.h>
#include <TFile.h>
#include <TStopwatch.h>

#include <string>
#include <vector>

#include "AliPythia8.h"

#include "functions.hxx"

using namespace Pythia8;

namespace {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
enum processes { kSoftQCD = 0, kHardQCD };
enum triggers { kMB = 0, kHighMultV0 };
}  // namespace

void MakeMomentumList(int nEvents=100000, triggers trigger=kMB, tunes tune=kCRMode2, processes process=kSoftQCD,
                     int hPdg=3122, int lPdg=-3122, int seed = 42, std::string outDir="") {
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
        pythia.ReadString("SoftQCD:nonDiffractive = on");
        // pythia.ReadString("SoftQCD:all = on");
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
    pythia.Initialize(2212, 2212, 13000);



    int fwMultThr = 0;
    if (trigger == kHighMultV0 && process == kSoftQCD && tune == kMonash)
        fwMultThr = 130;
    else if (trigger == kHighMultV0 && process == kHardQCD && tune == kMonash)
        fwMultThr = 133;
    else
        printf("error. mult thr not implemented");

    //__________________________________________________________
    // define outputs
    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    TTree* tEvents = new TTree("tEvents", "tEvents");
    tEvents->Branch("particles", "TClonesArray", &particles);

    int nCharm = 0;
    for (int iEvent = 0; iEvent < nEvents; iEvent++) {

        pythia.GenerateEvent();
        pythia.ImportParticles(particles, "All");

        for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
            if (GetMultV0(particles) > 130) {
                goto print;
            }
        }
    }

    print:
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (std::abs(particle->Eta()) < 2 && IsDetectable(std::abs(particle->GetPdgCode())))
            printf("%f, ", particle->Px());
    }
    printf("\n");
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (std::abs(particle->Eta()) < 2 && IsDetectable(std::abs(particle->GetPdgCode())))
            printf("%f, ", particle->Py());
    }
    printf("\n");
    
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (std::abs(particle->Eta()) < 2 && IsDetectable(std::abs(particle->GetPdgCode())))
            printf("%f, ", particle->Pz());
    }
    printf("\n");
    
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (std::abs(particle->Eta()) < 2 && IsDetectable(std::abs(particle->GetPdgCode())))
            printf("%d, ", 1 ? particle->GetPdgCode()>0 : -1);
    }
    printf("\n");

}