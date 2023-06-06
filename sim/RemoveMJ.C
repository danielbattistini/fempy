#include <TClonesArray.h>
#include <TFile.h>
#include <TStopwatch.h>

#include <string>
#include <vector>

#include "AliPythia8.h"
#include "functions.hxx"

using namespace Pythia8;

namespace {
enum processes { kSoftQCD = 0, kHardQCD };
enum triggers { kMB = 0, kHighMultV0 };
}  // namespace


void RemoveMJ(int nEvents=1000, triggers trigger=kMB, tunes tune=kCRMode2, processes process=kSoftQCD,
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

    SetTune(pythia, tune);
    SetProcess(pythia, process);

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
            if (std::abs(particle->GetPdgCode()) == hPdg) {
                goto print;
            }
        }
    }

    print:
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (std::abs(particle->Eta()) < 0.8 && IsDetectable(std::abs(particle->GetPdgCode())))
            printf("%f, ", particle->Px());
    }
    printf("\n");
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (std::abs(particle->Eta()) < 0.8 && IsDetectable(std::abs(particle->GetPdgCode())))
            printf("%f, ", particle->Py());
    }
    printf("\n");
    
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (std::abs(particle->Eta()) < 0.8 && IsDetectable(std::abs(particle->GetPdgCode())))
            printf("%f, ", particle->Pz());
    }

}