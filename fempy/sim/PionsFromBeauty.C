#include <map>
#include <string>
#include <vector>

#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TSystem.h"
#include "functions.hxx"
#include "AliPythia8.h"


void PionsFromBeauty(int nEvents = 5000, int seed = 1) {
    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    AliPythia8 pythia8;

    // select physics processes
    pythia8.ReadString("Tune:pp = 14");
    pythia8.ReadString("HardQCD:hardbbbar = on");
    // pythia8.ReadString("SoftQCD:all = on");

    // set seed for simulation
    pythia8.ReadString(Form("Random:seed %d", seed));
    pythia8.ReadString("Random:setSeed = on");
    pythia8.Initialize(2212, 2212, 13000);

    TH1D *hPions = new TH1D("hPions", "", 100, -0.5, 99.5);
    TH1D *hPionsFromBeauty = new TH1D("hPionsFromBeauty", "", 100, -0.5, 99.5);

    const std::string oFileName = Form("CountPions_%d.root", seed);
    TFile oFile(oFileName.data(), "recreate");

    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia8.GenerateEvent();
        pythia8.ImportParticles(particles, "All");
        pythia8.EventListing();

        TParticle* partl = dynamic_cast<TParticle*>(particles->At(0));

        printf("zero :---> %d\n", partl->GetPdgCode());

        if (Trigger(particles, "411_inAcc")) {
            // Count total number of pions
            int nPions = 0;
            for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
                TParticle* part = dynamic_cast<TParticle*>(particles->At(iPart));
                if (std::abs(part->GetPdgCode()) == 211 && IsDetectableInTPC(part)) nPions++;
            }
            hPions->Fill(nPions);

            // Count number of pions from beauty quark
            int countBeauty = 0;
            for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
                std::cout <<"\n <<<<< " << iPart<< std::endl;
                
                TParticle* part = dynamic_cast<TParticle*>(particles->At(iPart));

                if (std::abs(part->GetPdgCode()) != 5 || !IsLast(particles, iPart)) continue;
                countBeauty++;
                printf("iB = %d   pdg = %d\n", iPart,part->GetPdgCode()  );
                auto parts = SelectInDaughterTree(particles, iPart, {211}, true);
                hPionsFromBeauty->Fill(parts.size());
            }
            printf("cb: %d\n", countBeauty);
        }
    }

    // save root output file
    hPions->Write();
    hPionsFromBeauty->Write();
    oFile.Close();

    printf("Output saved in %s\n", oFileName.data());
}