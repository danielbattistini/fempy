#define ALIPYTHIA 1

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

#if ALIPYTHIA
#include "AliPythia8.h"
#else
#include "TPythia8.h"
#endif



void TestAcceptanceFunc(int nEvents = 5000, int pdg = 211) {
#if ALIPYTHIA
    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    AliPythia8 pythia8;
#else
    TPythia8 pythia8;
#endif
    // select physics processes

    pythia8.ReadString("Tune:pp = 14");
    pythia8.ReadString("SoftQCD:all = on");

    // set seed for simulation
    pythia8.ReadString(Form("Random:seed 1"));
    pythia8.ReadString("Random:setSeed = on");
    pythia8.ReadString("HardQCD:hardccbar = on");
    
    pythia8.Initialize(2212, 2212, 13000);

    TH1F *hEta = new TH1F("hEta", ":Eta", 200, -3, 3.);
    TH1F *hEtaInAcc = new TH1F("hEtaInAcc", ":Eta", 200, -3, 3.);
    const std::string oFileName = Form("TestDistr_AccFunc_%d.root", pdg);
    TFile oFile(oFileName.data(), "recreate");

    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia8.GenerateEvent();
        pythia8.ImportParticles(particles, "All");

        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));

            if (std::abs(particle->GetPdgCode()) == pdg) {
                if (IsDetectable(pdg)) {
                    hEta->Fill(particle->Eta());
                    if (IsInAcc(particles, iPart)) hEtaInAcc->Fill(particle->Eta());
                }
                else {
                    for (int iDau = particle->GetFirstDaughter(); iDau <= particle->GetLastDaughter(); iDau++) {
                        TParticle* dau = dynamic_cast<TParticle*>(particles->At(iDau));

                        hEta->Fill(dau->Eta());
                        if (IsInAcc(particles, iPart)) hEtaInAcc->Fill(dau->Eta());
                    }
                }
            }
        }
    }

    hEta->Write();
    hEtaInAcc->Write();

    oFile.Close();

    printf("Output saved in TestDistr_AccFunc_%d.root\n", pdg);
}