#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TTree.h>

#include <array>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include "functions.hxx"

//__________________________________________________________________________________________________
void PionPShape(std::string inFileName, bool selDauKinem, std::string oDir);
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2);


//__________________________________________________________________________________________________
void PionPShape(std::string inFileName, std::string oDir) {
    //__________________________________________________________
    // define outputs
    TH1D *hPtAllEvts = new TH1D("hPtAllEvts", "hPtAllEvts;pt;counts", 200, 0, 10);
    TH1D *hPtDmesonEvts = new TH1D("hPtDmesonEvts", "hPtDmesonEvts;pt;counts", 200, 0, 10);

    TClonesArray *particles = new TClonesArray("TParticle", 1000);

    TFile *inFile = TFile::Open(inFileName.data());
    if (!inFile) {
        printf("The file %s does not exist. Skip!", inFileName.data());
        return;
    }

    TTree *tEvents = (TTree *)inFile->Get("tEvents");
    tEvents->SetBranchStatus("particles", 1);
    tEvents->SetBranchAddress("particles", &particles);

    for (int iEvent = 0; iEvent < tEvents->GetEntries(); iEvent++) {
        tEvents->GetEntry(iEvent);

        bool thereIsD = ThereIsD(particles);

        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
            if (std::abs(particle->Y()) > 2)  // keep only midrapidity
                continue;

            double prodR = std::sqrt(particle->Vx() * particle->Vx() + particle->Vy() * particle->Vy() +
                                     particle->Vz() * particle->Vz());

            if (prodR > 1.)  // to remove products of strange decays
                continue;

            int pdg = particle->GetPdgCode();
            int absPdg = std::abs(pdg);

            if (absPdg == 211) {
                hPtAllEvts->Fill(particle->Pt());

                if (thereIsD) hPtDmesonEvts->Fill(particle->Pt());
            }
        }
    }

    // save root output file
    TFile outFile(Form("%s/AnalysisResults.root", oDir.data()), "recreate");
    hPtAllEvts->Scale(1./hPtAllEvts->GetEntries());
    hPtDmesonEvts->Scale(1./hPtDmesonEvts->GetEntries());
    hPtAllEvts->Write();
    hPtDmesonEvts->Write();

    outFile.Close();
}
