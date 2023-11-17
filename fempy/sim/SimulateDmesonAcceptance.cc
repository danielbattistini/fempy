#include <array>
#include <deque>
#include <map>
#include <string>
#include <vector>

#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3F.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>

#include "AliPythia8.h"
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "functions.hxx"

using namespace Pythia8;

void SimulateDmesonAcceptance(
    int nEvents = 1000000,
    int pdg = 411,
    int seed = 42,
    double minMult = 0,
    double maxMult = 1.e10,
    double minMajorSpheri = 0,
    double maxMajorSpheri = 1,
    double minMinorSpheri = 0,
    double maxMinorSpheri = 1,
    double minLambdaSpheriMax = 0,
    double maxLambdaSpheriMax = 1.e10,
    double minLambdaSpheriMid = 0,
    double maxLambdaSpheriMid = 1.e10,
    double minLambdaSpheriMin = 0,
    double maxLambdaSpheriMin = 1.e10,
    double minPrincipalAxisTheta = -TMath::Pi(),
    double maxPrincipalAxisTheta = +TMath::Pi(),
    std::string outFileNameRoot = "AnalysisResults.root") {

    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

    AliPythia8 pythia;
    // pythia.ReadString("SoftQCD:nonDiffractive = on");
    pythia.ReadString("HardQCD:hardccbar = on");
    pythia.ReadString("HardQCD:hardbbbar = on");

    pythia.ReadString(Form("Tune:pp = 14"));
    pythia.ReadString("Random:setSeed = on");
    pythia.ReadString(Form("Random:seed %d", seed));

    // force the char hadron decays
    pythia.ReadString("Random:setSeed = on");
    pythia.ReadString("411:onMode = off");
    pythia.ReadString("421:onMode = off");
    pythia.ReadString("413:onMode = off");
    pythia.ReadString("411:onIfMatch = 211 211 321");
    pythia.ReadString("421:onIfMatch = 211 321");
    pythia.ReadString("413:onIfMatch = 211 421");

    pythia.Initialize(2212, 2212, 13000);

    TFile oFile(outFileNameRoot.data(), "recreate");
    TH1F* hSel = new TH1F("hSel", ";#it{y};Counts", 300, -3, 3.);
    TH1F* hTot = new TH1F("hTot", ";#it{y};Counts", 300, -3, 3.);

    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia.GenerateEvent();
        pythia.ImportParticles(particles, "All");

        // Multiplicity selections
        int mult = GetMultTPC(particles);
        if (!(minMult <= mult && mult < maxMult)) continue;

        // Event-shape selections
        auto eivals = ComputeEigenValues(ComputeSphericityMatrix(particles));
        if (eivals[0] < eivals[1] && eivals[1] < eivals[2] && eivals[0] > 0) {  // sanity check on the eigenvalues
            double spheriMaj = ComputeMajorRelativeSphericity(particles);
            double spheriMin = ComputeMinorRelativeSphericity(particles);
            double PrincipalAxisTheta = ComputePrincipalAxisTheta(particles);
            if (!(minLambdaSpheriMax <= eivals[2] && eivals[2] < maxLambdaSpheriMax &&
                  minLambdaSpheriMid <= eivals[1] && eivals[1] < maxLambdaSpheriMid &&
                  minLambdaSpheriMin <= eivals[0] && eivals[0] < maxLambdaSpheriMin && minMajorSpheri <= spheriMaj &&
                  spheriMaj < maxMajorSpheri && minMinorSpheri <= spheriMin && spheriMin < maxMinorSpheri &&
                  minPrincipalAxisTheta <= PrincipalAxisTheta && PrincipalAxisTheta < maxPrincipalAxisTheta))
                continue;
        }  // todo : add option to reject events with bad event-shape

        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));

            if (pdg == particle->GetPdgCode()) {
                hTot->Fill(particle->Y());

                int firstDauIdx = particle->GetFirstDaughter();
                int lastDauIdx = particle->GetLastDaughter();

                bool allDauInAcc = true;
                for (int iDau = firstDauIdx; iDau <= lastDauIdx; iDau++) {
                    TParticle* dau = dynamic_cast<TParticle*>(particles->At(iDau));
                    if (!isInTPC(dau)) {
                        allDauInAcc = false;
                        break;
                    }
                }

                if (allDauInAcc) {
                    hSel->Fill(particle->Y());
                }
            } // todo: implement D*
        }
    }

    // save root output file
    hSel->Write();
    hTot->Write();
    oFile.Close();
}
