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

void BeautyBkg(int nEvents = 5000, int pdg1 = 411, int pdg2 = 211, int seed = 1, bool requireInAcc = false) {
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
    pythia8.ReadString(Form("Random:seed %d", seed));
    pythia8.ReadString("Random:setSeed = on");
    pythia8.ReadString("HardQCD:hardbbbar = on");
    pythia8.ReadString("411:onMode = off");
    pythia8.ReadString("411:onIfMatch = 211 211 321");
    pythia8.Initialize(2212, 2212, 13000);

    std::map<std::string, TH1*> hSEBeauty;
    hSEBeauty.insert({"p02_13", new TH1F("hSEBeauty_02_13", ";#it{k*} (GeV/c);Counts", 3000, 0., 3.)});
    hSEBeauty.insert({"p03_12", new TH1F("hSEBeauty_03_12", ";#it{k*} (GeV/c);Counts", 3000, 0., 3.)});
    hSEBeauty.insert({"counts_p02_13", new TH2F("hCountsBeauty_p02_13", ";Mult of D^{+};Mult of #pi^{+};Counts", 20,
                                                -0.5, 19.5, 20, -0.5, 19.5)});
    hSEBeauty.insert({"counts_p03_12", new TH2F("hCountsBeauty_p03_12", ";Mult of D^{+};Mult of #pi^{#minus};Counts",
                                                20, -0.5, 19.5, 20, -0.5, 19.5)});

    std::map<std::string, TH1*> hSETot;
    hSETot.insert({"p02_13", new TH1F("hSETot_02_13", ";#it{k*} (GeV/c);Counts", 3000, 0., 3.)});
    hSETot.insert({"p03_12", new TH1F("hSETot_03_12", ";#it{k*} (GeV/c);Counts", 3000, 0., 3.)});
    hSETot.insert({"counts_p02_13", new TH2F("hCountsTot_p02_13", ";Mult of D^{+};Mult of #pi^{+};Counts", 20, -0.5,
                                             19.5, 20, -0.5, 19.5)});
    hSETot.insert({"counts_p03_12", new TH2F("hCountsTot_p03_12", ";Mult of D^{+};Mult of #pi^{#minus};Counts", 20,
                                             -0.5, 19.5, 20, -0.5, 19.5)});

    const std::string oFileName = Form("Distr_BeautyBkg_%d_%d_reqInAcc-%s_seed%d.root", pdg1, pdg2, requireInAcc ? "Yes" : "No", seed);
    TFile oFile(oFileName.data(), "recreate");

    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        pythia8.GenerateEvent();
        pythia8.ImportParticles(particles, "All");

        // Correlated background: D and pi from same b quark
        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));

            if (std::abs(particle->GetPdgCode()) != 5 || !IsLast(particles, iPart)) continue;

            // Filter for particles in the TPC acceptance
            auto parts = SelectInDaughterTree(particles, iPart, {pdg1, pdg2}, requireInAcc);
            std::vector<int> daus = {};
            std::copy(parts.begin(), parts.end(), std::back_inserter(daus));

            // Compute correlated SE
            for (size_t iDau1 = 0; iDau1 < daus.size(); iDau1++) {
                TParticle* dau1 = dynamic_cast<TParticle*>(particles->At(daus[iDau1]));
                if (std::abs(dau1->GetPdgCode()) != pdg1) continue;

                for (size_t iDau2 = iDau1 + 1; iDau2 < daus.size(); iDau2++) {
                    TParticle* dau2 = dynamic_cast<TParticle*>(particles->At(daus[iDau2]));
                    if (std::abs(dau2->GetPdgCode()) != pdg2) continue;
                    if (ProductionRadius(dau2) > 1) continue;

                    std::string pair = dau1->GetPdgCode() * dau2->GetPdgCode() > 0 ? "p02_13" : "p03_12";
                    hSEBeauty[pair.data()]->Fill(RelativePairMomentum(dau1, dau2));
                }
            }

            // Count daus
            int count1Pos = 0;
            int count1Neg = 0;
            int count2Pos = 0;
            int count2Neg = 0;
            for (const int& iDau : daus) {
                TParticle* dau = dynamic_cast<TParticle*>(particles->At(iDau));

                if (dau->GetPdgCode() == +pdg1) count1Pos++;
                if (dau->GetPdgCode() == -pdg1) count1Neg++;
                if (dau->GetPdgCode() == +pdg2) count2Pos++;
                if (dau->GetPdgCode() == -pdg2) count2Neg++;
            }
            hSEBeauty["counts_p02_13"]->Fill(count1Pos, count2Pos);
            hSEBeauty["counts_p02_13"]->Fill(count1Neg, count2Neg);
            hSEBeauty["counts_p03_12"]->Fill(count1Pos, count2Neg);
            hSEBeauty["counts_p03_12"]->Fill(count1Neg, count2Pos);
        }

        // Uncorrelated background: all D and pions
        for (int iPart1 = 2; iPart1 < particles->GetEntriesFast(); iPart1++) {
            TParticle* part1 = dynamic_cast<TParticle*>(particles->At(iPart1));

            if (std::abs(part1->GetPdgCode()) != pdg1) continue;
            if (requireInAcc && !IsInAcc(particles, iPart1)) continue;

            for (int iPart2 = 2; iPart2 < particles->GetEntriesFast(); iPart2++) {
                TParticle* part2 = dynamic_cast<TParticle*>(particles->At(iPart2));
                if (std::abs(part2->GetPdgCode()) != pdg2 || !IsInAcc(particles, iPart2)) continue;
                if (ProductionRadius(part2) > 1) continue;

                if (part1->GetFirstDaughter() <= iPart2 && iPart2 <= part1->GetLastDaughter()) continue;

                std::string pair = part1->GetPdgCode() * part2->GetPdgCode() > 0 ? "p02_13" : "p03_12";
                hSETot[pair.data()]->Fill(RelativePairMomentum(part1, part2));
                hSETot[pair.data()]->Fill(RelativePairMomentum(part1, part2));
            }
        }

        // Count particles
        int count1Pos = 0;
        int count1Neg = 0;
        int count2Pos = 0;
        int count2Neg = 0;
        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle* part = dynamic_cast<TParticle*>(particles->At(iPart));

            if (IsInAcc(particles, iPart) && part->GetPdgCode() == +pdg1) count1Pos++;
            if (IsInAcc(particles, iPart) && part->GetPdgCode() == -pdg1) count1Neg++;
            if (IsInAcc(particles, iPart) && part->GetPdgCode() == +pdg2) count2Pos++;
            if (IsInAcc(particles, iPart) && part->GetPdgCode() == -pdg2) count2Neg++;
        }
        hSETot["counts_p02_13"]->Fill(count1Pos, count2Pos);
        hSETot["counts_p02_13"]->Fill(count1Neg, count2Neg);
        hSETot["counts_p03_12"]->Fill(count1Pos, count2Neg);
        hSETot["counts_p03_12"]->Fill(count1Neg, count2Pos);
    }

    // save root output file
    oFile.cd();
    for (const auto& pair : hSEBeauty) {
        pair.second->Write();
    }
    for (const auto& pair : hSETot) {
        pair.second->Write();
    }

    oFile.Close();

    printf("Output saved in %s\n", oFileName.data());
}