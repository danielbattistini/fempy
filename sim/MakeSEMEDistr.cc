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
#include "./functions.hxx"

std::array<int, 2> lightPDG{211, 321};  // pions, kaons, protons
std::array<int, 2> charmPDG{411, 413};  // D+, D0, Ds, D*, Lc

void MakeSEMEDistr(std::string inFileName, std::string oFileName, bool selDauKinem, bool align) {
    // Define outputs
    std::map<std::tuple<int, int, std::string, std::string>, TH2F*> hPairSE, hPairME;
    for (auto& pdgC : charmPDG) {
        for (auto& pdgL : lightPDG) {
            for (std::string origin : {"prompt", "nonprompt"}) {
                std::string name = Form("hPairSE_%d_%d_%s", pdgC, pdgL, origin.data());
                std::string title = ";#it{k}* (GeV/#it{c});#it{p}_{T}(H_{c}) (GeV/#it{c});pairs";
                hPairSE[{pdgC, pdgL, "part", origin}] = new TH2F(name, title, 2000, 0., 2., 10, 0, 10);

                name = Form("hPairSE_%d_%d_%s", pdgC, -pdgL, origin.data());
                title = ";#it{k}* (GeV/#it{c});#it{p}_{T}(H_{c}) (GeV/#it{c});pairs";
                hPairSE[{pdgC, pdgL, "antipart", origin}] = new TH2F(name, title, 2000, 0., 2., 10, 0, 10);

                name = Form("hPairME_%d_%d_%s", pdgC, pdgL, origin.data());
                title = ";#it{k}* (GeV/#it{c});#it{p}_{T}(H_{c}) (GeV/#it{c});pairs";
                hPairME[{pdgC, pdgL, "part", origin}] = new TH2F(name, title, 2000, 0., 2., 10, 0, 10);

                name = Form("hPairME_%d_%d_%s", pdgC, -pdgL, origin.data());
                title = ";#it{k}* (GeV/#it{c});#it{p}_{T}(H_{c}) (GeV/#it{c});pairs";
                hPairME[{pdgC, pdgL, "antipart", origin}] = new TH2F(name, title, 2000, 0., 2., 10, 0, 10);
            }
        }
    }

    std::vector<ROOT::Math::PxPyPzMVector> partLF{};
    std::map<std::string, std::vector<ROOT::Math::PxPyPzMVector>> partCharm{{"prompt", {}}, {"nonprompt", {}}};
    std::map<std::string, std::vector<ROOT::Math::PxPyPzMVector>> partCharmOriginal{{"prompt", {}}, {"nonprompt", {}}};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBufferLF{};
    std::map<std::string, std::deque<std::vector<ROOT::Math::PxPyPzMVector>>> partBufferCharm{
        {"prompt", {}},
        {"nonprompt", {}},
    };
    std::map<std::string, std::vector<int>> idxCharm{{"prompt", {}}, {"nonprompt", {}}};
    std::vector<std::vector<int>> motherLF{};
    std::vector<int> idxLF{};
    std::vector<int> pdgLF{};
    std::map<std::string, std::vector<int>> pdgCharm{{"prompt", {}}, {"nonprompt", {}}};
    std::deque<std::vector<int>> pdgBufferLF{};
    std::map<std::string, std::deque<std::vector<int>>> pdgBufferCharm{{"prompt", {}}, {"nonprompt", {}}};

    TClonesArray* particles = new TClonesArray("TParticle", 1000);

    TFile* inFile = TFile::Open(inFileName.data());
    if (!inFile) {
        printf("The file %s does not exist. Skip!", inFileName.data());
        return;
    }

    TTree* tEvents = (TTree*)inFile->Get("tEvents");
    tEvents->SetBranchStatus("particles", 1);
    tEvents->SetBranchAddress("particles", &particles);

    for (int iEvent = 0; iEvent < tEvents->GetEntries(); iEvent++) {
        tEvents->GetEntry(iEvent);

        partCharm["prompt"].clear();
        partCharm["nonprompt"].clear();
        pdgCharm["prompt"].clear();
        pdgCharm["nonprompt"].clear();
        idxCharm["prompt"].clear();
        idxCharm["nonprompt"].clear();
        partLF.clear();
        pdgLF.clear();
        motherLF.clear();

        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
            if (std::abs(particle->Y()) > 2)  // keep only midrapidity
                continue;

            double prodR = std::sqrt(particle->Vx() * particle->Vx() + particle->Vy() * particle->Vy() +
                                     particle->Vz() * particle->Vz());

            if (prodR > 1.)  // to remove products of strange decays
                continue;

            int pdg = particle->GetPdgCode();
            int absPdg = std::abs(pdg);
            if (std::find(charmPDG.begin(), charmPDG.end(), absPdg) != charmPDG.end()) {
                if (selDauKinem) {
                    std::multiset<int> dauPdgs = {};
                    // first selection on BR
                    if (particle->GetNDaughters() != decayChannels[pdg].size()) {
                        continue;
                    }

                    // species dependent selections
                    if (absPdg == 411) {  // select decay kinem of Dplus
                        for (int iDau = particle->GetFirstDaughter(); iDau <= particle->GetLastDaughter(); iDau++) {
                            TParticle* dau = dynamic_cast<TParticle*>(particles->At(iDau));

                            // kinem selections
                            if (!isInTPC(dau)) {  // same selections as in data
                                goto nextparticle;
                            }

                            // prepare BR selection
                            dauPdgs.insert(dau->GetPdgCode());
                        }
                        // BR selection
                        if (dauPdgs != decayChannels[pdg]) {  // select the interesting decay channel
                            goto nextparticle;
                        }
                    } else if (absPdg == 413) {  // selct decay kinem of Dstar
                        // first: Dstarplus -> D0 Pi+
                        for (int iDau = particle->GetFirstDaughter(); iDau <= particle->GetLastDaughter(); iDau++) {
                            TParticle* dau = dynamic_cast<TParticle*>(particles->At(iDau));
                            int dauPdg = dau->GetPdgCode();
                            if (std::abs(dauPdg) == 211 && (dau->Pt() < 0.3 || std::abs(dau->Eta())) > 0.8) {
                                goto nextparticle;
                            }

                            // kinem selections on D0 daus
                            if (std::abs(dauPdg) == 421) {
                                // first selection on BR
                                if (dau->GetNDaughters() != decayChannels[dauPdg].size()) {
                                    goto nextparticle;
                                }

                                // select kinem of D0 dau
                                std::multiset<int> D0dauPdgs = {};
                                for (int iD0Dau = dau->GetFirstDaughter(); iD0Dau <= dau->GetLastDaughter(); iD0Dau++) {
                                    TParticle* D0dau = dynamic_cast<TParticle*>(particles->At(iD0Dau));
                                    if (!isInTPC(D0dau)) {
                                        goto nextparticle;
                                    }
                                    dauPdgs.insert(D0dau->GetPdgCode());
                                }

                                // select decay channel
                                if (D0dauPdgs != decayChannels[dauPdg]) {  // select the interesting decay channel
                                    goto nextparticle;
                                }
                            } else if (std::abs(dauPdg) == 211) {
                                if (!isInTPC(dau)) {  // same selections as in data
                                    goto nextparticle;
                                }
                            } else {  // neither pion nor D0: skip
                                goto nextparticle;
                            }
                            dauPdgs.insert(dau->GetPdgCode());
                        }  // end of loop over Dstar daus
                    }

                    // selection on decay channel
                    if (dauPdgs != decayChannels[pdg]) goto nextparticle;
                }  // sel dau kinem

                bool isFromB = false;
                int motherIdx = particle->GetFirstMother();
                while (motherIdx > 1) {  // 0 and 1 protons
                    TParticle* mom = dynamic_cast<TParticle*>(particles->At(motherIdx));
                    int absPdgMom = std::abs(mom->GetPdgCode());
                    if (absPdgMom == 5 || absPdgMom / 100 == 5 || absPdgMom / 1000 == 5 ||
                        (absPdgMom - 10000) / 100 == 5 || (absPdgMom - 20000) / 100 == 5 ||
                        (absPdgMom - 30000) / 100 == 5 || (absPdgMom - 100000) / 100 == 5 ||
                        (absPdgMom - 200000) / 100 == 5 ||
                        (absPdgMom - 300000) / 100 == 5) {  // remove beauty feed-down
                        isFromB = true;
                        break;
                    }
                    motherIdx = mom->GetFirstMother();
                }
                double mass = TDatabasePDG::Instance()->GetParticle(absPdg)->Mass();
                ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(), mass);
                std::string origin = isFromB ? "nonprompt" : "prompt";
                partCharm[origin].push_back(part);
                pdgCharm[origin].push_back(pdg);
                idxCharm[origin].push_back(iPart);
                // isCharmFromB.push_back(isFromB);
            } else if (std::find(lightPDG.begin(), lightPDG.end(), absPdg) != lightPDG.end()) {
                std::vector<int> mothers{};
                int motherIdx = particle->GetFirstMother();

                if (std::abs(particle->Eta()) > 0.8) continue;
                if (std::abs(particle->Pt()) < 0.3) continue;

                while (motherIdx > 1) {  // 0 and 1 protons
                    mothers.push_back(motherIdx);
                    TParticle* mom = dynamic_cast<TParticle*>(particles->At(motherIdx));
                    motherIdx = mom->GetFirstMother();
                }
                motherLF.push_back(mothers);
                double mass = TDatabasePDG::Instance()->GetParticle(absPdg)->Mass();
                ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(), mass);
                partLF.push_back(part);
                pdgLF.push_back(pdg);
            }
        nextparticle:
            continue;
        }
        if (align) {
            // Copy the original particles otherwise the calculation of the pt(D) is wrong
            partCharmOriginal["prompt"] = partCharm["prompt"];
            partCharmOriginal["nonprompt"] = partCharm["nonprompt"];

            // Perform alignment of the events
            std::array<std::array<double, 3>, 3> rot = GetRotationSpheri3DFull(particles);
            partLF = Rotate(rot, partLF);
            partCharm["prompt"] = Rotate(rot, partCharm["prompt"]);
            partCharm["nonprompt"] = Rotate(rot, partCharm["nonprompt"]);
        }

        partBufferLF.push_back(partLF);
        partBufferCharm["prompt"].push_back(partCharm["prompt"]);
        partBufferCharm["nonprompt"].push_back(partCharm["nonprompt"]);
        pdgBufferLF.push_back(pdgLF);
        pdgBufferCharm["prompt"].push_back(pdgCharm["prompt"]);
        pdgBufferCharm["nonprompt"].push_back(pdgCharm["nonprompt"]);
        if (partBufferLF.size() > 10) {  // buffer full, let's kill the first entry
            partBufferLF.pop_front();
            partBufferCharm["prompt"].pop_front();
            partBufferCharm["nonprompt"].pop_front();
            pdgBufferLF.pop_front();
            pdgBufferCharm["prompt"].pop_front();
            pdgBufferCharm["nonprompt"].pop_front();
        }

        // same event
        for (std::string origin : {"prompt", "nonprompt"}) {
            for (size_t iCharm = 0; iCharm < partCharm[origin].size(); iCharm++) {
                float ptDmeson = align ? partCharmOriginal[origin][iCharm].Pt() : partCharm[origin][iCharm].Pt();
                for (size_t iLF = 0; iLF < partLF.size(); iLF++) {
                    auto moms = motherLF[iLF];
                    if (std::find(moms.begin(), moms.end(), idxCharm[origin][iCharm]) != moms.end()) {
                        continue;
                    }

                    double kStar = ComputeKstar(partCharm[origin][iCharm], partLF[iLF]);

                    if (pdgCharm[origin][iCharm] * pdgLF[iLF] > 0) {
                        int charmPdg = std::abs(pdgCharm[origin][iCharm]);
                        hPairSE[{charmPdg, std::abs(pdgLF[iLF]), "part", origin}]->Fill(kStar, ptDmeson);
                    } else {
                        int charmPdg = std::abs(pdgCharm[origin][iCharm]);
                        hPairSE[{charmPdg, std::abs(pdgLF[iLF]), "antipart", origin}]->Fill(kStar, ptDmeson);
                    }
                }
            }
        }

        // mixed event
        if (partBufferLF.size() < 2)  // to avoid repetitions
            continue;

        for (std::string origin : {"prompt", "nonprompt"}) {
            auto charmBuffer = partBufferCharm[origin];
            for (size_t iCharm = 0; iCharm < charmBuffer[charmBuffer.size() - 1].size(); iCharm++) {  // last only
                for (size_t iME = 0; iME < partBufferLF.size() - 1; iME++) {  // from 0 to last-1
                    for (size_t iLF = 0; iLF < partBufferLF[iME].size(); iLF++) {
                        double kStar =
                            ComputeKstar(charmBuffer[charmBuffer.size() - 1][iCharm], partBufferLF[iME][iLF]);
                        float ptDmeson = charmBuffer[charmBuffer.size() - 1][iCharm].Pt();

                        if (pdgBufferCharm[origin][charmBuffer.size() - 1][iCharm] * pdgBufferLF[iME][iLF] > 0) {
                            int charmPdg = std::abs(pdgBufferCharm[origin][charmBuffer.size() - 1][iCharm]);
                            int std::abs(pdgBufferLF[iME][iLF]);
                            hPairME[{charmPdg, lightPdg, "part", origin}]->Fill(kStar, ptDmeson);
                        } else {
                            int charmPdg = std::abs(pdgBufferCharm[origin][charmBuffer.size() - 1][iCharm]);
                            int std::abs(pdgBufferLF[iME][iLF]);
                            hPairME[{charmPdg, lightPdg, "antipart", origin}]->Fill(kStar, ptDmeson);
                        }
                    }
                }
            }
        }
    }

    // save root output file
    TFile outFile(oFileName.data(), "create");
    for (auto& pdgC : charmPDG) {
        for (auto& pdgL : lightPDG) {
            for (std::string origin : {"prompt", "nonprompt"}) {
                hPairSE[{pdgC, pdgL, "part", origin}]->Write();
                hPairSE[{pdgC, pdgL, "antipart", origin}]->Write();
                hPairME[{pdgC, pdgL, "part", origin}]->Write();
                hPairME[{pdgC, pdgL, "antipart", origin}]->Write();
            }
        }
    }
    outFile.Close();
}
