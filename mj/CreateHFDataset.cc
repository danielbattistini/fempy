#include <TClonesArray.h>
#include <TFile.h>
#include <TStopwatch.h>

#include <string>
#include <vector>

#include "AliPythia8.h"

using namespace Pythia8;

namespace {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
enum processes { kSoftQCD = 0, kHardQCD };
enum triggers { kMB = 0, kHighMultV0 };
enum kinem { kAny = 0, kDauInEta08 };
}  // namespace

std::map<int, std::multiset<int>> decayChannels = {
    // D
    {+411, {-321, +211, +211}},
    {-411, {+321, -211, -211}},
    // Dzero
    {+421, {-321, +211}},
    {-421, {+321, -211}},
    // Dstar
    {+413, {+421, +211}},
    {-413, {-421, -211}},
};

bool IsDetectable(const int &absPdg) {
    if (absPdg == 11 ||   // electrons
        absPdg == 13 ||   // muons
        absPdg == 211 ||  // pions
        absPdg == 321 ||  // kaons
        absPdg == 2212    // protons
    )
        return true;

    return false;
}

int GetV0Mult(const TClonesArray* particles) {
    int nCh = 0;
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        int pdg = std::abs(particle->GetPdgCode());
        int status = std::abs(particle->GetStatusCode());
        float eta = particle->Eta();

        // V0A and V0C acceptance
        if (IsDetectable(pdg) && ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1)) && status == 1) {
            nCh++;
        }
    }
    return nCh;
}

inline bool isInAcc(TParticle* p) { return p->Pt() > 0.3 && std::abs(p->Eta()) < 0.8; }

bool isParticleInEvent(const TClonesArray* particles, const int &absPdg, const bool &requireAllDauInAcc = false) {
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        int currentPdg = particle->GetPdgCode();
        if (absPdg == std::abs(currentPdg)) {
            if (!requireAllDauInAcc) return true;

            if (particle->GetNDaughters() != decayChannels[currentPdg].size()) continue;

            std::multiset<int> dauPdgs = {};
            if (absPdg == 411) {  // Dplus
                for (int iDau = particle->GetFirstDaughter(); iDau <= particle->GetLastDaughter(); iDau++) {
                    TParticle* dau = dynamic_cast<TParticle*>(particles->At(iDau));

                    // kinem selections
                    if (!isInAcc(dau)) goto nextparticle;

                    // prepare BR selection
                    dauPdgs.insert(dau->GetPdgCode());
                }
                // BR selection
                if (dauPdgs != decayChannels[currentPdg]) {  // select the interesting decay channel
                    continue;
                }
                return true;
            } else if (absPdg == 413) { // selct decay kinem of Dstar
                // todo: D* kinem to be fixed
                // first: Dstarplus -> D0 Pi+
                for (int iDau = particle->GetFirstDaughter(); iDau <= particle->GetLastDaughter(); iDau++){
                    TParticle* dau = dynamic_cast<TParticle*>(particles->At(iDau));
                    int dauPdg = dau->GetPdgCode();

                    // if (std::abs(dauPdg) != 211 && std::abs(dauPdg) != 421) goto nextparticle;
                    if (std::abs(dauPdg) == 211 && !isInAcc(dau)) goto nextparticle;
                    
                    // kinem selections on D0 daus
                    if (std::abs(dauPdg) == 421) {
                        // first selection on BR
                        if (dau->GetNDaughters() != decayChannels[dauPdg].size()) {
                            goto nextparticle;
                        }

                        // select kinem of D0 dau
                        std::multiset<int>D0dauPdgs = {};
                        for (int iD0Dau = dau->GetFirstDaughter(); iD0Dau <= dau->GetLastDaughter(); iD0Dau++){
                            TParticle* D0dau = dynamic_cast<TParticle*>(particles->At(iD0Dau));
                            if (!isInAcc(D0dau)) goto nextparticle;

                            dauPdgs.insert(D0dau->GetPdgCode());
                        }

                        // select decay channel
                        if (D0dauPdgs != decayChannels[dauPdg]) { // select the interesting decay channel
                            goto nextparticle;
                        }
                    } else if (std::abs(dauPdg) == 211) {
                        if (!isInAcc(dau)) { // same selections as in data
                            goto nextparticle;
                        }
                    } else { // neither pion nor D0: skip 
                        goto nextparticle;
                    }
                    dauPdgs.insert(dau->GetPdgCode());
                    
                } // end of loop over Dstar daus
            }
        }
    nextparticle:
        continue;
    }
    return false;
}

void CreateHFDataset(int nEvents, double maxRunTime, int energy, triggers trigger, tunes tune, processes process,
                     int hfPdg, kinem kinemSel, int seed, std::string outDir) {
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
    pythia.Initialize(2212, 2212, energy);



    int fwMultThr = 0;
    if (trigger == kHighMultV0 && process == kSoftQCD && tune == Monash)
        fwMultThr = 130;
    else if (trigger == kHighMultV0 && process == kHardQCD && tune == Monash)
        fwMultThr = 133;
    else
        print("error. mult thr not implemented")

    //__________________________________________________________
    // define outputs
    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    TTree* tEvents = new TTree("tEvents", "tEvents");
    tEvents->Branch("particles", "TClonesArray", &particles);

    int nCharm = 0;
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


        // trigger selection
        if (trigger == kHighMultV0 && GetV0Mult(particles) < fwMultThr) continue;

        // particle (and kinem) selection
        if (!isParticleInEvent(particles, hfPdg, kinemSel == kDauInEta08)) continue;

        nCharm++;
        tEvents->Fill();
    }

    // save root output file
    std::string oFileName = "AnalysisResults.root";
    TFile outFile(Form("%s/%s", outDir.data(), oFileName.data()), "recreate");
    tEvents->Write();
    outFile.Close();

    printf("---> simulated %d interesting events.\n", nCharm);
}