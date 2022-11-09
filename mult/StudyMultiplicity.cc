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

#include <array>
#include <deque>
#include <map>
#include <string>
#include <vector>

#include "AliPythia8.h"
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"


using namespace Pythia8;

namespace {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };

enum processes { kSoftQCD = 0, kHardQCD };

}  // namespace

void StudyMultiplicity(int nEvents, float maxRunTime, tunes tune, processes process, int seed, std::string oDir);
bool IsDetectable(int absPdg);


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

void StudyMultiplicity(int nEvents, float maxRunTime, tunes tune, processes process, int seed, std::string oDir) {
    TStopwatch timer;
    timer.Start();

    int energy = 13000;
    int collidingParticle1Pdg = 2212;
    int collidingParticle2Pdg = 2212;

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

    // save root output file
    std::string oFileName = "MultDistr";
    if (collidingParticle1Pdg == 2212 && collidingParticle2Pdg == 2212) {
        oFileName += "_pp";
    } else {
        printf("Error: colliding system not implemented! Exit!\n");
        return;
    }

    if (energy == 13000) {
        oFileName += "13TeV";
    } else {
        printf("Error: energy not implemented! Exit!\n");
        return;
    }

    oFileName += ".root";
    //__________________________________________________________
    // define outputs
    TH2I *hMult = new TH2I("hMult", "hMult;particles in V0A or V0C;Partices in |#eta|<0.8;Entries", 201, -0.5, 200.5, 201, -0.5, 200.5);
    TH2I *hEtaMultForward = new TH2I("hEtahEtaMultForward", "hEtaMultForward;#eta;Particles in V0A or V0C;Entries", 100, -4, 6, 201, -0.5, 200.5);
    TH2I *hEtaMultMid = new TH2I("hEtahEtaMultMid", "hEtaMultMid;#eta;Particles in V0A or V0C;Entries", 100, -1, 1, 201, -0.5, 200.5);

    TClonesArray* particles = new TClonesArray("TParticle", 1000);

    std::vector<float> etasForward;
    std::vector<float> etasMid;

    for (int iEvent = 0; iEvent < nEvents; iEvent++) {
        printf("%f, %f\n", timer.RealTime(), maxRunTime*3600);
        if (timer.RealTime() > maxRunTime*3600){
            int time = int(timer.RealTime()); // in seconds
            int nHours = int(time/3600);
            int nMins = int((time%3600)/60);

            printf("Reached max run time: %d:%d\n", nHours, nMins);
            break;
        } else {
            timer.Continue();
        }

        etasForward.clear();
        etasMid.clear();
        
        pythia.GenerateEvent();
        pythia.ImportParticles(particles, "All");

        // evaluate multiplicity at forward rapidity
        int nChForward = 0;
        int nChMid = 0;
        for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
            int pdg = std::abs(particle->GetPdgCode());
            int status = std::abs(particle->GetStatusCode());
            float eta = particle->Eta();

            if (IsDetectable(pdg) && std::abs(eta) < 0.8 && status == 1) {
                etasMid.push_back(eta);
                nChMid++;
            } else if (IsDetectable(pdg) && ((-3.7 < eta && eta < -1.7) || (2.8 < eta && eta < 5.1))  && status == 1) { // V0A and V0C acceptance
                etasForward.push_back(eta);
                nChForward++;
            }
        }
        if (nChForward > 2 || nChMid > 2) {
            for (float eta : etasForward) {
                hEtaMultForward->Fill(eta, nChForward);
            }
            for (float eta : etasMid) {
                hEtaMultMid->Fill(eta, nChMid);
            }
            hMult->Fill(nChForward, nChMid);
        }
    }

    TFile outFile(Form("%s/%s", oDir.data(), oFileName.data()), "recreate");
    hMult->Write();
    hEtaMultForward->Write();
    hEtaMultMid->Write();
    outFile.Close();


}