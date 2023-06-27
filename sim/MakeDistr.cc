// Only works if compiled

#include <Riostream.h>
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

#include "AliPythia8.h"
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "functions.hxx"
#include "yaml-cpp/yaml.h"

enum AlignMode {
    kNo = 0,
    kId,
    kSpheriFull,
    kLeading,
};

using namespace Pythia8;

void MakeDistrTest(std::string inFileName = "", std::string oFileName = "Distr.root",
                   std::string cfgFile = "cfg_distr.yml", std::string cfgSel1File = "sel_lf.yml",
                   std::string cfgSel2File = "sel_lf.yml", int seed = 42) {
    // Load simulation settings
    YAML::Node cfg = YAML::LoadFile(cfgFile.data());
    unsigned int md = cfg["mixdepth"].as<int>();
    bool rejBadShape = cfg["rejbadshape"].as<bool>();
    bool cleanPairs = cfg["cleanpairs"].as<bool>();
    bool useEvtsWithPairs = cfg["useevtswithpairs"].as<bool>();
    std::pair<double, double> y = cfg["y"].as<std::pair<double, double>>();

    std::string sAlign = cfg["align"].as<std::string>();
    AlignMode align;
    if (sAlign == "kNo")
        align = AlignMode::kNo;
    else if (sAlign == "kId")
        align = AlignMode::kId;
    else if (sAlign == "kSpheriFull")
        align = AlignMode::kSpheriFull;
    else if (sAlign == "kLeading")
        align = AlignMode::kLeading;
    else {
        std::cerr << "Align mode not implemented. Exit!" << std::endl;
        return;
    }

    // Load pythia libraries
    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    AliPythia8 pythia;

    TClonesArray *particles = new TClonesArray("TParticle", 1000);
    TTree *tEvents;

    size_t nEvents;
    std::string exeType = cfg["execution"]["type"].as<std::string>();
    if (exeType == "generate") {
        nEvents = cfg["execution"]["generate"]["nevts"].as<int>();

        // Set tune
        std::string sTune = cfg["execution"]["generate"]["tune"].as<std::string>();
        tunes tune;
        if (sTune == "kMonash")
            tune = tunes::kMonash;
        else if (sTune == "kCRMode0")
            tune = tunes::kCRMode0;
        else if (sTune == "kCRMode2")
            tune = tunes::kCRMode2;
        else if (sTune == "kCRMode2")
            tune = tunes::kCRMode2;
        else {
            std::cerr << "Pythia tune '" << sTune << "' not implemented. Exit!" << std::endl;
            exit(1);
        }
        SetTune(pythia, tune);

        // Set process
        std::string sProcess = cfg["execution"]["generate"]["process"].as<std::string>();
        processes process;
        if (sProcess == "kSoftQCD")
            process = processes::kSoftQCD;
        else if (sProcess == "kHardQCD")
            process = processes::kHardQCD;
        else {
            std::cerr << "Pythia process '" << sProcess << "' not implemented. Exit!" << std::endl;
            exit(1);
        }
        SetProcess(pythia, process);

        // Init pythia
        pythia.ReadString(Form("Tune:pp = 14"));
        pythia.ReadString("Random:setSeed = on");
        pythia.ReadString(Form("Random:seed %d", seed));
        pythia.Initialize(2212, 2212, 13000);
    } else if (exeType == "analyze") {
        // Load dataset file
        TFile *inFile = TFile::Open(inFileName.data());
        if (!inFile) {
            cerr << "The file '" << inFileName << "' does not exist. Exit!" << std::endl;
            exit(1);
        }

        // Load the events
        tEvents = (TTree *)inFile->Get("tEvents");
        // tEvents->SetDirectory(0);
        // inFile->Close();

        tEvents->SetBranchStatus("particles", 1);
        tEvents->SetBranchAddress("particles", &particles);
        nEvents = tEvents->GetEntries();
    } else {
        std::cerr << "execution type '" << exeType << "' not implemented. Exit!" << std::endl;
        return;
    }

    // Load selections for particle 1
    YAML::Node cfgSel1 = YAML::LoadFile(cfgSel1File.data());
    auto pdg1 = cfgSel1["pdg"].as<int>();
    auto selPt1 = cfgSel1["pt"].as<std::pair<double, double>>();
    auto selEta1 = cfgSel1["eta"].as<std::pair<double, double>>();
    auto selY1 = cfgSel1["y"].as<std::pair<double, double>>();
    auto selInTPC1 = cfgSel1["inTPC"].as<bool>();
    auto selDauInTPC1 = cfgSel1["dauinTPC"].as<bool>();

    // Load selections for particle 2
    YAML::Node cfgSel2 = YAML::LoadFile(cfgSel2File.data());
    auto pdg2 = cfgSel2["pdg"].as<int>();
    auto selPt2 = cfgSel2["pt"].as<std::pair<double, double>>();
    auto selEta2 = cfgSel2["eta"].as<std::pair<double, double>>();
    auto selY2 = cfgSel2["y"].as<std::pair<double, double>>();
    auto selInTPC2 = cfgSel2["inTPC"].as<bool>();
    auto selDauInTPC2 = cfgSel2["dauinTPC"].as<bool>();

    std::map<std::string, TH1D *> hSE = {
        {"sc", new TH1D("hSE_sc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
        {"oc", new TH1D("hSE_oc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
    };

    std::map<std::string, TH1D *> hME = {
        {"sc", new TH1D("hME_sc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
        {"oc", new TH1D("hME_oc", ";#it{k}* (GeV/#it{c});pairs", 2000, 0., 2.)},
    };

    std::vector<FemtoParticle> part1{};
    std::vector<FemtoParticle> part2{};
    std::deque<std::vector<FemtoParticle>> partBuffer2{};

    for (size_t iEvent = 0; iEvent < nEvents; iEvent++) {
        part1.clear();
        part2.clear();

        if (exeType == "generate") {
            pythia.GenerateEvent();
            pythia.ImportParticles(particles, "All");
        } else {
            tEvents->GetEntry(iEvent);
            printf("ievt: %d    nparticles: %d\n", iEvent, particles->GetEntriesFast());
        }

        m3 mSphi = ComputeSphericityMatrix(particles);

        // check event shape
        if (rejBadShape) {
            auto eivals = ComputeEigenValues(mSphi);
            if (!(0 < eivals[0] && eivals[0] < eivals[1] && eivals[1] < eivals[2])) continue;
        }

        for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
            TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));

            int pdg = particle->GetPdgCode();
            int absPdg = std::abs(pdg);
            if (absPdg != pdg1 && absPdg != pdg2) continue;

            // Keep only mid rapidity
            if (!(y.first < particle->Y() && particle->Y() < y.second)) continue;
            // Remove products of strange decays
            if (ProductionRadius(particle) > 1.) continue;

            double px = particle->Px();
            double py = particle->Py();
            double pz = particle->Pz();
            double mass = TDatabasePDG::Instance()->GetParticle(absPdg)->Mass();
            ROOT::Math::PxPyPzMVector part(px, py, pz, mass);
            int d1 = particle->GetFirstDaughter();
            int d2 = particle->GetLastDaughter();

            if (absPdg == pdg1) {
                if (
                    !In(selPt1, particle->Pt()) ||
                    !In(selEta1, particle->Eta()) ||
                    !In(selY1, particle->Y()) ||
                    (selInTPC1 && !isInTPC(particle))
                ) continue;

                // todo: daughters in TPC
                part1.push_back({part, pdg, iPart, {d1, d2}});
            }
            if (absPdg == pdg2) {
                if (
                    !In(selPt2, particle->Pt()) ||
                    !In(selEta2, particle->Eta()) ||
                    !In(selY2, particle->Y()) ||
                    (selInTPC2 && !isInTPC(particle))
                ) continue;

                // todo: daughters in TPC
                part2.push_back({part, pdg, iPart, {d1, d2}});
            }
        }

        if (useEvtsWithPairs && (part1.size() == 0 || part2.size() == 0)) continue;

        // Perform alignment of the events
        if (align == kSpheriFull) {
            m3 rot = GetRotationSpheri3DFull(mSphi);
            part2 = Rotate(rot, part2);
            part1 = Rotate(rot, part1);
        } else if (align == kLeading) {
            TParticle pLeading(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
            for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
                TParticle *p = dynamic_cast<TParticle *>(particles->At(iPart));
                if (p->P() > pLeading.P() && isInTPC(p) && IsDetectable(std::abs(p->GetPdgCode()))) {
                    pLeading = *p;
                }
            }
            double phiLeading = atan(pLeading.Py() / pLeading.Px());
            if (pLeading.Px() < 0) phiLeading += TMath::Pi();

            m3 rot;
            rot[0] = {+cos(phiLeading), +sin(phiLeading), 0};
            rot[1] = {-sin(phiLeading), -cos(phiLeading), 0};
            rot[2] = {0, 0, 1};

            part2 = Rotate(rot, part2);
            part1 = Rotate(rot, part1);

            double thetaLeading = atan(pLeading.Pt() / pLeading.Pz());
            if (pLeading.Pz() < 0) thetaLeading += TMath::Pi();

            rot[0] = {+cos(thetaLeading), 0, -sin(thetaLeading)};
            rot[1] = {0, +1, 0};
            rot[2] = {+sin(thetaLeading), 0, +cos(thetaLeading)};

            part2 = Rotate(rot, part2);
            part1 = Rotate(rot, part1);
        }

        // same event
        for (size_t i1 = 0; i1 < part1.size(); i1++) {
            const auto p1 = part1[i1];

            // don't pair twice in case of same-particle femto
            int start = pdg1 == pdg2 ? i1 + 1 : 0;
            for (size_t i2 = start; i2 < part2.size(); i2++) {
                const auto p2 = part2[i2];

                if (cleanPairs && !IsPairClean(p1, p2)) continue;

                double kStar = ComputeKstar(p1.p, p2.p);
                std::string pair = p1.pdg * p2.pdg > 0 ? "sc" : "oc";
                hSE[pair]->Fill(kStar);
            }
        }

        // mixed event
        for (size_t i1 = 0; i1 < part1.size(); i1++) {
            const auto p1 = part1[i1];

            for (size_t iME = 0; iME < partBuffer2.size(); iME++) {
                for (size_t i2 = 0; i2 < partBuffer2[iME].size(); i2++) {
                    const auto p2 = partBuffer2[iME][i2];

                    double kStar = ComputeKstar(p1.p, p2.p);
                    std::string pair = p1.pdg * p2.pdg > 0 ? "sc" : "oc";

                    hME[pair]->Fill(kStar);
                }
            }
        }

        partBuffer2.push_back(part2);
        if (partBuffer2.size() > md) partBuffer2.pop_front();
    }

    TFile *oFile = new TFile(oFileName.data(), "recreate");
    for (auto pair : {"sc", "oc"}) {
        oFile->mkdir(pair);
        oFile->cd(pair);
        hSE[pair]->Write("hSE");
        hME[pair]->Write("hME");
    }
    oFile->Close();
}