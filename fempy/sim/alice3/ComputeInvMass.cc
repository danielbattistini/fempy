#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include "yaml-cpp/yaml.h"

#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH2F.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "ActsFatras/EventData/Barcode.hpp"
#include "/home/ktas/ge86rim/phsw/fempy/fempy/sim/functions.hxx"

std::map<int, double> masses = {
    {211, TDatabasePDG::Instance()->GetParticle(211)->Mass()},  //  pion+
    {321, TDatabasePDG::Instance()->GetParticle(321)->Mass()},  //  kaon+
    {411, TDatabasePDG::Instance()->GetParticle(411)->Mass()},  //  D+
    {421, TDatabasePDG::Instance()->GetParticle(421)->Mass()},  //  D0
    {413, TDatabasePDG::Instance()->GetParticle(413)->Mass()},  //  D*+(2010)
};

struct particle {
    int pdg;                        // pdg code
    uint64_t id;                    // index of the particle in the event
    double t_x;                     // production vertex x
    double t_y;                     // production vertex y
    double t_z;                     // production vertex z
    ROOT::Math::PxPyPzMVector p;    // reconstructed quadrimomentum
    ROOT::Math::PxPyPzMVector t_p;  // true quadrimomentum
    int mother1idx;                 // index of mother 1
    int mother2idx;                 // index of mother 2
    int mother1pdg;                 // pdg code of mother 1
};

// Returns true if the particles have adjacent indeces.
bool AreSiblings(particle p1, particle p2) { return std::abs(static_cast<int>(p1.id - p2.id)) == 1; }

// Returns true if the mass is compatible with a given hypothesis. The precision is in GeV.
bool IsMassCorrect(int pdg, double mass, double precision = 0.0001) { return std::abs(mass - masses[pdg]) < precision; }

/*Build the mother of the particles p1 and p2.
It is assumed that the particles are ordered in descending mass order and that the PDG codes are correct.*/
particle BuildMother(int pdg, particle p1, particle p2) {
    // Check if particles have the same mother. Only check for weak decays. In this case m1>0 and and m2=0
    // see: https://pythia.org/latest-manual/ParticleProperties.html
    int abspdg = std::abs(pdg);

    if (abspdg == 421) {
        if (std::abs(p1.mother1pdg) == abspdg &&
            std::abs(p2.mother1pdg) == abspdg &&
            p1.pdg * p2.pdg < 0 &&
            p1.mother1idx == p2.mother1idx &&
            p1.mother1idx > 0 &&
            p2.mother1idx > 0 &&
            p1.mother2idx == 0 &&
            p2.mother2idx == 0 &&
            AreSiblings(p1, p2)) {
            particle mother({p1.mother1pdg, (u_int64_t)p1.mother1idx, 0, 0, 0, p1.p + p2.p, p1.t_p + p2.t_p, 0, 0, 0});
            return mother;
        }
    } else if (abspdg == 413) {
        particle Dzero = p1;
        particle piSoft = p2;

        if (std::abs(piSoft.mother1pdg) == 413 || piSoft.mother1idx > 0 || piSoft.mother2idx == 0) {
            particle mother(
                {413, (u_int64_t)piSoft.mother1idx, 0, 0, 0, Dzero.p + piSoft.p, Dzero.t_p + piSoft.t_p, 0, 0, 0});
            return mother;
        }
    }

    return particle({});
}

void ComputeInvMass(const char *configFile) {
    // consts
    const double Pi = TMath::Pi();

    YAML::Node cfg = YAML::LoadFile(configFile);
    std::string inFileName = cfg["infile"].as<std::string>();
    std::string oFileName = cfg["ofile"].as<std::string>();
    int pdg1 = cfg["pdg1"].as<int>();
    int pdg2 = cfg["pdg2"].as<int>();
    int minNHits = cfg["min_hits"].as<int>();
    int md = cfg["mixing_depth"].as<int>();
    bool useMCTruth = cfg["use_mc_truth"].as<bool>();
    bool pairOnlyPrim = cfg["pair_only_prim"].as<bool>();

    TChain *tree;
    if (inFileName.find(".chain") != std::string::npos) {
        std::cout << "Loading the TChain..." << std::endl;

        tree = new TChain("tracksummary");

        std::ifstream fileList(inFileName.data());
        std::string file;
        while (fileList >> file) {
            std::cout << file << std::endl;
            tree->Add(file.data());
        }
        std::cout << "Done!" << std::endl;
    } else {
        tree = new TChain("tracksummary");
        tree->AddFile(inFileName.data());
    }

    // init to 0 otherwise it breaks
    std::vector<int> *particle_type = 0;
    std::vector<int> *mother1_particle_id = 0;
    std::vector<int> *mother2_particle_id = 0;
    std::vector<int> *mother1_pdg = 0;
    std::vector<int> *mother2_pdg = 0;
    std::vector<int> *nMeasurements = 0;
    std::vector<float> *px = 0;
    std::vector<float> *py = 0;
    std::vector<float> *pz = 0;

    std::vector<float> *qop = 0;  // charge over momentum
    std::vector<float> *phi = 0;
    std::vector<float> *theta = 0;

    std::vector<float> *vx = 0;
    std::vector<float> *vy = 0;
    std::vector<float> *vz = 0;

    // MC truth variables
    std::vector<float> *t_vx = 0;
    std::vector<float> *t_vy = 0;
    std::vector<float> *t_vz = 0;
    std::vector<float> *t_px = 0;
    std::vector<float> *t_py = 0;
    std::vector<float> *t_pz = 0;
    std::vector<int> *t_charge = 0;
    std::vector<uint64_t> *t_majorityParticleId = 0;

    tree->SetBranchAddress("eQOP_fit", &qop);
    tree->SetBranchAddress("ePHI_fit", &phi);
    tree->SetBranchAddress("eTHETA_fit", &theta);

    // MC truth vaariables
    tree->SetBranchAddress("particle_type", &particle_type);
    tree->SetBranchAddress("mother1_particle_id", &mother1_particle_id);
    tree->SetBranchAddress("mother2_particle_id", &mother2_particle_id);
    tree->SetBranchAddress("mother1_pdg", &mother1_pdg);
    tree->SetBranchAddress("mother2_pdg", &mother2_pdg);
    tree->SetBranchAddress("nMeasurements", &nMeasurements);
    tree->SetBranchAddress("t_vx", &t_vx);
    tree->SetBranchAddress("t_vy", &t_vy);
    tree->SetBranchAddress("t_vz", &t_vz);
    tree->SetBranchAddress("t_px", &t_px);
    tree->SetBranchAddress("t_py", &t_py);
    tree->SetBranchAddress("t_pz", &t_pz);
    tree->SetBranchAddress("t_charge", &t_charge);
    tree->SetBranchAddress("majorityParticleId", &t_majorityParticleId);

    TFile *oFile = new TFile(oFileName.data(), "recreate");
    double massMin;
    double massMax;
    double trueMassMin;
    double trueMassMax;
    TString titleMassDzero = "#it{M}(K,#pi) (GeV/#it{c}^{2})";
    TString titleMassDstar = "#it{M}(K,#pi,#pi) (GeV/#it{c}^{2})";
    TString titlePt = "#it{p}_{T} (GeV/#it{c})";

    std::vector<int> lightPdg = {211, 321};
    std::vector<int> charmPdg = {421, 413};
    std::vector<int> allPdg = lightPdg;
    allPdg.insert(allPdg.end(), charmPdg.begin(), charmPdg.end());

    // names
    std::map<int, const char *> pdg2name = {
        {321, "K"},
        {211, "Pi"},
        {421, "Dzero"},
        {413, "Dstar"},
    };

    std::map<int, std::vector<particle>> particles = {
        {321, std::vector<particle>()},
        {211, std::vector<particle>()},
        {421, std::vector<particle>()},
        {413, std::vector<particle>()},
    };

    // this map also includes TH2*
    std::map<int, std::map<std::string, TH1 *>> hPartProp = {
        {321, {}},
        {211, {}},
        {421, {}},
        {413, {}},
    };

    std::map<std::string, std::map<std::string, TH1 *>> hFemto = {
        {"p02_13", {}},
        {"p12_23", {}},
    };

    // Mixing buffer
    std::deque<std::vector<particle>> partBuffer2{};

    TString name;
    TString title;
    double massBinWidth = 0.1;  // MeV/c2

    for (auto abspdg : {211, 321, 421, 413}) {
        // p
        name = Form("hP_%s", pdg2name[abspdg]);
        title = ";#it{p} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hP", new TH1D(name, title, 500, 0, 50)});

        // pt
        name = Form("hPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hPt", new TH1D(name, title, 500, 0, 50)});

        // eta
        name = Form("hEta_%s", pdg2name[abspdg]);
        title = ";#eta;Counts";
        hPartProp[abspdg].insert({"hEta", new TH1D(name, title, 500, -5, 5)});

        // phi
        name = Form("hPhi_%s", pdg2name[abspdg]);
        title = ";#varphi (rad);Counts";
        hPartProp[abspdg].insert({"hPhi", new TH1D(name, title, 500, -Pi, Pi)});

        // phi vs eta
        name = Form("hPhiVsEta_%s", pdg2name[abspdg]);
        title = ";#eta;#varphi (rad);Counts";
        hPartProp[abspdg].insert({"hPhiVsEta", new TH2F(name, title, 200, -5, 5, 200, -Pi, Pi)});

        // phi vs rProd
        name = Form("hPhiVsRprod_%s", pdg2name[abspdg]);
        title = ";#it{r}_{prod} (mm);#varphi (rad);Counts";
        hPartProp[abspdg].insert({"hPhiVsRprod", new TH2F(name, title, 200, 0, 5, 200, -Pi, Pi)});

        // P resolution
        name = Form("hResolutionP_%s", pdg2name[abspdg]);
        title = ";#it{p}^{true} (GeV/#it{c});#it{p}^{reco} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionP", new TH2F(name, title, 100, 0, 10, 100, 0, 10)});

        // P resolution delta
        name = Form("hResolutionDeltaP_%s", pdg2name[abspdg]);
        title = ";#it{p}^{true} (GeV/#it{c});#it{p}^{reco} - #it{p}^{true} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaP", new TH2F(name, title, 100, 0, 10, 100, -0.5, 0.5)});

        // P resolution percentage
        name = Form("hResolutionPercP_%s", pdg2name[abspdg]);
        title = ";#it{p}^{true} (GeV/#it{c});(#it{p}^{reco} - #it{p}^{true})/#it{p}^{true}) (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercP", new TH2F(name, title, 100, 0, 10, 100, -5, 5)});

        // Pt resolution
        name = Form("hResolutionPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T}^{true} (GeV/#it{c});#it{p}_{T}^{reco} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionPt", new TH2F(name, title, 100, 0, 10, 100, 0, 10)});

        // Pt resolution delta
        name = Form("hResolutionDeltaPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T}^{true} (GeV/#it{c});#it{p}_{T}^{reco} - #it{p}_{T}^{true} (GeV/#it{c});Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaPt", new TH2F(name, title, 100, 0, 10, 100, -0.5, 0.5)});

        // Pt resolution percentage
        name = Form("hResolutionPercPt_%s", pdg2name[abspdg]);
        title = ";#it{p}_{T}^{true} (GeV/#it{c});(#it{p}_{T}^{reco} - #it{p}_{T}^{true})/#it{p}_{T}^{true} (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercPt", new TH2F(name, title, 100, 0, 10, 100, -5, 5)});

        // Eta resolution
        name = Form("hResolutionEta_%s", pdg2name[abspdg]);
        title = ";#eta^{true};#eta^{reco};Counts";
        hPartProp[abspdg].insert({"hResolutionEta", new TH2F(name, title, 200, -5, 5, 200, -5, 5)});

        // Eta resolution delta
        name = Form("hResolutionDeltaEta_%s", pdg2name[abspdg]);
        title = ";#eta^{true};#eta^{reco} - #eta^{true};Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaEta", new TH2F(name, title, 200, -5, 5, 200, -0.1, 0.1)});

        // Eta resolution percentage
        name = Form("hResolutionPercEta_%s", pdg2name[abspdg]);
        title = ";#eta^{true};(#eta^{reco} - #eta^{true})/#eta^{true} (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercEta", new TH2F(name, title, 200, -5, 5, 200, -5, 5)});

        // Phi resolution
        name = Form("hResolutionPhi_%s", pdg2name[abspdg]);
        title = ";#phi^{true};#phi^{reco};Counts";
        hPartProp[abspdg].insert({"hResolutionPhi", new TH2F(name, title, 200, -Pi, Pi, 200, -Pi, Pi)});

        // Phi resolution delta
        name = Form("hResolutionDeltaPhi_%s", pdg2name[abspdg]);
        title = ";#phi^{true};#phi^{reco} - #phi^{true};Counts";
        hPartProp[abspdg].insert({"hResolutionDeltaPhi", new TH2F(name, title, 200, -Pi, Pi, 200, -Pi, Pi)});

        // Phi resolution percentage
        name = Form("hResolutionPercPhi_%s", pdg2name[abspdg]);
        title = ";#phi^{true};(#phi^{reco} - #phi^{true})/#phi^{true} (%);Counts";
        hPartProp[abspdg].insert({"hResolutionPercPhi", new TH2F(name, title, 200, -Pi, Pi, 200, -Pi, Pi)});

        if (std::find(charmPdg.begin(), charmPdg.end(), abspdg) != charmPdg.end()) {
            TString titleMass;
            int nMassBins;
            if (abspdg == 421) {
                titleMass = titleMassDzero;
                massMin = 1.5;
                massMax = 2.2;
                trueMassMin = 1.860;
                trueMassMax = 1.870;
                nMassBins = (int)std::round((massMax - massMin) * 1000 / massBinWidth);
            } else if (abspdg == 413) {
                titleMass = titleMassDstar;
                massMin = 1.6;
                massMax = 2.4;
                trueMassMin = 2.005;
                trueMassMax = 2.015;
                nMassBins = (int)std::round((massMax - massMin) * 1000 / massBinWidth);
            } else {
                exit(1);
            }

            // invariant mass
            name = Form("hInvMass_%s", pdg2name[abspdg]);
            title = ";" + titleMass + ";Counts";
            hPartProp[abspdg].insert({"hInvMass", new TH1D(name, title, nMassBins, massMin, massMax)});

            // True invariant mass
            name = Form("hTrueInvMass_%s", pdg2name[abspdg]);
            title = ";" + titleMass + ";Counts";
            hPartProp[abspdg].insert({"hTrueInvMass", new TH1D(name, title, nMassBins, trueMassMin, trueMassMax)});

            // invariant mass vs pt
            name = Form("hInvMassVsPt_%s", pdg2name[abspdg]);
            title = ";" + titlePt + ";" + titleMass + ";Counts";
            hPartProp[abspdg].insert({"hInvMassVsPt", new TH2F(name, title, 100, 0, 10, nMassBins, massMin, massMax)});
        }
    }

    // SE
    hFemto["p02_13"].insert({"hSE", new TH1D("hSE_02_13", ";#it{k}* (GeV/#it{c});Counts", 3000, 0, 3)});
    hFemto["p02_13"].insert({"hME", new TH1D("hME_02_13", ";#it{k}* (GeV/#it{c});Counts", 3000, 0, 3)});

    // ME
    hFemto["p12_23"].insert({"hSE", new TH1D("hSE_12_23", ";#it{k}* (GeV/#it{c});Counts", 3000, 0, 3)});
    hFemto["p12_23"].insert({"hME", new TH1D("hME_12_23", ";#it{k}* (GeV/#it{c});Counts", 3000, 0, 3)});

    for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
        tree->GetEntry(iEvent);

        particles[321].clear();
        particles[211].clear();
        particles[421].clear();
        particles[413].clear();

        size_t nPart = particle_type->size();
        for (size_t iPart = 0; iPart < nPart; iPart++) {
            int pdg = (*particle_type)[iPart];
            int abspdg = std::abs(pdg);
            if (abspdg != 211 && abspdg != 321) continue;

            int nHits = (*nMeasurements)[iPart];
            if (nHits < minNHits) continue;

            double t_ppx = (*t_px)[iPart];
            double t_ppy = (*t_py)[iPart];
            double t_ppz = (*t_pz)[iPart];
            uint64_t t_pmajorityParticleId = (*t_majorityParticleId)[iPart];

            double pqop = (*qop)[iPart];
            double pphi = (*phi)[iPart];
            double ptheta = (*theta)[iPart];

            double ppx = 1. / std::abs(pqop) * sin(ptheta) * cos(pphi);
            double ppy = 1. / std::abs(pqop) * sin(ptheta) * sin(pphi);
            double ppz = 1. / std::abs(pqop) * cos(ptheta);
            double pm = masses[abspdg];

            uint64_t id = ActsFatras::Barcode(t_pmajorityParticleId).particle();

            auto p = ROOT::Math::PxPyPzMVector(ppx, ppy, ppz, pm);
            auto t_p = ROOT::Math::PxPyPzMVector(t_ppx, t_ppy, t_ppz, pm);

            double x = (*t_vx)[iPart];  // u. m. = mm (should be)
            double y = (*t_vy)[iPart];  // u. m. = mm (should be)
            double z = (*t_vz)[iPart];  // u. m. = mm (should be)

            // reject decay products of strange decays. This selection also removes the peak in eta-phi for pions
            double rProd = pow(x * x + y * y + z * z, 0.5);

            int m1idx = (*mother1_particle_id)[iPart];
            int m2idx = (*mother2_particle_id)[iPart];

            int m1pdg = (*mother1_pdg)[iPart];
            particle part = particle({pdg, id, x, y, z, p, t_p, m1idx, m2idx, m1pdg});
            particles[abspdg].push_back(part);

            // Fill hists
            hPartProp[abspdg]["hP"]->Fill(p.P());
            hPartProp[abspdg]["hPt"]->Fill(p.Pt());
            hPartProp[abspdg]["hEta"]->Fill(p.Eta());
            hPartProp[abspdg]["hPhi"]->Fill(p.Phi());
            hPartProp[abspdg]["hPhiVsEta"]->Fill(p.Eta(), p.Phi());
            hPartProp[abspdg]["hPhiVsRprod"]->Fill(rProd, p.Phi());

            hPartProp[abspdg]["hResolutionEta"]->Fill(t_p.Eta(), p.Eta());
            hPartProp[abspdg]["hResolutionDeltaEta"]->Fill(t_p.Eta(), p.Eta() - t_p.Eta());
            hPartProp[abspdg]["hResolutionPercEta"]->Fill(t_p.Eta(), (p.Eta() - t_p.Eta()) / t_p.Eta());

            hPartProp[abspdg]["hResolutionPhi"]->Fill(t_p.Phi(), p.Phi());
            hPartProp[abspdg]["hResolutionDeltaPhi"]->Fill(t_p.Phi(), p.Phi() - t_p.Phi());
            hPartProp[abspdg]["hResolutionPercPhi"]->Fill(t_p.Phi(), (p.Phi() - t_p.Phi()) / t_p.Phi());

            hPartProp[abspdg]["hResolutionP"]->Fill(t_p.P(), p.P());
            hPartProp[abspdg]["hResolutionDeltaP"]->Fill(t_p.P(), p.P() - t_p.P());
            hPartProp[abspdg]["hResolutionPercP"]->Fill(t_p.P(), (p.P() - t_p.P()) / t_p.P() * 100);

            hPartProp[abspdg]["hResolutionPt"]->Fill(t_p.Pt(), p.Pt());
            hPartProp[abspdg]["hResolutionDeltaPt"]->Fill(t_p.Pt(), p.Pt() - t_p.Pt());
            hPartProp[abspdg]["hResolutionPercPt"]->Fill(t_p.Pt(), (p.Pt() - t_p.Pt()) / t_p.Pt() * 100);
        }

        // Reconstruct Dzero
        auto kaons = particles[321];
        auto pions = particles[211];

        for (size_t iK = 0; iK < kaons.size(); iK++) {
            auto K = kaons[iK];

            for (size_t iPi = 0; iPi < pions.size(); iPi++) {
                auto Pi = pions[iPi];
                auto Dzero = BuildMother(421, K, Pi);
                if (std::abs(Dzero.pdg) != 421 || !IsMassCorrect(421, Dzero.t_p.M())) continue;
                particles[421].push_back(Dzero);

                double m = Dzero.p.M();
                double eta = Dzero.p.Eta();
                double phi = Dzero.p.Phi();
                double p = Dzero.p.P();
                double pt = Dzero.p.Pt();

                double t_m = Dzero.t_p.M();
                double t_eta = Dzero.t_p.Eta();
                double t_phi = Dzero.t_p.Phi();
                double t_p = Dzero.t_p.P();
                double t_pt = Dzero.t_p.Pt();

                // kinematics
                hPartProp[421]["hP"]->Fill(p);
                hPartProp[421]["hPt"]->Fill(pt);
                hPartProp[421]["hEta"]->Fill(eta);
                hPartProp[421]["hPhi"]->Fill(phi);
                hPartProp[421]["hPhiVsEta"]->Fill(eta, phi);

                // Invariant mass
                hPartProp[421]["hInvMass"]->Fill(m);
                hPartProp[421]["hTrueInvMass"]->Fill(t_m);
                hPartProp[421]["hInvMassVsPt"]->Fill(pt, m);

                // Resolution
                hPartProp[421]["hResolutionEta"]->Fill(t_eta, eta);
                hPartProp[421]["hResolutionDeltaEta"]->Fill(t_eta, eta - t_eta);
                hPartProp[421]["hResolutionPercEta"]->Fill(t_eta, (eta - t_eta) / t_eta * 100);

                hPartProp[421]["hResolutionPhi"]->Fill(t_phi, phi);
                hPartProp[421]["hResolutionDeltaPhi"]->Fill(t_phi, phi - t_phi);
                hPartProp[421]["hResolutionPercPhi"]->Fill(t_phi, (phi - t_phi) / t_phi * 100);

                hPartProp[421]["hResolutionP"]->Fill(t_p, p);
                hPartProp[421]["hResolutionDeltaP"]->Fill(t_p, p - t_p);
                hPartProp[421]["hResolutionPercP"]->Fill(t_p, (p - t_p) / t_p * 100);

                hPartProp[421]["hResolutionPt"]->Fill(t_pt, pt);
                hPartProp[421]["hResolutionDeltaPt"]->Fill(t_pt, pt - t_pt);
                hPartProp[421]["hResolutionPercPt"]->Fill(t_pt, (pt - t_pt) / t_pt * 100);

                // Reconstruct D*(2010)
                for (size_t iPi2 = iPi + 1; iPi2 < pions.size(); iPi2++) {
                    auto Pi2 = pions[iPi2];

                    auto Dstar = BuildMother(413, Dzero, Pi2);
                    if (std::abs(Dstar.pdg) != 413 || !IsMassCorrect(413, Dstar.t_p.M())) continue;

                    particles[413].push_back(Dstar);

                    double m = Dstar.p.M();
                    double eta = Dstar.p.Eta();
                    double phi = Dstar.p.Phi();
                    double p = Dstar.p.P();
                    double pt = Dstar.p.Pt();

                    double t_m = Dstar.t_p.M();
                    double t_eta = Dstar.t_p.Eta();
                    double t_phi = Dstar.t_p.Phi();
                    double t_p = Dstar.t_p.P();
                    double t_pt = Dstar.t_p.Pt();

                    // kinematics
                    hPartProp[413]["hP"]->Fill(p);
                    hPartProp[413]["hPt"]->Fill(pt);
                    hPartProp[413]["hEta"]->Fill(eta);
                    hPartProp[413]["hPhi"]->Fill(phi);
                    hPartProp[413]["hPhiVsEta"]->Fill(eta, phi);

                    // Invariant mass
                    hPartProp[413]["hInvMass"]->Fill(m);
                    hPartProp[413]["hTrueInvMass"]->Fill(t_m);
                    hPartProp[413]["hInvMassVsPt"]->Fill(pt, m);

                    // Resolution
                    hPartProp[413]["hResolutionEta"]->Fill(t_eta, eta);
                    hPartProp[413]["hResolutionDeltaEta"]->Fill(t_eta, eta - t_eta);
                    hPartProp[413]["hResolutionPercEta"]->Fill(t_eta, (eta - t_eta) / t_eta * 100);

                    hPartProp[413]["hResolutionPhi"]->Fill(t_phi, phi);
                    hPartProp[413]["hResolutionDeltaPhi"]->Fill(t_phi, phi - t_phi);
                    hPartProp[413]["hResolutionPercPhi"]->Fill(t_phi, (phi - t_phi) / t_phi * 100);

                    hPartProp[413]["hResolutionP"]->Fill(t_p, p);
                    hPartProp[413]["hResolutionDeltaP"]->Fill(t_p, p - t_p);
                    hPartProp[413]["hResolutionPercP"]->Fill(t_p, (p - t_p) / t_p * 100);

                    hPartProp[413]["hResolutionPt"]->Fill(t_pt, pt);
                    hPartProp[413]["hResolutionDeltaPt"]->Fill(t_pt, pt - t_pt);
                    hPartProp[413]["hResolutionPercPt"]->Fill(t_pt, (pt - t_pt) / t_pt * 100);
                }  // loop over pions2
            }      // loop over pions
        }          // loop over kaons

        if (pairOnlyPrim) {
            for (auto abspdg : {211, 321}) {
                particles[abspdg].erase(std::remove_if(
                    particles[abspdg].begin(),
                    particles[abspdg].end(),
                    [](const particle &p) {
                        return p.mother1idx <= 0 || p.mother2idx <= 0 || p.mother1idx > p.mother2idx;
                    }), particles[abspdg].end());
            }
        }

        // Same event
        for (size_t i1 = 0; i1 < particles[pdg1].size(); i1++) {
            const auto p1 = particles[pdg1][i1];

            // Don't pair twice in case of same-particle femto
            int start = pdg1 == pdg2 ? i1 + 1 : 0;
            for (size_t i2 = start; i2 < particles[pdg2].size(); i2++) {
                const auto p2 = particles[pdg2][i2];

                // todo implement pair cleaner here

                double kStar = useMCTruth ? ComputeKstar(p1.t_p, p2.t_p) : ComputeKstar(p1.p, p2.p);
                std::string pair = p1.pdg * p2.pdg > 0 ? "p02_13" : "p12_23";
                hFemto[pair]["hSE"]->Fill(kStar);
            }
        }

        // Mixed event
        for (size_t i1 = 0; i1 < particles[pdg1].size(); i1++) {
            const auto p1 = particles[pdg1][i1];

            for (size_t iME = 0; iME < partBuffer2.size(); iME++) {
                for (size_t i2 = 0; i2 < partBuffer2[iME].size(); i2++) {
                    const auto p2 = partBuffer2[iME][i2];

                    double kStar = useMCTruth ? ComputeKstar(p1.t_p, p2.t_p) : ComputeKstar(p1.p, p2.p);
                    std::string pair = p1.pdg * p2.pdg > 0 ? "p02_13" : "p12_23";

                    hFemto[pair]["hME"]->Fill(kStar);
                }
            }
        }
        partBuffer2.push_back(particles[pdg2]);
        if (partBuffer2.size() > md) partBuffer2.pop_front();
    }  // event loop

    // Write histograms
    for (auto pdg : allPdg) {
        oFile->mkdir(pdg2name[pdg]);
        oFile->cd(pdg2name[pdg]);

        // change the name of the qa histograms and save them to file
        for (auto const &[name, hist] : hPartProp[pdg]) {
            hist->SetName(name.data());
            hist->Write();
        }
    }

    for (auto const &[name, hist] : hFemto) {
        oFile->mkdir(name.data());
        oFile->cd(name.data());

        hFemto[name]["hSE"]->SetName("hSE");
        hFemto[name]["hSE"]->Write();
        hFemto[name]["hME"]->SetName("hME");
        hFemto[name]["hME"]->Write();
    }

    oFile->Close();
    printf("Output saved in %s\n", oFileName.data());
}
