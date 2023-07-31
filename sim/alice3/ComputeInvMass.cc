#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "ActsFatras/EventData/Barcode.hpp"

template <typename T>
void print(std::vector<T> vec) {
    for (const auto &v : vec) {
        std::cout << v << "  ";
    }
    std::cout << std::endl;
}

template <typename T>
void print(std::vector<T> *vec) {
    for (const auto &v : *vec) {
        std::cout << v << "  ";
    }
    std::cout << std::endl;
}

struct particle {
    int pdg;
    uint64_t id;
    ROOT::Math::PxPyPzMVector p;
    ROOT::Math::PxPyPzMVector t_p;
};

bool AreSiblings(particle p1, particle p2) {
    return std::abs(static_cast<int>(p1.id - p2.id)) == 1;
}

static auto pdgdb = TDatabasePDG::Instance();

void ComputeInvMass(const char *inFileName, const char *oFileName, int pdg1, int pdg2, double rProdMin, double rProdMax) {
    TFile *inFile = new TFile(inFileName);

    TList *keys = (TList *)inFile->GetListOfKeys();
    if (keys->GetEntries() > 1) {
        printf("Warning: more than 1 tree found. Check\n");
    }

    // init to 0 otherwise it breaks
    std::vector<int> *particle_type = 0;
    std::vector<double> *px = 0;
    std::vector<double> *py = 0;
    std::vector<double> *pz = 0;
    
    std::vector<double> *qop = 0; // charge over momentum
    std::vector<double> *phi = 0;
    std::vector<double> *theta = 0;

    std::vector<double> *vx = 0;
    std::vector<double> *vy = 0;
    std::vector<double> *vz = 0;

    // MC truth variables
    std::vector<double> *t_vx = 0;
    std::vector<double> *t_vy = 0;
    std::vector<double> *t_vz = 0;
    std::vector<double> *t_px = 0;
    std::vector<double> *t_py = 0;
    std::vector<double> *t_pz = 0;
    std::vector<int> *t_charge = 0;
    std::vector<uint64_t> *t_majorityParticleId = 0;

    auto treeName = std::string(keys->At(0)->GetName());
    TTree *tree = (TTree *)inFile->Get(treeName.data());
    tree->SetBranchAddress("eQOP_fit", &qop);
    tree->SetBranchAddress("ePHI_fit", &phi);
    tree->SetBranchAddress("eTHETA_fit", &theta);

    // MC truth vaariables
    tree->SetBranchAddress("particle_type", &particle_type);  
    tree->SetBranchAddress("t_vx", &t_vx);
    tree->SetBranchAddress("t_vy", &t_vy);
    tree->SetBranchAddress("t_vz", &t_vz);
    tree->SetBranchAddress("t_px", &t_px);
    tree->SetBranchAddress("t_py", &t_py);
    tree->SetBranchAddress("t_pz", &t_pz);
    tree->SetBranchAddress("t_charge", &t_charge);
    tree->SetBranchAddress("majorityParticleId", &t_majorityParticleId);

    TFile *oFile = new TFile(oFileName, "recreate");
    double massMin = 0.4;
    double massMax = 2.5;
    const char *titleMass = "#it{M}(K,#pi) (GeV/#it{c}^{2})";
    const char *titlePt = "#it{p}_{T} (GeV/#it{c})";

    TH1D *hInvMass = new TH1D("hInvMass", Form(";%s;Counts", titleMass), (massMax - massMin) * 1000 * 10, massMin, massMax);
    TH1D *hInvMassSiblings = new TH1D("hInvMassSiblings", Form(";%s;Counts", titleMass), (massMax - massMin) * 1000 * 10, massMin, massMax);
    TH2F *hInvMassVsPt = new TH2F("hInvMassVsPt", Form(";%s;%s;Counts", titlePt, titleMass), 100, 0, 10, (massMax - massMin) * 1000 * 10, massMin, massMax);
    TH1D *t_hInvMass = new TH1D("t_hInvMass", Form(";%s;Counts", titleMass), (massMax - massMin) * 1000 * 10, massMin, massMax);

    for (int iEvent = 0; iEvent < tree->GetEntries(); iEvent++) {
        tree->GetEntry(iEvent);

        std::vector<particle> p1 = {};
        std::vector<particle> p2 = {};
        
        size_t nPart = particle_type->size();
        for (size_t iPart = 0; iPart < nPart; iPart++) {
            int pdg = (*particle_type)[iPart];
            if (std::abs(pdg) != pdg1 && std::abs(pdg) != pdg2) continue;
            
            double t_pvx = (*t_vx)[iPart];
            double t_pvy = (*t_vy)[iPart];
            double t_pvz = (*t_vz)[iPart];
            double rProd = pow(t_pvx * t_pvx + t_pvy * t_pvy + t_pvz * t_pvz, 0.5);
            if (!(rProdMin < rProd && rProd < rProdMax)) continue;

            double t_ppx = (*t_px)[iPart];
            double t_ppy = (*t_py)[iPart];
            double t_ppz = (*t_pz)[iPart];
            int t_pcharge = (*t_charge)[iPart];
            uint64_t t_pmajorityParticleId = (*t_majorityParticleId)[iPart];
            
            double pqop = (*qop)[iPart];
            double pphi = (*phi)[iPart];
            double ptheta = (*theta)[iPart];

            double ppx = 1. / std::abs(pqop) * sin(ptheta) * cos(pphi);
            double ppy = 1. / std::abs(pqop) * sin(ptheta) * sin(pphi);
            double ppz = 1. / std::abs(pqop) * cos(ptheta);
            double pm = pdgdb->GetParticle(std::abs(pdg))->Mass();
            
            uint64_t id = ActsFatras::Barcode(t_pmajorityParticleId).particle();

            auto p = ROOT::Math::PxPyPzMVector(ppx, ppy, ppz, pm);
            auto t_p = ROOT::Math::PxPyPzMVector(t_ppx, t_ppy, t_ppz, pm);
            if (std::abs(pdg) ==  pdg1) {
                p1.push_back(particle({pdg, id, p, t_p}));
            } else {
                p2.push_back(particle({pdg, id, p, t_p}));
            }
            
        }

        // compute the invariant mass
        for (const auto &pp1 : p1) {
            for (const auto &pp2 : p2) {
                if (pp1.pdg * pp2.pdg < 0) {
                    // reconstructed
                    hInvMass->Fill((pp1.p + pp2.p).M());
                    if (AreSiblings(pp1, pp2)) {
                        hInvMassSiblings->Fill((pp1.p + pp2.p).M());
                    }
                    hInvMassVsPt->Fill((pp1.p + pp2.p).Pt(), (pp1.p + pp2.p).M());

                    // MC truth
                    t_hInvMass->Fill((pp1.t_p + pp2.t_p).M());
                }
            }
        }
    }
    hInvMass->Write();
    hInvMassSiblings->Write();
    hInvMassVsPt->Write();
    t_hInvMass->Write();
    oFile->Close();
}