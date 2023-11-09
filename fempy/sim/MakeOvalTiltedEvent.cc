#include "MiniJetFactorization.cc"  // NOLINT

void EventDisplay(TClonesArray *particles, const char *name, TVirtualPad *pad) {
    // TCanvas *cEvent = new TCanvas(Form("cEvent_%s", name), "", 1000, 1000);
    // cEvent->cd();
    pad->cd();
    TH3F *hEvent = new TH3F(Form("hEvent_%s", name), ";x;y;z", 10, -3, 3, 10, -3, 3, 10, -3, 3);
    hEvent->DrawClone();

    TPolyLine3D *lAxisX = new TPolyLine3D(2);
    lAxisX->SetPoint(0, 0, 0, 0);
    lAxisX->SetPoint(1, 3, 0, 0);
    lAxisX->SetLineWidth(3);
    lAxisX->DrawClone("same");
    delete lAxisX;
    TPolyLine3D *lAxisY = new TPolyLine3D(2);
    lAxisY->SetPoint(0, 0, 0, 0);
    lAxisY->SetPoint(1, 0, 3, 0);
    lAxisY->SetLineWidth(3);
    lAxisY->DrawClone("same");
    delete lAxisY;
    TPolyLine3D *lAxisZ = new TPolyLine3D(2);
    lAxisZ->SetPoint(0, 0, 0, 0);
    lAxisZ->SetPoint(1, 0, 0, 3);
    lAxisZ->SetLineWidth(3);
    lAxisZ->DrawClone("same");
    delete lAxisZ;

    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        TPolyLine3D *lParticle = new TPolyLine3D(2);
        lParticle->SetPoint(0, 0, 0, 0);
        lParticle->SetPoint(1, particle->Px(), particle->Py(), particle->Pz());
        lParticle->DrawClone("same");
        delete lParticle;
    }

    gStyle->SetOptStat(0);
}

void print(std::array<std::array<double, 3>, 3> mat) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%.3f  ", mat[i][j]);
        }
        printf("\n");
    }
}

void print(std::array<double, 3> vec) {
    for (int i = 0; i < 3; i++) {
        printf("%.3f  ", vec[i]);
    }
    printf("\n");
}

std::vector<ROOT::Math::PxPyPzMVector> GetLorentzVectors(TClonesArray *particles) {
    std::vector<ROOT::Math::PxPyPzMVector> part{};

    for (int iPart = 0; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        part.push_back(ROOT::Math::PxPyPzMVector(particle->Px(), particle->Py(), particle->Pz(), particle->GetMass()));
    }
    return part;
}

void MakeOvalTiltedEvent() {
    // define outputs
    // char* name = (char*)"";
    // char* title = (char*)";#it{k}* (GeV/#it{c});pairs";

    // TH1F *hSE = new TH1F("hSE", title, 3000, 0., 3.);
    // TH1F *hME = new TH1F("hME", title, 3000, 0., 3.);
    // TH1F *hSphi = new TH1F("hSphi", ";Sphericity;Counts", 100, 0., 1);
    // TH1F *hMult = new TH1F("hMult", ";Multiplicity;Counts", 50, 0., 50);
    // TH1F *hPx = new TH1F("hPx", ";Px;Counts", 200, -2., 2);
    // TH1F *hPy = new TH1F("hPy", ";Py;Counts", 200, -2., 2);
    // TH1F *hPz = new TH1F("hPz", ";Pz;Counts", 200, -20., 20);
    // TH1F *hPrincipalAxisTheta = new TH1F("hPrincipalAxisTheta", ";PrincipalAxisTheta;Counts", 200, -1, 1);
    // //__________________________________________________________
    // // perform the simulation
    // std::vector<ROOT::Math::PxPyPzMVector> part1{};
    // std::vector<ROOT::Math::PxPyPzMVector> part2{};
    // std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBuffer1{};
    // std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBuffer2{};

    // TFile outFile("AnalysisResultsDebug.root", "recreate");

    TClonesArray *particles = new TClonesArray("TParticle", 1000);

    int nPart = 100;
    for (int iPart = 0; iPart < nPart; iPart++) {
        double frac = static_cast<double>(iPart) / nPart;
        (*particles)[iPart] = new TParticle(211, 1, 0, 0, 0, 0, 3 * cos(6.28 * frac), sin(6.28 * frac),
                                            gRandom->Gaus() * 0.05, 0, 0, 0, 0, 0);
    }

    TCanvas *cBefore = new TCanvas("cBefore", "Before", 2000, 1000);
    cBefore->Divide(2, 1);
    auto pad = cBefore->cd(1);
    EventDisplay(particles, "before", pad);
    // rotate
    auto parts = GetLorentzVectors(particles);
    parts = Rotate(Mul(RotZ(40 * 3.14 / 180), Mul(RotX(45 * 3.14 / 180), RotY(25 * 3.14 / 180))), parts);

    TClonesArray *particlesAfter = new TClonesArray("TParticle", 1000);
    for (int iPart = 0; iPart < nPart; iPart++) {
        (*particlesAfter)[iPart] =
            new TParticle(211, 1, 0, 0, 0, 0, parts[iPart].Px(), parts[iPart].Py(), parts[iPart].Pz(), 0, 0, 0, 0, 0);
        printf(" %.3f, %.3f, %.3f\n", parts[iPart].Px(), parts[iPart].Py(), parts[iPart].Pz());
    }

    pad = cBefore->cd(2);
    cBefore->Update();
    EventDisplay(particlesAfter, "after", pad);

    print(ComputeSphericityMatrixAll(particles));

    // save root output file
    // hSE->Write();
    // hME->Write();
    // hSphi->Write();
    // hMult->Write();
    // hPx->Write();
    // hPy->Write();
    // hPz->Write();
    // hPrincipalAxisTheta->Write();

    // outFile.Close();
}
