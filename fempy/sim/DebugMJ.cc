#include "./functions.hxx"

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

void DebugMJ() {
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
    int iPart = 0;
    (*particles)[iPart++] = new TParticle(211, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.545, 2.423, -0.864, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.461, 2.467, -0.864, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.468, 2.474, -0.776, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.395, 2.492, -0.749, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.378, 2.483, -0.666, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.292, 2.483, -0.636, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.202, 2.472, -0.602, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.169, 2.434, -0.513, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.103, 2.394, -0.445, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.015, 2.350, -0.392, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.004, 2.272, -0.266, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.837, 2.232, -0.272, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.787, 2.146, -0.172, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.649, 2.077, -0.145, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.573, 1.981, -0.061, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.468, 1.885, -0.000, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.453, 1.754, 0.140, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.302, 1.656, 0.163, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.253, 1.521, 0.276, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.119, 1.405, 0.314, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.111, 1.246, 0.461, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.110, 1.145, 0.421, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.120, 0.977, 0.563, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.277, 0.849, 0.574, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.363, 0.697, 0.644, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.472, 0.549, 0.690, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.554, 0.392, 0.755, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.571, 0.215, 0.872, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.737, 0.081, 0.854, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.741, -0.100, 0.971, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.897, -0.234, 0.949, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.992, -0.385, 0.973, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.996, -0.560, 1.070, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.120, -0.696, 1.054, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.135, -0.860, 1.126, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.208, -1.003, 1.137, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.316, -1.129, 1.111, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.370, -1.266, 1.122, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.420, -1.397, 1.128, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.455, -1.526, 1.137, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.529, -1.635, 1.104, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.516, -1.762, 1.136, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.539, -1.869, 1.127, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.560, -1.968, 1.110, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.618, -2.047, 1.052, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.511, -2.164, 1.127, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.586, -2.217, 1.034, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.576, -2.285, 1.006, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.597, -2.332, 0.941, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.578, -2.381, 0.902, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.583, -2.411, 0.833, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.428, -2.477, 0.894, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.422, -2.487, 0.818, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.387, -2.495, 0.759, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.340, -2.495, 0.702, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.336, -2.470, 0.601, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.259, -2.456, 0.556, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.202, -2.425, 0.487, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.200, -2.367, 0.364, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.077, -2.333, 0.341, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -1.001, -2.276, 0.272, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.884, -2.220, 0.235, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.725, -2.167, 0.230, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.721, -2.059, 0.086, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.696, -1.947, -0.043, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.533, -1.869, -0.052, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.362, -1.785, -0.057, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.364, -1.642, -0.214, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.225, -1.534, -0.247, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.150, -1.400, -0.337, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, -0.064, -1.264, -0.416, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.047, -1.131, -0.472, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.149, -0.990, -0.534, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.254, -0.847, -0.591, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.387, -0.709, -0.620, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.468, -0.553, -0.690, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.544, -0.394, -0.761, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.680, -0.252, -0.773, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.739, -0.087, -0.849, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.791, 0.080, -0.924, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 0.929, 0.220, -0.918, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.026, 0.370, -0.941, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.086, 0.529, -0.989, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.089, 0.701, -1.079, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.199, 0.837, -1.067, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.224, 0.994, -1.122, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.286, 1.134, -1.135, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.314, 1.279, -1.170, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.445, 1.386, -1.105, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.461, 1.521, -1.131, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.431, 1.661, -1.188, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.565, 1.744, -1.093, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.548, 1.864, -1.120, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.551, 1.969, -1.119, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.557, 2.063, -1.105, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.588, 2.139, -1.060, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.557, 2.224, -1.060, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.623, 2.269, -0.966, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.552, 2.344, -0.981, 0, 0, 0, 0, 0);
    (*particles)[iPart++] = new TParticle(221, 1, 0, 0, 0, 0, 1.538, 2.392, -0.938, 0, 0, 0, 0, 0);

    TCanvas *cBefore = new TCanvas("cBefore", "Before", 2000, 1000);
    cBefore->Divide(2, 1);
    auto pad = cBefore->cd(1);
    EventDisplay(particles, "before", pad);

    // rotate
    std::vector<ROOT::Math::PxPyPzMVector> parts = GetLorentzVectors(particles);

    parts = Rotate(GetRotationSpheri3DFull(particles), parts);

    TClonesArray *particlesAfter = new TClonesArray("TParticle", 1000);
    for (int iPart = 0; iPart < nPart; iPart++) {
        (*particlesAfter)[iPart] =
            new TParticle(211, 1, 0, 0, 0, 0, parts[iPart].Px(), parts[iPart].Py(), parts[iPart].Pz(), 0, 0, 0, 0, 0);
    }

    print(ComputeSphericityMatrixAll(particlesAfter));

    pad = cBefore->cd(2);
    cBefore->Update();
    EventDisplay(particlesAfter, "after", pad);
}
