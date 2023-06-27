#include <TClonesArray.h>
#include <TParticle.h>
#include "AliPythia8.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

namespace {
enum tunes { kMonash = 0, kCRMode0, kCRMode2, kCRMode3 };
enum processes { kSoftQCD = 0, kHardQCD };

}

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


struct FemtoParticle {
    ROOT::Math::PxPyPzMVector p;
    int pdg;
    int idx;
    std::pair<int, int> daus;
};


using v2 = std::array<double, 2>;
using v3 = std::array<double, 3>;
using m2 = std::array<std::array<double, 2>, 2>;
using m3 = std::array<std::array<double, 3>, 3>;

void print(m3 mat) {
    for (int i = 0; i<3; i++) {
        for (int j = 0; j<3; j++) {
            printf("%.3f  ", mat[i][j]);
        }
        printf("\n");
    }
}

void print(v3 vec) {
    for (int i = 0; i<3; i++) {
        printf("%.3f  ", vec[i]);
    }
    printf("\n");
}


// #############################################################################
// Math ########################################################################
// #############################################################################

// signum function: return the signum of value, for value = 0 returns 0
int sgn(double value) { return (value > 0) - (value < 0); }


// 
template <typename Pair, typename Value>
bool In(Pair selection, Value value) {
    return selection.first <= value && value < selection.second;
}
// #############################################################################
// Linear algebra ##############################################################
// #############################################################################

// Returns a 3d array of zeros
v3 Zero() { return {0, 0, 0}; }


m3 One() {
    return m3({
        v3({1, 0, 0}),
        v3({0, 1, 0}),
        v3({0, 0, 1})
    });
}


// Returns the generator of the rotations along the X axis
m3 RotX(double angle) {
    double c = cos(angle);
    double s = sin(angle);

    m3 rot;
    rot[0] = {1, 0, 0};
    rot[1] = {0, c, -s};
    rot[2] = {0, s, c};

    return rot;
}

// Returns the generator of the rotations along the Y axis
m3 RotY(double angle) {
    double c = cos(angle);
    double s = sin(angle);

    m3 rot;
    rot[0] = {c, 0, s};
    rot[1] = {0, 1, 0};
    rot[2] = {-s, 0, c};

    return rot;
}

// Returns the generator of the rotations along the Z axis
m3 RotZ(double angle) {
    double c = cos(angle);
    double s = sin(angle);

    m3 rot;
    rot[0] = {c, -s, 0};
    rot[1] = {s, c, 0};
    rot[2] = {0, 0, 1};

    return rot;
}

// Multiply a 3x3 matrix and a 3d array
v3 Mul(m3 left, v3 right) {
    v3 ans = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ans[i] += left[i][j] * right[j];
        }
    }

    return ans;
}

// Multiply two 3x3 matrix
m3 Mul(m3 left, m3 right) {
    m3 ans{Zero(), Zero(), Zero()};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                ans[i][j] += left[i][k] * right[k][j];
            }
        }
    }

    return ans;
}

// Compute the eigenvector of a 3x3 matrix given the eigenvalue
v3 ComputeEigenVector(m3 mat, double lambda) {
    double b = (lambda * mat[1][2] + mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) /
               (lambda * mat[0][2] + mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);
    double c = (lambda - mat[0][0] - mat[0][1] * b) / mat[0][2];

    v3 eivec = {1, b, c};
    return eivec;
}

// Compute the eigenvector of a 3x3 matrix given the eigenvalue
v3 ComputeEigenVector(double mat[3][3], double lambda) {
    m3 matVec = {v3({mat[0][0], mat[0][1], mat[0][2]}), v3({mat[1][0], mat[1][1], mat[1][2]}),
                 v3({mat[2][0], mat[2][1], mat[2][2]})};
    return ComputeEigenVector(matVec, lambda);
}

// Diagonalize a 3x3 symmetric matrix and returns an array with its eigenvalues, serted in increasing order
v3 ComputeEigenValues(m3 mat) {
    double b = mat[0][0] + mat[1][1] + mat[2][2];
    double c = mat[0][0] * mat[1][1] + mat[0][0] * mat[2][2] + mat[1][1] * mat[2][2] - mat[0][1] * mat[0][1] -
               mat[0][2] * mat[0][2] - mat[1][2] * mat[1][2];
    double d = mat[0][0] * mat[1][2] * mat[1][2] + mat[1][1] * mat[0][2] * mat[0][2] +
               mat[2][2] * mat[0][1] * mat[0][1] - mat[0][0] * mat[1][1] * mat[2][2] -
               2 * mat[0][1] * mat[0][2] * mat[1][2];
    double p = b * b - 3 * c;
    double q = 2 * b * b * b - 9 * b * c - 27 * d;
    double delta = acos(q / 2 / pow(p, 1.5));

    // compute the eigenvalues
    v3 eivals;
    eivals[0] = 1. / 3 * (b + 2 * pow(p, 0.5) * cos(delta / 3));
    eivals[1] = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta + 2 * TMath::Pi()) / 3));
    eivals[2] = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta - 2 * TMath::Pi()) / 3));
    std::sort(eivals.begin(), eivals.end());
    return eivals;
}

double Determinant(m2 mat) {
    return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
}

double Determinant(m3 mat) {
    m2 smi = m2({v2({mat[1][1], mat[1][2]}),
                 v2({mat[2][1], mat[2][2]})});
    
    m2 smj = m2({v2({mat[1][0], mat[1][2]}),
                 v2({mat[2][0], mat[2][2]})});
    
    m2 smk = m2({v2({mat[1][0], mat[1][1]}),
                 v2({mat[2][0], mat[2][1]})});
    
    return mat[0][0] * Determinant(smi) - mat[0][1] * Determinant(smj) + mat[0][2] * Determinant(smk);
}


template <typename T>
bool IsZero(T mat) {
    for (const auto & row : mat)
        for (const auto & val : row)
            if (val != 0) return false;
    return true;
}

// Apply rotation to a list of lorentz vector in 2 spatial dimensions
std::vector<ROOT::Math::PxPyPzMVector> RotateXY(std::array<std::array<double, 2>, 2> rot,
                                                std::vector<ROOT::Math::PxPyPzMVector> particles) {
    exit(1);
    std::vector<ROOT::Math::PxPyPzMVector> rotated = {};
    for (const auto &particle : particles) {
        rotated.push_back(ROOT::Math::PxPyPzMVector(rot[0][0] * particle.Px() + rot[0][1] * particle.Py(),
                                                    rot[1][0] * particle.Px() + rot[1][1] * particle.Py(),
                                                    particle.Pz(), particle.M()));
    }
    return rotated;
}

// Apply rotation to a list of lorentz vector in 3 spatial dimensions
std::vector<ROOT::Math::PxPyPzMVector> Rotate(m3 rot, std::vector<ROOT::Math::PxPyPzMVector> particles) {
    std::vector<ROOT::Math::PxPyPzMVector> rotated = {};
    for (const auto &particle : particles) {
        rotated.push_back(ROOT::Math::PxPyPzMVector(
            rot[0][0] * particle.Px() + rot[0][1] * particle.Py() + rot[0][2] * particle.Pz(),
            rot[1][0] * particle.Px() + rot[1][1] * particle.Py() + rot[1][2] * particle.Pz(),
            rot[2][0] * particle.Px() + rot[2][1] * particle.Py() + rot[2][2] * particle.Pz(), particle.M()));
    }
    return rotated;
}

// Apply rotation to a list of lorentz vector in 3 spatial dimensions
std::vector<FemtoParticle> Rotate(m3 rot, std::vector<FemtoParticle> particles) {
    std::vector<FemtoParticle> rotated = {};
    for (const auto &particle : particles) {
        ROOT::Math::PxPyPzMVector pRot(
            rot[0][0] * particle.p.Px() + rot[0][1] * particle.p.Py() + rot[0][2] * particle.p.Pz(),
            rot[1][0] * particle.p.Px() + rot[1][1] * particle.p.Py() + rot[1][2] * particle.p.Pz(),
            rot[2][0] * particle.p.Px() + rot[2][1] * particle.p.Py() + rot[2][2] * particle.p.Pz(), particle.p.M());
        FemtoParticle particleRot = particle;
        particleRot.p = pRot;
        rotated.push_back(particleRot);
    }
    return rotated;
}


// Rotate the momentum vector of particle
FemtoParticle Rotate(m3 rot, const FemtoParticle &particle) {
    ROOT::Math::PxPyPzMVector pRot(
        rot[0][0] * particle.p.Px() + rot[0][1] * particle.p.Py() + rot[0][2] * particle.p.Pz(),
        rot[1][0] * particle.p.Px() + rot[1][1] * particle.p.Py() + rot[1][2] * particle.p.Pz(),
        rot[2][0] * particle.p.Px() + rot[2][1] * particle.p.Py() + rot[2][2] * particle.p.Pz(),
        particle.p.M());
    FemtoParticle particleRot = particle;
    particleRot.p = pRot;
    return particleRot;
}


// Rotate the momentum vector of particle
ROOT::Math::PxPyPzMVector Rotate(m3 rot, const ROOT::Math::PxPyPzMVector &particle) {
    ROOT::Math::PxPyPzMVector pRot(
        rot[0][0] * particle.Px() + rot[0][1] * particle.Py() + rot[0][2] * particle.Pz(),
        rot[1][0] * particle.Px() + rot[1][1] * particle.Py() + rot[1][2] * particle.Pz(),
        rot[2][0] * particle.Px() + rot[2][1] * particle.Py() + rot[2][2] * particle.Pz(),
        particle.M());
    return pRot;
}


// #############################################################################
// Physics #####################################################################
// #############################################################################

// return true if the particle is in the TPC acceptance
inline bool isInTPC(TParticle *p) { return p->GetStatusCode() == 1 && p->Pt() > 0.3 && std::abs(p->Eta()) < 0.8; }

// Return true if the particle is charged and long lived and in the final state
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

// Get the multiplicity of charged particles in the V0 acceptance (same as ALICE high mult trigger)
int GetMultV0(const TClonesArray *particles) {
    int nCh = 0;
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
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

// Compute the multiplicity of charged particles in the TPC acceptance
int GetMultTPC(TClonesArray *particles) {
    // evaluate multiplicity at forward rapidity
    int nParticles = 0;
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        int pdg = std::abs(particle->GetPdgCode());
        int status = std::abs(particle->GetStatusCode());
        float eta = particle->Eta();

        if (IsDetectable(pdg) && std::abs(eta) < 0.8 && status == 1) {  // V0A and V0C acceptance
            nParticles++;
        }
    }
    return nParticles;
}

// Compute the relative momentum in the center of mass reference frame
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2) {
    ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
    ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
    ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);

    ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;
    float kStar = 0.5 * trackRelK.P();
    return kStar;
}

// Compute the relative momentum in the center of mass reference frame
float ComputeKstarFD(ROOT::Math::PxPyPzMVector PartOne, ROOT::Math::PxPyPzMVector PartTwo) {
    ROOT::Math::PxPyPzMVector trackSum = PartOne + PartTwo;

    float beta = trackSum.Beta();
    float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
    float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
    float betaz = beta * cos(trackSum.Theta());

    TLorentzVector PartOneCMS = TLorentzVector(PartOne.Px(), PartOne.Py(), PartOne.Pz(), PartOne.E());
    TLorentzVector PartTwoCMS = TLorentzVector(PartTwo.Px(), PartTwo.Py(), PartTwo.Pz(), PartTwo.E());

    PartOneCMS.Boost(-betax, -betay, -betaz);
    PartTwoCMS.Boost(-betax, -betay, -betaz);

    TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

    return 0.5 * trackRelK.P();
}

// Compute the 3x3 sphericity matrix with the particles in the TPC
m3 ComputeSphericityMatrix(TClonesArray *particles) {
    // Compute the sphericity matrix
    m3 sphi = {(v3){0, 0, 0}, (v3){0, 0, 0}, (v3){0, 0, 0}};
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        if (isInTPC(particle) && IsDetectable(std::abs(particle->GetPdgCode()))) {
            sphi[0][0] += particle->Px() * particle->Px();
            sphi[0][1] += particle->Px() * particle->Py();
            sphi[0][2] += particle->Px() * particle->Pz();
            sphi[1][0] += particle->Py() * particle->Px();
            sphi[1][1] += particle->Py() * particle->Py();
            sphi[1][2] += particle->Py() * particle->Pz();
            sphi[2][0] += particle->Pz() * particle->Px();
            sphi[2][1] += particle->Pz() * particle->Py();
            sphi[2][2] += particle->Pz() * particle->Pz();
        }
    }
    return sphi;
}

// Compute the 3x3 sphericity matrix with all the particles in the event
m3 ComputeSphericityMatrixAll(TClonesArray *particles) {
    // Compute the sphericity matrix
    m3 sphi = {(v3){0, 0, 0}, (v3){0, 0, 0}, (v3){0, 0, 0}};
    for (int iPart = 0; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        sphi[0][0] += particle->Px() * particle->Px();
        sphi[0][1] += particle->Px() * particle->Py();
        sphi[0][2] += particle->Px() * particle->Pz();
        sphi[1][0] += particle->Py() * particle->Px();
        sphi[1][1] += particle->Py() * particle->Py();
        sphi[1][2] += particle->Py() * particle->Pz();
        sphi[2][0] += particle->Pz() * particle->Px();
        sphi[2][1] += particle->Pz() * particle->Py();
        sphi[2][2] += particle->Pz() * particle->Pz();
    }
    return sphi;
}

// Compute the sphericity using the minimum eigenvalue of the sphericity matrix
double ComputeMajorSphericity(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    return 3 * eivals[0] / (eivals[0] + eivals[1] + eivals[2]);
}

// Compute the sphericity for the biggest and mid direction
double ComputeMajorRelativeSphericity(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    return 2 * eivals[1] / (eivals[1] + eivals[2]);
}

// Compute the sphericity using the intermediate eigenvalue of the sphericity matrix
double ComputeMinorSphericity(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    return 2 * eivals[1] / (eivals[0] + eivals[1] + eivals[2]);
}

// Compute the sphericity for the mid and smallest direction
double ComputeMinorRelativeSphericity(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    return 2 * eivals[0] / (eivals[0] + eivals[1]);
}

// Compute transverse sphericity
double ComputeTransverseSphericity(TClonesArray *particles) {
    // Compute the sphericity matrix
    double sphi[2][2] = {{0, 0}, {0, 0}};
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        sphi[0][0] += particle->Px() * particle->Px();
        sphi[0][1] += particle->Px() * particle->Py();
        sphi[1][0] += particle->Py() * particle->Px();
        sphi[1][1] += particle->Py() * particle->Py();
    }

    // Compute the eigenvalues
    double ecc = pow((sphi[0][0] - sphi[1][1]) * (sphi[0][0] - sphi[1][1]) + 4 * sphi[0][1] * sphi[0][1], 0.5);
    double eval1 = 0.5 * (sphi[0][0] + sphi[1][1] + sgn(sphi[0][0] - sphi[1][1]) * ecc);
    double eval2 = 0.5 * (sphi[0][0] + sphi[1][1] - sgn(sphi[0][0] - sphi[1][1]) * ecc);

    return 2 * std::min(eval1, eval2) / (eval1 + eval2);
}

// Compute the rotation matrix for the event based on the sphericity (no minor axis rotation)
m3 GetRotationSpheri3D(TClonesArray *particles) {
    double sphi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        if (isInTPC(particle)) {
            sphi[0][0] += particle->Px() * particle->Px();
            sphi[0][1] += particle->Px() * particle->Py();
            sphi[0][2] += particle->Px() * particle->Pz();
            sphi[1][0] += particle->Py() * particle->Px();
            sphi[1][1] += particle->Py() * particle->Py();
            sphi[1][2] += particle->Py() * particle->Pz();
            sphi[2][0] += particle->Pz() * particle->Px();
            sphi[2][1] += particle->Pz() * particle->Py();
            sphi[2][2] += particle->Pz() * particle->Pz();
        }
    }

    // printf("sphericity matrix:\n\n");
    // for (int i = 0; i<3; i++){
    //     for (int j = 0; j<3; j++){
    //         printf("%.3f  ", sphi[i][j]);
    //     }
    //     printf("\n");
    // }
    // auxiliary variables to solve the cubic characteristic eq of the matrix
    double b = sphi[0][0] + sphi[1][1] + sphi[2][2];
    double c = sphi[0][0] * sphi[1][1] + sphi[0][0] * sphi[2][2] + sphi[1][1] * sphi[2][2] - sphi[0][1] * sphi[0][1] -
               sphi[0][2] * sphi[0][2] - sphi[1][2] * sphi[1][2];
    double d = sphi[0][0] * sphi[1][2] * sphi[1][2] + sphi[1][1] * sphi[0][2] * sphi[0][2] +
               sphi[2][2] * sphi[0][1] * sphi[0][1] - sphi[0][0] * sphi[1][1] * sphi[2][2] -
               2 * sphi[0][1] * sphi[0][2] * sphi[1][2];
    double p = b * b - 3 * c;
    double q = 2 * b * b * b - 9 * b * c - 27 * d;
    double delta = acos(q / 2 / pow(p, 1.5));

    // compute the eigenvalues
    double l1 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos(delta / 3));
    double l2 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta + 2 * TMath::Pi()) / 3));
    double l3 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta - 2 * TMath::Pi()) / 3));

    // printf("eigenvectors: %.3f  %.3f  %.3f\n", l1, l2, l3);

    // compute the largest eigenvector -> principal axis of the event
    double lMax = std::max(std::max(l1, l2), l3);

    double evecMax[3];
    evecMax[0] = 1;
    evecMax[1] = (sphi[1][0] + (sphi[1][2] * (lMax - sphi[0][0])) / (sphi[0][2])) /
                 (lMax + (sphi[1][2] * sphi[0][1]) / (sphi[0][2]) - sphi[1][1]);
    evecMax[2] = (lMax - sphi[0][0] - sphi[0][1] * evecMax[1]) / (sphi[0][2]);

    float phi = atan(evecMax[1]);
    float theta = atan(pow(evecMax[0] * evecMax[0] + evecMax[1] * evecMax[1], 0.5) / evecMax[2]);

    float ct = cos(theta);
    float st = sin(theta);
    float cp = cos(phi);
    float sp = sin(phi);

    m3 rot;
    rot[0] = {ct * cp, ct * sp, -st};
    rot[1] = {-sp, cp, 0};
    rot[2] = {st * cp, st * sp, ct};

    return rot;
}

// Compute the rotation matrix for the event based on the sphericity (with minor axis included)
m3 GetRotationSpheri3DFull(TClonesArray *particles) {
    double sphi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        if (isInTPC(particle) && IsDetectable(std::abs(particle->GetPdgCode()))) {
            sphi[0][0] += particle->Px() * particle->Px();
            sphi[0][1] += particle->Px() * particle->Py();
            sphi[0][2] += particle->Px() * particle->Pz();
            sphi[1][0] += particle->Py() * particle->Px();
            sphi[1][1] += particle->Py() * particle->Py();
            sphi[1][2] += particle->Py() * particle->Pz();
            sphi[2][0] += particle->Pz() * particle->Px();
            sphi[2][1] += particle->Pz() * particle->Py();
            sphi[2][2] += particle->Pz() * particle->Pz();
        }
    }

    // printf("sphericity matrix:\n\n");
    // for (int i = 0; i<3; i++){
    //     for (int j = 0; j<3; j++){
    //         printf("%.3f  ", sphi[i][j]);
    //     }
    //     printf("\n");
    // }
    // auxiliary variables to solve the cubic characteristic eq of the matrix
    double b = sphi[0][0] + sphi[1][1] + sphi[2][2];
    double c = sphi[0][0] * sphi[1][1] + sphi[0][0] * sphi[2][2] + sphi[1][1] * sphi[2][2] - sphi[0][1] * sphi[0][1] -
               sphi[0][2] * sphi[0][2] - sphi[1][2] * sphi[1][2];
    double d = sphi[0][0] * sphi[1][2] * sphi[1][2] + sphi[1][1] * sphi[0][2] * sphi[0][2] +
               sphi[2][2] * sphi[0][1] * sphi[0][1] - sphi[0][0] * sphi[1][1] * sphi[2][2] -
               2 * sphi[0][1] * sphi[0][2] * sphi[1][2];
    double p = b * b - 3 * c;
    double q = 2 * b * b * b - 9 * b * c - 27 * d;
    double delta = acos(q / 2 / pow(p, 1.5));

    // compute the eigenvalues
    double l1 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos(delta / 3));
    double l2 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta + 2 * TMath::Pi()) / 3));
    double l3 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta - 2 * TMath::Pi()) / 3));

    // printf("eigenvectors: %.3f  %.3f  %.3f\n", l1, l2, l3);

    // compute the largest eigenvector -> principal axis of the event
    // sort eigenvalues
    std::vector<double> eigenvalues = {l1, l2, l3};
    std::sort(eigenvalues.begin(), eigenvalues.end());
    double lMed = eigenvalues[1];
    double lMax = eigenvalues[2];

    double evecMax[3];
    evecMax[0] = 1;
    evecMax[1] = (sphi[1][0] + (sphi[1][2] * (lMax - sphi[0][0])) / (sphi[0][2])) /
                 (lMax + (sphi[1][2] * sphi[0][1]) / (sphi[0][2]) - sphi[1][1]);
    evecMax[2] = (lMax - sphi[0][0] - sphi[0][1] * evecMax[1]) / (sphi[0][2]);

    float phi = atan(evecMax[1]);
    float theta = atan(pow(evecMax[0] * evecMax[0] + evecMax[1] * evecMax[1], 0.5) / evecMax[2]);

    // float ct = cos(theta);
    // float st = sin(theta);
    // float cp = cos(phi);
    // float sp = sin(phi);

    m3 rot = Mul(RotY(-theta), RotZ(-phi));
    // rot[0] = { ct*cp, ct*sp, -st};
    // rot[1] = { -sp  , cp   , 0  };
    // rot[2] = { st*cp, st*sp, ct };

    v3 evecMed = ComputeEigenVector(sphi, lMed);

    v3 evecMedRot = Mul(rot, evecMed);
    // evecMedRot[0] = rot[0][0] * evecMed[0] + rot[0][1] * evecMed[1] + rot[0][2] * evecMed[2];
    // evecMedRot[1] = rot[1][0] * evecMed[0] + rot[1][1] * evecMed[1] + rot[1][2] * evecMed[2];
    // evecMedRot[2] = rot[2][0] * evecMed[0] + rot[2][1] * evecMed[1] + rot[2][2] * evecMed[2];

    double phiSecondaryAxis = atan(evecMedRot[1] / evecMedRot[0]);
    m3 rotSecondaryAxis = RotZ(-phiSecondaryAxis);

    return Mul(rotSecondaryAxis, rot);
}


// Compute the rotation matrix for the event based on the sphericity (with minor axis included)
m3 GetRotationSpheri3DFull(m3 sphi) {
    double b = sphi[0][0] + sphi[1][1] + sphi[2][2];
    double c = sphi[0][0] * sphi[1][1] + sphi[0][0] * sphi[2][2] + sphi[1][1] * sphi[2][2] - sphi[0][1] * sphi[0][1] -
               sphi[0][2] * sphi[0][2] - sphi[1][2] * sphi[1][2];
    double d = sphi[0][0] * sphi[1][2] * sphi[1][2] + sphi[1][1] * sphi[0][2] * sphi[0][2] +
               sphi[2][2] * sphi[0][1] * sphi[0][1] - sphi[0][0] * sphi[1][1] * sphi[2][2] -
               2 * sphi[0][1] * sphi[0][2] * sphi[1][2];
    double p = b * b - 3 * c;
    double q = 2 * b * b * b - 9 * b * c - 27 * d;
    double delta = acos(q / 2 / pow(p, 1.5));

    // Compute the eigenvalues
    double l1 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos(delta / 3));
    double l2 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta + 2 * TMath::Pi()) / 3));
    double l3 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta - 2 * TMath::Pi()) / 3));

    // Sort the eigenvalues
    std::vector<double> eigenvalues = {l1, l2, l3};
    std::sort(eigenvalues.begin(), eigenvalues.end());
    double lMed = eigenvalues[1];
    double lMax = eigenvalues[2];

    double evecMax[3];
    evecMax[0] = 1;
    evecMax[1] = (sphi[1][0] + (sphi[1][2] * (lMax - sphi[0][0])) / (sphi[0][2])) /
                 (lMax + (sphi[1][2] * sphi[0][1]) / (sphi[0][2]) - sphi[1][1]);
    evecMax[2] = (lMax - sphi[0][0] - sphi[0][1] * evecMax[1]) / (sphi[0][2]);

    float phi = atan(evecMax[1]);
    float theta = atan(pow(evecMax[0] * evecMax[0] + evecMax[1] * evecMax[1], 0.5) / evecMax[2]);

    m3 rot = Mul(RotY(-theta), RotZ(-phi));
    v3 evecMed = ComputeEigenVector(sphi, lMed);
    v3 evecMedRot = Mul(rot, evecMed);

    double phiSecondaryAxis = atan(evecMedRot[1] / evecMedRot[0]);
    m3 rotSecondaryAxis = RotZ(-phiSecondaryAxis);

    return Mul(rotSecondaryAxis, rot);
}


// Get the rotation 2x2 matrix in the trasnverse plane using the siagonalized trasnverse sphericity
std::array<std::array<double, 2>, 2> GetRotationTransverseSpheri(TClonesArray *particles) {
    // Compute the sphericity matrix
    double sphi[2][2] = {{0, 0}, {0, 0}};
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        if (particle->GetStatusCode() == 1) {
            sphi[0][0] += particle->Px() * particle->Px();
            sphi[0][1] += particle->Px() * particle->Py();
            sphi[1][0] += particle->Py() * particle->Px();
            sphi[1][1] += particle->Py() * particle->Py();
        }
    }

    // Compute the eigenvalues
    double ecc = pow((sphi[0][0] - sphi[1][1]) * (sphi[0][0] - sphi[1][1]) + 4 * sphi[0][1] * sphi[0][1], 0.5);
    double eval1 = 0.5 * (sphi[0][0] + sphi[1][1] + sgn(sphi[0][0] - sphi[1][1]) * ecc);
    double eval2 = 0.5 * (sphi[0][0] + sphi[1][1] - sgn(sphi[0][0] - sphi[1][1]) * ecc);

    // Compute the rotation angle in the azimutal direction
    double phi = atan((std::max(eval1, eval2) - sphi[0][0]) / sphi[0][1]);

    // printf("phiii  %.2f\n", phi * 180 / 3.14);
    // Compute the inverse rotation matrix
    std::array<std::array<double, 2>, 2> rot;
    rot[0] = {cos(phi), sin(phi)};
    rot[1] = {-sin(phi), cos(phi)};

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            // printf("%.3f  ", rot[i][j]);
        }
        // printf("\n");
    }
    return rot;
}

// Compute the rotation matrix using the leading track
std::array<std::array<double, 2>, 2> GetRotationXYLeadingTrack(TClonesArray *particles) {
    // Compute the sphericity matrix
    double leadingPhi = 0;
    double leadingPxy = 0;
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));

        if (particle->GetStatusCode() == 1) {
            if (pow(particle->Px() * particle->Px() + particle->Py() * particle->Py(), 0.5) > leadingPxy) {
                leadingPhi = atan(particle->Py() / particle->Px());
            }
        }
    }

    // Compute the inverse rotation matrix
    std::array<std::array<double, 2>, 2> rot;
    rot[0] = {cos(leadingPhi), -sin(leadingPhi)};
    rot[1] = {sin(leadingPhi), cos(leadingPhi)};

    return rot;
}

// Compute the eigenvectors using the 3D sphericity matrix
m3 ComputeEigenVectorsSpheri3D(TClonesArray *particles) {
    double sphi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        if (isInTPC(particle)) {
            sphi[0][0] += particle->Px() * particle->Px();
            sphi[0][1] += particle->Px() * particle->Py();
            sphi[0][2] += particle->Px() * particle->Pz();
            sphi[1][0] += particle->Py() * particle->Px();
            sphi[1][1] += particle->Py() * particle->Py();
            sphi[1][2] += particle->Py() * particle->Pz();
            sphi[2][0] += particle->Pz() * particle->Px();
            sphi[2][1] += particle->Pz() * particle->Py();
            sphi[2][2] += particle->Pz() * particle->Pz();
        }
    }

    // printf("sphericity matrix:\n\n");
    // for (int i = 0; i<3; i++){
    //     for (int j = 0; j<3; j++){
    //         printf("%.3f  ", sphi[i][j]);
    //     }
    //     printf("\n");
    // }
    // auxiliary variables to solve the cubic characteristic eq of the matrix
    double b = sphi[0][0] + sphi[1][1] + sphi[2][2];
    double c = sphi[0][0] * sphi[1][1] + sphi[0][0] * sphi[2][2] + sphi[1][1] * sphi[2][2] - sphi[0][1] * sphi[0][1] -
               sphi[0][2] * sphi[0][2] - sphi[1][2] * sphi[1][2];
    double d = sphi[0][0] * sphi[1][2] * sphi[1][2] + sphi[1][1] * sphi[0][2] * sphi[0][2] +
               sphi[2][2] * sphi[0][1] * sphi[0][1] - sphi[0][0] * sphi[1][1] * sphi[2][2] -
               2 * sphi[0][1] * sphi[0][2] * sphi[1][2];
    double p = b * b - 3 * c;
    double q = 2 * b * b * b - 9 * b * c - 27 * d;
    double delta = acos(q / 2 / pow(p, 1.5));

    // compute the eigenvalues
    double l1 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos(delta / 3));
    double l2 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta + 2 * TMath::Pi()) / 3));
    double l3 = 1. / 3 * (b + 2 * pow(p, 0.5) * cos((delta - 2 * TMath::Pi()) / 3));

    // printf("eigenvectors: %.3f  %.3f  %.3f\n", l1, l2, l3);

    // compute the largest eigenvector -> principal axis of the event
    std::vector<double> eivals = {l1, l2, l3};
    std::sort(eivals.begin(), eivals.end());
    double lMin = eivals[0];
    double lMed = eivals[1];
    double lMax = eivals[2];

    v3 evecMin = ComputeEigenVector(sphi, lMin);
    v3 evecMed = ComputeEigenVector(sphi, lMed);
    v3 evecMax = ComputeEigenVector(sphi, lMax);

    m3 vec = {evecMax, evecMed, evecMin};

    // printf("Eigen vec:\n");
    // printf("lMin (%.3f): %.3f  %.3f  %.3f  \n", eivals[0], evecMin[0], evecMin[1], evecMin[2]);
    // printf("lMed (%.3f): %.3f  %.3f  %.3f  \n", eivals[1], evecMed[0], evecMed[1], evecMed[2]);
    // printf("lMax (%.3f): %.3f  %.3f  %.3f  \n", eivals[2], evecMax[0], evecMax[1], evecMax[2]);

    return vec;
}

// Compute the direction along which the event develops
double ComputePrincipalAxisTheta(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    v3 principalAxis = ComputeEigenVector(sphi, std::max(std::max(eivals[0], eivals[1]), eivals[2]));
    double theta =
        atan(principalAxis[2] / pow(principalAxis[0] * principalAxis[0] + principalAxis[1] * principalAxis[1], 0.5));
    return theta;
}

bool ThereIsD(TClonesArray *particles) {
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++) {
        TParticle *particle = dynamic_cast<TParticle *>(particles->At(iPart));
        int pdg = particle->GetPdgCode();
        // if (std::abs(particle->GetPdgCode()) == 411) return true;

        if (std::abs(pdg) == 413) {  // select decay kinem of Dplus
            return true;
            std::multiset<int> dauPdgs = {};
            for (int iDau = particle->GetFirstDaughter(); iDau <= particle->GetLastDaughter(); iDau++) {
                TParticle *dau = dynamic_cast<TParticle *>(particles->At(iDau));

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
            return true;
        }

    nextparticle:
        continue;
    }
    return false;
}

// #############################################################################
// Pythia ######################################################################
// #############################################################################

// Set the pythia tune
void SetTune(AliPythia8 &pythia, tunes tune) {
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
}

// Set the pythia process
void SetProcess(AliPythia8 &pythia, processes process) {
    if (process == kSoftQCD) {
        pythia.ReadString("SoftQCD:nonDiffractive = on");
    } else if (process == kHardQCD) {
        pythia.ReadString("HardQCD:hardccbar = on");
        pythia.ReadString("HardQCD:hardbbbar = on");
    }
}

inline double ProductionRadius(TParticle * p) {
    return std::sqrt(p->Vx() * p->Vx() + p->Vy() * p->Vy() + p->Vz() * p->Vz());
}

bool IsPairClean(FemtoParticle charm, FemtoParticle light) {
    for (int iPart = charm.daus.first; iPart <= charm.daus.second; iPart++) {
        if (iPart == light.idx) return false;
    }
    return true;
}

std::vector<ROOT::Math::PxPyPzMVector> GetLorentzVectors(TClonesArray *particles) {
    std::vector<ROOT::Math::PxPyPzMVector> part{};

    for (int iPart = 0; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        part.push_back(ROOT::Math::PxPyPzMVector(particle->Px(), particle->Py(), particle->Pz(), particle->GetMass()));
    }
    return part;
}