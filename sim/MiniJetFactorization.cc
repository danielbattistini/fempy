#include <array>
#include <string>
#include <vector>
#include <map>
#include <deque>
// #include <cmath>

#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH3F.h>
#include <TLorentzVector.h>

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"

#include <TClonesArray.h>
#include <TParticle.h>
#include "AliPythia8.h"

int sgn(double value) { return (value > 0) - (value < 0); }

using namespace Pythia8;

namespace 
{
    const std::array<int, 2> lightPDG{211, 321}; // pions, kaons, protons
    const std::array<int, 2> charmPDG{411, 413}; // D+, D*
}

enum FactorizationMode {
    kNo = 0,
    kXY,
    kLeadingXY,
    kLeading,
    kSpheri3D,
    kSpheri3DFull,
};

enum ComputeKStarMode {
    kCustom = 0,
    kFD,
};

// Linear algebra ##############################################################
// Returns a 3d array of zeros
std::array<double, 3> Zero() {
    return {0, 0, 0};
}
// Generator of the rotations along the X axis
std::array<std::array<double, 3>, 3> RotX(double angle) {
    double c = cos(angle);
    double s = sin(angle);

    std::array<std::array<double, 3>, 3> rot;
    rot[0] = {1, 0,  0};
    rot[1] = {0, c, -s}; 
    rot[2] = {0, s,  c};

    return rot;
}

// Generator of the rotations along the Y axis
std::array<std::array<double, 3>, 3> RotY(double angle) {
    double c = cos(angle);
    double s = sin(angle);

    std::array<std::array<double, 3>, 3> rot;
    rot[0] = { c, 0, s};
    rot[1] = { 0, 1, 0};
    rot[2] = {-s, 0, c};

    return rot;
}

// Generator of the rotations along the Z axis
std::array<std::array<double, 3>, 3> RotZ(double angle) {
    double c = cos(angle);
    double s = sin(angle);

    std::array<std::array<double, 3>, 3> rot;
    rot[0] = {c, -s, 0};
    rot[1] = {s,  c, 0};
    rot[2] = {0,  0, 1};

    return rot;
}

// Multiply a 3x3 matrix and a 3d array
std::array<double, 3> Mul(std::array<std::array<double, 3>, 3> left, std::array<double, 3> right) {
    std::array<double, 3> ans = {0, 0, 0};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            ans[i] += left[i][j]*right[j];
        }
    }

    return ans;
}

// Multiply two 3x3 matrix
std::array<std::array<double, 3>, 3> Mul(std::array<std::array<double, 3>, 3> left, std::array<std::array<double, 3>, 3> right) {
    std::array<std::array<double, 3>, 3> ans {Zero(), Zero(), Zero()};

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                ans[i][j] += left[i][k]*right[k][j];
            }
        }
    }

    return ans;
}

// Compute the eigenvector of a 3x3 matrix given the eigenvalue
std::array<double, 3> ComputeEigenVector(double mat[3][3], double lambda) {
    double b = (lambda * mat[1][2] + mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2])/(lambda * mat[0][2] + mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1]);
    double c = (lambda - mat[0][0] - mat[0][1]*b)/mat[0][2];

    std::array<double, 3> eivec = {1, b, c};
    return eivec;
}

// Compute the eigenvector of a 3x3 matrix given the eigenvalue
std::array<double, 3> ComputeEigenVector(std::array<std::array<double, 3>, 3> mat, double lambda) {
    double b = (lambda * mat[1][2] + mat[0][2]*mat[1][0] - mat[0][0]*mat[1][2])/(lambda * mat[0][2] + mat[0][1]*mat[1][2] - mat[0][2]*mat[1][1]);
    double c = (lambda - mat[0][0] - mat[0][1]*b)/mat[0][2];

    std::array<double, 3> eivec = {1, b, c};
    return eivec;
}

std::vector<ROOT::Math::PxPyPzMVector> RotateXY(std::array<std::array<double, 2>, 2> rot, std::vector<ROOT::Math::PxPyPzMVector> particles){
    exit(1);
    std::vector<ROOT::Math::PxPyPzMVector> rotated = {};
    for (const auto &particle : particles) {
        rotated.push_back(ROOT::Math::PxPyPzMVector(
            rot[0][0] * particle.Px() + rot[0][1] * particle.Py(),
            rot[1][0] * particle.Px() + rot[1][1] * particle.Py(),
            particle.Pz(),
            particle.M()
        ));
    }
    return rotated;
}

std::vector<ROOT::Math::PxPyPzMVector> Rotate(std::array<std::array<double, 3>, 3> rot, std::vector<ROOT::Math::PxPyPzMVector> particles){
    std::vector<ROOT::Math::PxPyPzMVector> rotated = {};
    for (const auto &particle : particles) {
        rotated.push_back(ROOT::Math::PxPyPzMVector(
            rot[0][0] * particle.Px() + rot[0][1] * particle.Py() + rot[0][2] * particle.Pz(),
            rot[1][0] * particle.Px() + rot[1][1] * particle.Py() + rot[1][2] * particle.Pz(),
            rot[2][0] * particle.Px() + rot[2][1] * particle.Py() + rot[2][2] * particle.Pz(),
            particle.M()
        ));
    }
    return rotated;
}

// Get the rotation 2x2 matrix in the trasnverse plane using the siagonalized trasnverse sphericity
std::array<std::array<double, 2>, 2> GetRotationTransverseSpheri(TClonesArray *particles) {
    // Compute the sphericity matrix
    double sphi[2][2] = {{0, 0}, {0, 0}};
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (particle->GetStatusCode() == 1) {
            sphi[0][0] += particle->Px() * particle->Px();
            sphi[0][1] += particle->Px() * particle->Py();
            sphi[1][0] += particle->Py() * particle->Px();
            sphi[1][1] += particle->Py() * particle->Py();
        }
    }

    // Compute the eigenvalues
    double ecc = pow((sphi[0][0] - sphi[1][1])*(sphi[0][0] - sphi[1][1]) + 4 * sphi[0][1]*sphi[0][1], 0.5);
    double eval1 = 0.5 * (sphi[0][0] + sphi[1][1] + sgn(sphi[0][0] - sphi[1][1]) * ecc);
    double eval2 = 0.5 * (sphi[0][0] + sphi[1][1] - sgn(sphi[0][0] - sphi[1][1]) * ecc);

    // Compute the rotation angle in the azimutal direction
    double phi = atan((std::max(eval1, eval2) - sphi[0][0])/sphi[0][1]);

    // printf("phiii  %.2f\n", phi * 180 / 3.14);
    // Compute the inverse rotation matrix
    std::array<std::array<double, 2>, 2> rot;
    rot[0] = {cos(phi), sin(phi)};
    rot[1] = {-sin(phi), cos(phi)};

    for (int i= 0; i<2; i++){
        for (int j= 0; j<2; j++){
            // printf("%.3f  ", rot[i][j]);
        }
        // printf("\n");
    }
    return rot;
}

// Compute the rotation
std::array<std::array<double, 2>, 2> GetRotationXYLeadingTrack(TClonesArray *particles) {
    // Compute the sphericity matrix
    double leadingPhi = 0;
    double leadingPxy = 0;
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        
        if (particle->GetStatusCode() == 1) {
            if (pow(particle->Px() * particle->Px() + particle->Py() *particle->Py(), 0.5) > leadingPxy) {
                leadingPhi = atan(particle->Py()/particle->Px());
            }
        }
    }

    // Compute the inverse rotation matrix
    std::array<std::array<double, 2>, 2> rot;
    rot[0] = {cos(leadingPhi), -sin(leadingPhi)};
    rot[1] = {sin(leadingPhi), cos(leadingPhi)};

    return rot;
}

// Pythia ######################################################################
// returns true if the particle is stable and charged
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

// Compute the multiplicity of charged particles in the final state at mid rapidity
int GetMult(TClonesArray *particles) {
    // evaluate multiplicity at forward rapidity
    int nParticles = 0;
    for (auto iPart = 2; iPart < particles->GetEntriesFast(); ++iPart) {
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        int pdg = std::abs(particle->GetPdgCode());
        int status = std::abs(particle->GetStatusCode());
        float eta = particle->Eta();

        if (IsDetectable(pdg) && std::abs(eta) < 0.8 && status == 1) {  // V0A and V0C acceptance
            nParticles++;
        }
    }
    return nParticles;
}

// compute transverse sphericity
double ComputeTransverseSphericity(TClonesArray *particles) {
    // Compute the sphericity matrix
    double sphi[2][2] = {{0, 0}, {0, 0}};
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        sphi[0][0] += particle->Px() * particle->Px();
        sphi[0][1] += particle->Px() * particle->Py();
        sphi[1][0] += particle->Py() * particle->Px();
        sphi[1][1] += particle->Py() * particle->Py();
    }

    // Compute the eigenvalues
    double ecc = pow((sphi[0][0] - sphi[1][1])*(sphi[0][0] - sphi[1][1]) + 4 * sphi[0][1]*sphi[0][1], 0.5);
    double eval1 = 0.5 * (sphi[0][0] + sphi[1][1] + sgn(sphi[0][0] - sphi[1][1]) * ecc);
    double eval2 = 0.5 * (sphi[0][0] + sphi[1][1] - sgn(sphi[0][0] - sphi[1][1]) * ecc);

    return 2*std::min(eval1, eval2)/(eval1 + eval2);
}

// Compute k*
float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2);

// Compute k* as in FD
float ComputeKstarFD(ROOT::Math::PxPyPzMVector PartOne, ROOT::Math::PxPyPzMVector PartTwo);

// Compute the rotation matrix for the event based on the sphericity (no minor axis rotation)
std::array<std::array<double, 3>, 3> GetRotationSpheri3D(TClonesArray *particles) {
    double sphi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (particle->GetStatusCode() == 1 && std::abs(particle->Eta()) < 0.8 && particle->Pt() > 0.3) {
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
    double c = sphi[0][0]*sphi[1][1] + sphi[0][0]*sphi[2][2] + sphi[1][1]*sphi[2][2] - sphi[0][1]*sphi[0][1] - sphi[0][2]*sphi[0][2] - sphi[1][2]*sphi[1][2];
    double d = sphi[0][0]*sphi[1][2]*sphi[1][2] + sphi[1][1]*sphi[0][2]*sphi[0][2] + sphi[2][2]*sphi[0][1]*sphi[0][1] - sphi[0][0]*sphi[1][1]*sphi[2][2] - 2*sphi[0][1]*sphi[0][2]*sphi[1][2];
    double p = b*b - 3*c;
    double q = 2*b*b*b - 9*b*c - 27*d;
    double delta = acos(q/2/pow(p, 1.5));

    // compute the eigenvalues
    double l1 = 1./3 * ( b+ 2* pow(p, 0.5) * cos(delta/3));
    double l2 = 1./3 * ( b+ 2* pow(p, 0.5) * cos((delta + 2*TMath::Pi())/3));
    double l3 = 1./3 * ( b+ 2* pow(p, 0.5) * cos((delta - 2*TMath::Pi())/3));

    // printf("eigenvectors: %.3f  %.3f  %.3f\n", l1, l2, l3);

    // compute the largest eigenvector -> principal axis of the event
    double lMax = std::max(std::max(l1, l2), l3);
    
    double evecMax[3];
    evecMax[0] = 1;
    evecMax[1] = (sphi[1][0] + (sphi[1][2]*(lMax - sphi[0][0]))/(sphi[0][2]))/(lMax + (sphi[1][2] * sphi[0][1])/(sphi[0][2]) - sphi[1][1]);
    evecMax[2] = (lMax - sphi[0][0] - sphi[0][1]*evecMax[1])/(sphi[0][2]);

    float phi = atan(evecMax[1]);
    float theta = atan(pow(evecMax[0]*evecMax[0] + evecMax[1]*evecMax[1], 0.5) / evecMax[2]);

    float ct = cos(theta);
    float st = sin(theta);
    float cp = cos(phi);
    float sp = sin(phi);

    std::array<std::array<double, 3>, 3> rot;
    rot[0] = { ct*cp, ct*sp, -st};
    rot[1] = { -sp  , cp   , 0  };
    rot[2] = { st*cp, st*sp, ct };

    return rot;
}

// Compute the rotation matrix for the event based on the sphericity (with minor axis included)
std::array<std::array<double, 3>, 3> GetRotationSpheri3DFull(TClonesArray *particles) {
    double sphi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (particle->GetStatusCode() == 1 && std::abs(particle->Eta()) < 0.8 && particle->Pt() > 0.3) {
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
    double c = sphi[0][0]*sphi[1][1] + sphi[0][0]*sphi[2][2] + sphi[1][1]*sphi[2][2] - sphi[0][1]*sphi[0][1] - sphi[0][2]*sphi[0][2] - sphi[1][2]*sphi[1][2];
    double d = sphi[0][0]*sphi[1][2]*sphi[1][2] + sphi[1][1]*sphi[0][2]*sphi[0][2] + sphi[2][2]*sphi[0][1]*sphi[0][1] - sphi[0][0]*sphi[1][1]*sphi[2][2] - 2*sphi[0][1]*sphi[0][2]*sphi[1][2];
    double p = b*b - 3*c;
    double q = 2*b*b*b - 9*b*c - 27*d;
    double delta = acos(q/2/pow(p, 1.5));

    // compute the eigenvalues
    double l1 = 1./3 * ( b+ 2* pow(p, 0.5) * cos(delta/3));
    double l2 = 1./3 * ( b+ 2* pow(p, 0.5) * cos((delta + 2*TMath::Pi())/3));
    double l3 = 1./3 * ( b+ 2* pow(p, 0.5) * cos((delta - 2*TMath::Pi())/3));

    // printf("eigenvectors: %.3f  %.3f  %.3f\n", l1, l2, l3);

    // compute the largest eigenvector -> principal axis of the event
    // sort eigenvalues
    std::vector<double> eigenvalues = {l1, l2, l3};
    std::sort(eigenvalues.begin(), eigenvalues.end());
    double lMed = eigenvalues[1];
    double lMax = eigenvalues[2];

    double evecMax[3];
    evecMax[0] = 1;
    evecMax[1] = (sphi[1][0] + (sphi[1][2]*(lMax - sphi[0][0]))/(sphi[0][2]))/(lMax + (sphi[1][2] * sphi[0][1])/(sphi[0][2]) - sphi[1][1]);
    evecMax[2] = (lMax - sphi[0][0] - sphi[0][1]*evecMax[1])/(sphi[0][2]);

    float phi = atan(evecMax[1]);
    float theta = atan(pow(evecMax[0]*evecMax[0] + evecMax[1]*evecMax[1], 0.5) / evecMax[2]);

    // float ct = cos(theta);
    // float st = sin(theta);
    // float cp = cos(phi);
    // float sp = sin(phi);

    std::array<std::array<double, 3>, 3> rot = Mul(RotY(-theta), RotZ(-phi));
    // rot[0] = { ct*cp, ct*sp, -st};
    // rot[1] = { -sp  , cp   , 0  };
    // rot[2] = { st*cp, st*sp, ct };


    std::array<double, 3> evecMed = ComputeEigenVector(sphi, lMed);

    std::array<double, 3> evecMedRot = Mul(rot, evecMed);
    // evecMedRot[0] = rot[0][0] * evecMed[0] + rot[0][1] * evecMed[1] + rot[0][2] * evecMed[2];
    // evecMedRot[1] = rot[1][0] * evecMed[0] + rot[1][1] * evecMed[1] + rot[1][2] * evecMed[2];
    // evecMedRot[2] = rot[2][0] * evecMed[0] + rot[2][1] * evecMed[1] + rot[2][2] * evecMed[2];

    double phiSecondaryAxis = atan(evecMedRot[1]/evecMedRot[0]);
    std::array<std::array<double, 3>, 3> rotSecondaryAxis = RotZ(-phiSecondaryAxis);
    
    return Mul(rotSecondaryAxis, rot);
}

// Compute the eigenvectors using the 3D sphericity matrix
std::array<std::array<double, 3>, 3> ComputeEigenVectorsSpheri3D(TClonesArray *particles) {
    double sphi[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (particle->GetStatusCode() == 1 && std::abs(particle->Eta()) < 0.8 && particle->Pt() > 0.3) {
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
    double c = sphi[0][0]*sphi[1][1] + sphi[0][0]*sphi[2][2] + sphi[1][1]*sphi[2][2] - sphi[0][1]*sphi[0][1] - sphi[0][2]*sphi[0][2] - sphi[1][2]*sphi[1][2];
    double d = sphi[0][0]*sphi[1][2]*sphi[1][2] + sphi[1][1]*sphi[0][2]*sphi[0][2] + sphi[2][2]*sphi[0][1]*sphi[0][1] - sphi[0][0]*sphi[1][1]*sphi[2][2] - 2*sphi[0][1]*sphi[0][2]*sphi[1][2];
    double p = b*b - 3*c;
    double q = 2*b*b*b - 9*b*c - 27*d;
    double delta = acos(q/2/pow(p, 1.5));

    // compute the eigenvalues
    double l1 = 1./3 * ( b+ 2* pow(p, 0.5) * cos(delta/3));
    double l2 = 1./3 * ( b+ 2* pow(p, 0.5) * cos((delta + 2*TMath::Pi())/3));
    double l3 = 1./3 * ( b+ 2* pow(p, 0.5) * cos((delta - 2*TMath::Pi())/3));

    // printf("eigenvectors: %.3f  %.3f  %.3f\n", l1, l2, l3);

    // compute the largest eigenvector -> principal axis of the event
    std::vector<double> eivals = {l1, l2, l3};
    std::sort(eivals.begin(), eivals.end());
    double lMin = eivals[0];
    double lMed = eivals[1];
    double lMax = eivals[2];
    
    std::array<double, 3> evecMin = ComputeEigenVector(sphi, lMin);
    std::array<double, 3> evecMed = ComputeEigenVector(sphi, lMed);
    std::array<double, 3> evecMax = ComputeEigenVector(sphi, lMax);

    std::array<std::array<double, 3>, 3> vec = {evecMax, evecMed, evecMin};

    // printf("Eigen vec:\n");
    // printf("lMin (%.3f): %.3f  %.3f  %.3f  \n", eivals[0], evecMin[0], evecMin[1], evecMin[2]);
    // printf("lMed (%.3f): %.3f  %.3f  %.3f  \n", eivals[1], evecMed[0], evecMed[1], evecMed[2]);
    // printf("lMax (%.3f): %.3f  %.3f  %.3f  \n", eivals[2], evecMax[0], evecMax[1], evecMax[2]);

    return vec;
}

// Compute the 3x3 sphericity matrix
std::array<std::array<double, 3>, 3> ComputeSphericityMatrix(TClonesArray *particles) {
    // Compute the sphericity matrix
    std::array<std::array<double, 3>, 3> sphi = {(std::array<double, 3>){0, 0, 0}, (std::array<double, 3>){0, 0, 0}, (std::array<double, 3>){0, 0, 0}};
    for (int iPart = 2; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
        if (particle->GetStatusCode() == 1 && std::abs(particle->Eta()) < 0.8 && particle->Pt() > 0.3 && IsDetectable(std::abs(particle->GetPdgCode()))) {
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

std::array<std::array<double, 3>, 3> ComputeSphericityMatrixAll(TClonesArray *particles) {
    // Compute the sphericity matrix
    std::array<std::array<double, 3>, 3> sphi = {(std::array<double, 3>){0, 0, 0}, (std::array<double, 3>){0, 0, 0}, (std::array<double, 3>){0, 0, 0}};
    for (int iPart = 0; iPart < particles->GetEntriesFast(); iPart++){
        TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));
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

// Diagonalize a 3x3 symmetric matrix and returns an array with its eigenvalues
std::array<double, 3> ComputeEigenValues(std::array<std::array<double, 3>, 3> mat) {
    double b = mat[0][0] + mat[1][1] + mat[2][2];
    double c = mat[0][0]*mat[1][1] + mat[0][0]*mat[2][2] + mat[1][1]*mat[2][2] - mat[0][1]*mat[0][1] - mat[0][2]*mat[0][2] - mat[1][2]*mat[1][2];
    double d = mat[0][0]*mat[1][2]*mat[1][2] + mat[1][1]*mat[0][2]*mat[0][2] + mat[2][2]*mat[0][1]*mat[0][1] - mat[0][0]*mat[1][1]*mat[2][2] - 2*mat[0][1]*mat[0][2]*mat[1][2];
    double p = b*b - 3*c;
    double q = 2*b*b*b - 9*b*c - 27*d;
    double delta = acos(q/2/pow(p, 1.5));

    // compute the eigenvalues
    std::array<double, 3> eivals;
    eivals[0] = 1./3 * (b + 2* pow(p, 0.5) * cos(delta/3));
    eivals[1] = 1./3 * (b + 2* pow(p, 0.5) * cos((delta + 2*TMath::Pi())/3));
    eivals[2] = 1./3 * (b + 2* pow(p, 0.5) * cos((delta - 2*TMath::Pi())/3));
    std::sort(eivals.begin(), eivals.end());
    return eivals;
}

// Compute the sphericity using the minimum eigenvalue of the sphericity matrix
double ComputeMajorSphericity(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    return 3 * eivals[0] / (eivals[0] + eivals[1] + eivals[2]);
}

// sphericity for the biggest and mid direction
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

double ComputeMinorRelativeSphericity(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    return 2 * eivals[0] / (eivals[0] + eivals[1]);
}

// Compute the direction along which the event develops
double ComputePrincipalAxisTheta(TClonesArray *particles) {
    auto sphi = ComputeSphericityMatrix(particles);
    auto eivals = ComputeEigenValues(sphi);
    std::array<double, 3> principalAxis = ComputeEigenVector(sphi, std::max(std::max(eivals[0], eivals[1]), eivals[2]));
    double theta = atan(principalAxis[2]/pow(principalAxis[0]*principalAxis[0] + principalAxis[1]*principalAxis[1], 0.5));
    return theta;
}


void MiniJetFactorization(
    int nEvents=1001,
    int pdg1=2212,
    int pdg2=-2212,
    int seed=42,
    FactorizationMode alignment = kNo,
    ComputeKStarMode computeKStarMode = kCustom,
    double minMult = 1,
    double maxMult = 1.e10,
    double minMajorSpheri = 0,
    double maxMajorSpheri = 1,
    double minMinorSpheri = 0,
    double maxMinorSpheri = 1,
    double minLambdaSpheriMax = 0,
    double maxLambdaSpheriMax = 1.e10,
    double minLambdaSpheriMid = 0,
    double maxLambdaSpheriMid = 1.e10,
    double minLambdaSpheriMin = 0,
    double maxLambdaSpheriMin = 1.e10,
    double minPrincipalAxisTheta = -3.14,
    double maxPrincipalAxisTheta = +3.14,
    std::string outFileNameRoot="AnalysisResults.root")
{
    //_________
    gSystem->Load("liblhapdf.so");
    gSystem->Load("libpythia8.so");
    gSystem->Load("libAliPythia8.so");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH", gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

    AliPythia8 pythia;
    // pythia.ReadString("SoftQCD:nonDiffractive = on");
    pythia.ReadString("SoftQCD:all = on");
    pythia.ReadString(Form("Tune:pp = 14"));
    pythia.ReadString("Random:setSeed = on");
    pythia.ReadString(Form("Random:seed %d", seed));
    pythia.Initialize(2212, 2212, 13000);

    // define outputs
    char* name = (char*)"";
    char* title = (char*)";#it{k}* (GeV/#it{c});pairs";

    TH1F *hSE = new TH1F("hSE", title, 3000, 0., 3.); 
    TH1F *hME = new TH1F("hME", title, 3000, 0., 3.); 
    TH2F *hSphi = new TH2F("hSphi", ";Major Sphericity;Minor Sphericity;Counts", 100, 0., 1, 100, 0, 1); 
    TH3F *hEivals = new TH3F("hEivals", ";max;mid;min;Counts", 100, 0., 100, 100, 0, 50,100, 0, 10); 
    TH1F *hMult = new TH1F("hMult", ";Multiplicity;Counts", 50, 0., 50); 
    TH2F *hEvtPartMult = new TH2F("hEvtPartMult", ";Event Multiplicity;Particle multiplicity;Counts", 50, 0., 50, 50, 0., 50); 
    TH1F *hPx = new TH1F("hPx", ";Px;Counts", 200, -2., 2); 
    TH1F *hPy = new TH1F("hPy", ";Py;Counts", 200, -2., 2); 
    TH1F *hPz = new TH1F("hPz", ";Pz;Counts", 200, -20., 20); 
    TH1F *hPrincipalAxisTheta = new TH1F("hPrincipalAxisTheta", ";PrincipalAxisTheta;Counts", 200, -1, 1); 
    //__________________________________________________________
    // perform the simulation
    std::vector<ROOT::Math::PxPyPzMVector> part1{};
    std::vector<ROOT::Math::PxPyPzMVector> part2{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBuffer1{};
    std::deque<std::vector<ROOT::Math::PxPyPzMVector>> partBuffer2{};


    float (*ComputeKStarFunction) (ROOT::Math::PxPyPzMVector PartOne, ROOT::Math::PxPyPzMVector PartTwo);
    if (computeKStarMode == ComputeKStarMode::kCustom) {
        ComputeKStarFunction = ComputeKstar;
    } else if (computeKStarMode == ComputeKStarMode::kFD) {
        ComputeKStarFunction = ComputeKstarFD;
    } else {
        printf("KStar function not implemented. Exit\n");
        exit(1);
    }


    TFile outFile(outFileNameRoot.data(), "recreate");

    TClonesArray* particles = new TClonesArray("TParticle", 1000);
    for (int iEvent=0; iEvent<nEvents; iEvent++)
    {
        part1.clear();
        part2.clear();
        
        pythia.GenerateEvent();
        pythia.ImportParticles(particles, "All");

        // Compute event properties
        int mult = GetMult(particles);
        // if (mult < 3) continue;
        double spheriMaj = ComputeMajorRelativeSphericity(particles);
        double spheriMin = ComputeMinorRelativeSphericity(particles);
        double principalAxisTheta = ComputePrincipalAxisTheta(particles);

        auto eivals = ComputeEigenValues(ComputeSphericityMatrix(particles));

        // Event selection
        if (
            !(minMult < mult && mult < maxMult &&
              minMajorSpheri < spheriMaj && spheriMaj < maxMajorSpheri &&             
              minMinorSpheri < spheriMin && spheriMin < maxMinorSpheri &&
              minLambdaSpheriMax < eivals[2] && eivals[2] < maxLambdaSpheriMax &&             
              minLambdaSpheriMid < eivals[1] && eivals[1] < maxLambdaSpheriMid &&
              minLambdaSpheriMin < eivals[0] && eivals[0] < maxLambdaSpheriMin &&
              minPrincipalAxisTheta < principalAxisTheta && principalAxisTheta < maxPrincipalAxisTheta)
        ) continue;
        // if ( mult < 30 || mult > 40 || spheri < 0.2 || spheri > 0.4 || std::abs(principalAxisTheta) > 0.2) continue;
        hSphi->Fill(spheriMaj, spheriMin);
        hEivals->Fill(eivals[2], eivals[1], eivals[0]);

        int iPoint = 0;
        int iPointEvt = 0;

        double phiLeading = 0;
        double thetaLeading = 0;
        double pLeading = 0;

        TParticle *leadingParticle = nullptr;
        // TCanvas *cBefore = new TCanvas("cBefore", "", 600, 600);
        // TH3F* hEvent = new TH3F("hEvent", "", 10, -3, 3, 10, -3, 3, 10, -3, 3);
        // hEvent->Draw();
        // TPolyLine3D *lAxisX = new TPolyLine3D(2);       
        // lAxisX->SetPoint(0,0,0,0);
        // lAxisX->SetPoint(1,3, 0, 0);
        // lAxisX->SetLineWidth(3);
        // lAxisX->Draw();
        // TPolyLine3D *lAxisY = new TPolyLine3D(2);       
        // lAxisY->SetPoint(0,0,0,0);
        // lAxisY->SetPoint(1,0, 3, 0);
        // lAxisY->SetLineWidth(3);
        // lAxisY->Draw();
        // TPolyLine3D *lAxisZ = new TPolyLine3D(2);       
        // lAxisZ->SetPoint(0,0,0,0);
        // lAxisZ->SetPoint(1,0, 0, 3);
        // lAxisZ->SetLineWidth(3);
        // lAxisZ->Draw();

        for(int iPart=2; iPart<particles->GetEntriesFast(); iPart++) {
            
            TParticle* particle = dynamic_cast<TParticle*>(particles->At(iPart));

            if(std::abs(particle->Eta()) > 0.8) continue; // keep only midrapidity

            double prodR = std::sqrt(particle->Vx()*particle->Vx() + 
                                     particle->Vy()*particle->Vy() + 
                                     particle->Vz()*particle->Vz());

            int pdg = particle->GetPdgCode();
            // if(pdg == 3122) {
            //     printf("lallero %d\n", particle->GetStatusCode());
            // }
            // if(prodR > 1. ) continue;

            if(prodR > 1.) continue;
            // if(prodR > 1. || particle->GetStatusCode() != 1) continue;

            double p = particle->P();

            if (p > pLeading) {
                leadingParticle = particle;
                pLeading = p;
            }

            // TPolyLine3D *lP = new TPolyLine3D(2);       
            // lP->SetPoint(0,0,0,0);
            // lP->SetPoint(1,particle->Px(),particle->Py(),particle->Pz());

            if(pdg == pdg1) {
                double mass = TDatabasePDG::Instance()->GetParticle(std::abs(pdg))->Mass();
                ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(), mass);
                part1.push_back(part);
                // lP->SetLineColor(kBlue);
                // lP->SetLineWidth(2);
            }
            // printf("lam p pdg: %d\n", pdg);

            else if(pdg == pdg2) {
                double mass = TDatabasePDG::Instance()->GetParticle(std::abs(pdg))->Mass();
                ROOT::Math::PxPyPzMVector part(particle->Px(), particle->Py(), particle->Pz(), mass);
                part2.push_back(part);
                // printf("lam\n");
                // lP->SetLineColor(kBlue);
                // lP->SetLineWidth(2);
            } else {
                // lP->SetLineColor(kGray+2);
            }
            // lP->DrawClone();
            // delete lP;
        }
        // if (part1.size() + part2.size() < 5) continue;
        hEvtPartMult->Fill(mult, part1.size() + part2.size());

        if (leadingParticle) {
            // TPolyLine3D *lLeading = new TPolyLine3D(2);       
            // lLeading->SetPoint(0,0,0,0);
            // lLeading->SetPoint(1,leadingParticle->Px(),leadingParticle->Py(),leadingParticle->Pz());
            // lLeading->SetLineColor(kRed);
            // lLeading->DrawClone();
            // delete lLeading;
        } else continue;




        // hSphi->Fill(spheri);
        hMult->Fill(mult);
        hPrincipalAxisTheta->Fill(principalAxisTheta);
        // printf("n1 %d\n", part1.size());
        // printf("n2 %d\n", part2.size());
        if (alignment == FactorizationMode::kXY) {
            std::array<std::array<double, 2>, 2> rot = GetRotationTransverseSpheri(particles);
            
            part1 = RotateXY(rot, part1);
            part2 = RotateXY(rot, part2);
        } else if (alignment == FactorizationMode::kLeadingXY) {
            std::array<std::array<double, 2>, 2> rot = GetRotationXYLeadingTrack(particles);
            
            part1 = RotateXY(rot, part1);
            part2 = RotateXY(rot, part2);

        } else if (alignment == FactorizationMode::kLeading) {
            phiLeading = atan(leadingParticle->Py()/leadingParticle->Px());
            if(leadingParticle->Px() < 0) phiLeading += TMath::Pi();
            std::array<std::array<double, 3>, 3> rot;

            rot[0] = { cos(phiLeading),  sin(phiLeading), 0};
            rot[1] = { -sin(phiLeading),  cos(phiLeading), 0};
            rot[2] = {      0         ,        0        , 1};

            part1 = Rotate(rot, part1);
            part2 = Rotate(rot, part2);

            thetaLeading = atan(leadingParticle->Pt()/leadingParticle->Pz());
            if(leadingParticle->Pz() < 0) thetaLeading += TMath::Pi();

            rot[0] = {  cos(thetaLeading),   0, -sin(thetaLeading) };
            rot[1] = {        0        ,     1,        0        };
            rot[2] = { sin(thetaLeading),    0, cos(thetaLeading) };

            part1 = Rotate(rot, part1);
            part2 = Rotate(rot, part2);
        } else if (alignment == FactorizationMode::kSpheri3D) {
            auto eivec = ComputeEigenVectorsSpheri3D(particles);

            // TPolyLine3D *lEv1 = new TPolyLine3D(2);       
            // lEv1->SetPoint(0,0,0,0);
            // lEv1->SetPoint(1,eivec[0][0],eivec[0][1],eivec[0][2]);
            // lEv1->SetLineColor(kGreen);
            // lEv1->SetLineWidth(2);
            // lEv1->DrawClone();
            // delete lEv1;

            // TPolyLine3D *lEv2 = new TPolyLine3D(2);       
            // lEv2->SetPoint(0,0,0,0);
            // lEv2->SetPoint(1,eivec[1][0],eivec[1][1],eivec[1][2]);
            // lEv2->SetLineColor(kGreen);
            // lEv2->SetLineWidth(2);
            // lEv2->DrawClone();
            // delete lEv2;

            // TPolyLine3D *lEv3 = new TPolyLine3D(2);
            // lEv3->SetPoint(0,0,0,0);
            // lEv3->SetPoint(1,eivec[2][0],eivec[2][1],eivec[2][2]);
            // lEv3->SetLineColor(kGreen);
            // lEv3->SetLineWidth(2);
            // lEv3->DrawClone();
            // delete lEv3;

            std::array<std::array<double, 3>, 3> rot = GetRotationSpheri3D(particles);
            part1 = Rotate(rot, part1);
            part2 = Rotate(rot, part2);
        } else if (alignment == FactorizationMode::kSpheri3DFull) {
            std::array<std::array<double, 3>, 3> rot = GetRotationSpheri3DFull(particles);
            part1 = Rotate(rot, part1);
            part2 = Rotate(rot, part2);
        }


        
        // cBefore->Update();
        // cBefore->Write();
        
        // TCanvas *cAfter = new TCanvas("cAfter", "", 600, 600);
        // hEvent->Draw();
        // lAxisX->Draw();
        // lAxisY->Draw();
        // lAxisZ->SetLineWidth(1);
        // lAxisZ->Draw();
        // for (const auto &p : part1) {
        //     TPolyLine3D *lParticle = new TPolyLine3D(2);       
        //     lParticle->SetPoint(0,0,0,0);
        //     lParticle->SetPoint(1,p.Px(), p.Py(), p.Pz());
        //     lParticle->SetLineColor(kBlue);
        //     lParticle->DrawClone();   
        //     delete lParticle;
        // }
        // for (const auto &p : part2) {
        //     TPolyLine3D *lParticle = new TPolyLine3D(2);       
        //     lParticle->SetPoint(0,0,0,0);
        //     lParticle->SetPoint(1,p.Px(), p.Py(), p.Pz());
        //     lParticle->SetLineColor(kBlue);
        //     lParticle->DrawClone();
        //     delete lParticle;
        // }
        // cAfter->Update();
        // cAfter->Write();

        // delete hEvent;
        // delete cBefore;
        // delete cAfter;
        

        // break;
        for (const auto &p1 : part1) {
            hPx->Fill(p1.Px());
            hPy->Fill(p1.Py());
            hPz->Fill(p1.Pz());
        }
        partBuffer1.push_back(part1);
        partBuffer2.push_back(part2);

        if (partBuffer1.size() > 50) { //buffer full, let's kill the first entry
            partBuffer1.pop_front();
            partBuffer2.pop_front();
        }

        // same event
        for(int iP1 = 0; iP1 < part1.size(); iP1++) {
            if (pdg1 == pdg2) {
                for(int iP2 = iP1+1; iP2 < part1.size(); iP2++) {
                    double kStar = ComputeKStarFunction(part1[iP1], part1[iP2]);
                    hSE->Fill(kStar);
                }    
            } else {
                for(int iP2 = 0; iP2 < part2.size(); iP2++) {
                    double kStar = ComputeKStarFunction(part1[iP1], part2[iP2]);
                    hSE->Fill(kStar);
                }
            }
        }

        // mixed event
        if(partBuffer1.size() < 2) // to avoid repetitions
            continue;

        if (pdg1 == pdg2) {
            for(size_t ip2=0; ip2<partBuffer1[partBuffer1.size()-1].size(); ip2++) // last only
            {
                for(size_t iME=0; iME<partBuffer1.size()-1; iME++) // from 0 to last-1
                {
                    for(size_t ip1=0; ip1<partBuffer1[iME].size(); ip1++)
                    {
                        hME->Fill(ComputeKStarFunction(partBuffer1[partBuffer1.size()-1][ip2], partBuffer1[iME][ip1]));
                    }
                }
            }
        } else {
            for(size_t ip2=0; ip2<partBuffer2[partBuffer2.size()-1].size(); ip2++) // last only
            {
                for(size_t iME=0; iME<partBuffer1.size()-1; iME++) // from 0 to last-1
                {
                    for(size_t ip1=0; ip1<partBuffer1[iME].size(); ip1++)
                    {
                        hME->Fill(ComputeKStarFunction(partBuffer2[partBuffer2.size()-1][ip2], partBuffer1[iME][ip1]));
                    }
                }
            }
        }
    }

    // save root output file
    hSE->Write();
    hME->Write();
    hSphi->Write();
    hEivals->Write();
    hMult->Write();
    hEvtPartMult->Write();
    hPx->Write();
    hPy->Write();
    hPz->Write();
    hPrincipalAxisTheta->Write();
    outFile.Close();
}

float ComputeKstar(ROOT::Math::PxPyPzMVector part1, ROOT::Math::PxPyPzMVector part2)
{
    ROOT::Math::PxPyPzMVector trackSum = part1 + part2;
    ROOT::Math::Boost boostv12{trackSum.BoostToCM()};
    ROOT::Math::PxPyPzMVector part1CM = boostv12(part1);
    ROOT::Math::PxPyPzMVector part2CM = boostv12(part2);

    ROOT::Math::PxPyPzMVector trackRelK = part1CM - part2CM;
    float kStar = 0.5 * trackRelK.P();
    return kStar;
}


float ComputeKstarFD(ROOT::Math::PxPyPzMVector PartOne, ROOT::Math::PxPyPzMVector PartTwo) {
  ROOT::Math::PxPyPzMVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartOneCMS = TLorentzVector(PartOne.Px(),PartOne.Py(),PartOne.Pz(), PartOne.E());
  TLorentzVector PartTwoCMS = TLorentzVector(PartTwo.Px(),PartTwo.Py(),PartTwo.Pz(), PartTwo.E());

  PartOneCMS.Boost(-betax, -betay, -betaz);
  PartTwoCMS.Boost(-betax, -betay, -betaz);

  TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;

  return 0.5 * trackRelK.P();
}
