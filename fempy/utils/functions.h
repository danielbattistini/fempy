#include <gsl/gsl_sf_dawson.h>

const double Pi(3.141592653589793);
const double AlphaFS(0.0072973525664);  // fine structure constant
const double RevSqrt2(0.7071067811865475);
const std::complex<double> i(0, 1);
const double EulerConst = 0.57721566490153;
const double FmToNu(5.067731237e-3);  // fm-1 MeV-1
const double NuToFm(197.3269602);     // fm MeV
const double hbarc(197.3269602);      // fm MeV
const double RadToDeg(57.295779513082);
const double DegToRad(0.017453292519943);

double CoulombEta(const double &Momentum, const double &RedMass, const double &Q1Q2) {
    if (!Momentum) return 0;

    return AlphaFS * RedMass * Q1Q2 / Momentum;
}

double CoulombPenetrationFactor(const double &eta) { return eta ? 2. * Pi * eta / (exp(2. * Pi * eta) - 1.) : 1; }

double CoulombEuler(const double &eta) {
    if (!eta) return 0;
    double RESULT = 0;
    double ADD;
    const double eta2 = eta * eta;
    for (double dIter = 1; dIter <= 11; dIter++) {
        ADD = 1. / (dIter * (dIter * dIter + eta2));
        RESULT += ADD;
        if (fabs(ADD / RESULT) < 1e-7) break;
    }
    RESULT *= eta2;
    RESULT -= log(eta2) + EulerConst;

    return RESULT;
}

double GeneralCoulombLednicky(const double &Momentum, const double &GaussR, const double &ScattLenSin,
                              const double &EffRangeSin, const bool &QS, const double &RedMass, const double &Q1Q2) {
    if (GaussR != GaussR) {
        std::cerr << "\033[1;33mWARNING:\033[0m GeneralCoulombLednicky got a bad "
                     "value for the Radius (nan). Returning default value of 1.\n"
                  << std::endl;
        return 1;
    }

    const double Momentum2 = Momentum * Momentum;
    const double Radius = GaussR * FmToNu;
    const double Rho = Radius * Momentum;
    const double Rho2 = Rho * Rho;
    const double sLen1 = ScattLenSin * FmToNu;
    const double eRan1 = EffRangeSin * FmToNu;

    double eta = CoulombEta(Momentum, RedMass, Q1Q2);
    double A_c = CoulombPenetrationFactor(eta);

    complex<double> ScattAmplSin =
        pow(1. / sLen1 + 0.5 * eRan1 * Momentum2 - i * Momentum * A_c - 2. * Momentum * eta * CoulombEuler(eta), -1.);

    double CkValue =
        0.25 * pow(abs(ScattAmplSin) / Radius, 2) *
            (1 - (eRan1) / (2 * sqrt(Pi) * Radius) + 0.5 * pow(A_c - 1., 2) * (1. - exp(-4. * Rho2))) +
        Momentum * real(ScattAmplSin) * gsl_sf_dawson(2. * Rho) / (2. * sqrt(Pi) * Rho2) -
        Momentum * imag(ScattAmplSin) * (0.25 * (1. - exp(-4. * Rho2)) / Rho2 + (A_c - 1.) * cos(Rho) * exp(-Rho2));

    if (QS) {  // quantum statistics
        CkValue -= 0.5 * exp(-4. * Rho2);
    } else {
        CkValue *= 2.;
    }

    return A_c * (CkValue + 1.);
}
lalala break 
double GeneralCoulombLednickySecond(const double &Momentum, const double &GaussR, const double &ScattLenSin,
                                    const double &EffRangeSin, const double &ScattLenTri, const double &EffRangeTri,
                                    const bool &QS, const double &RedMass, const double &Q1Q2) {
    return 0.25 * GeneralCoulombLednicky(Momentum, GaussR, ScattLenSin, EffRangeSin, QS, RedMass, Q1Q2) +
           0.75 * GeneralCoulombLednicky(Momentum, GaussR, ScattLenTri, EffRangeTri, QS, RedMass, Q1Q2);
}

double GeneralCoulombLednicky(double *x, double *pars) {
    return GeneralCoulombLednicky(x[0], pars[0], pars[1], pars[2], bool(pars[3]), pars[4], pars[5]);
}

double GeneralCoulombLednickySecond(double *x, double *pars) {
    return GeneralCoulombLednickySecond(x[0], pars[0], pars[1], pars[2], pars[3], pars[4], bool(pars[5]), pars[6],
                                        pars[7]);
}

// void TidyCats::GetCatsPionDstar(CATS *cats, int momBins, double kMin, double kMax, TidyCats::lightDmesonPot pot,
//                                 TidyCats::Sources source, int chargecombi) {
//     const auto pdgDatabase = TDatabasePDG::Instance();
//     const double massPion = pdgDatabase->GetParticle(211)->Mass() * 1000;
//     const double massDstar = pdgDatabase->GetParticle(chargecombi * 413)->Mass() * 1000;
//     CATSparameters *cPars;

//     cPars = new CATSparameters(CATSparameters::tSource, 1, true);
//     cPars->SetParameter(0, 0.8);
//     cats->SetAnaSource(GaussSource, *cPars);
//     cats->SetUseAnalyticSource(true);
//     cats->SetMomentumDependentSource(false);
//     cats->SetThetaDependentSource(false);
//     cats->SetMomBins(momBins, kMin, kMax);
//     cats->SetNumChannels(1);
//     cats->SetQuantumStatistics(0);
//     cats->SetNumPW(0, 1);
//     cats->SetSpin(0, 0);
//     cats->SetChannelWeight(0, 1.);

//     cats->SetQ1Q2(chargecombi);
//     cats->SetPdgId(321, (chargecombi * 413));
//     cats->SetRedMass((massPion * massDstar) / (massPion + massDstar));
//     DLM_Histo<complex<double>> ***ExternalWF = nullptr;

//     switch (pot) {
//         case TidyCats::lightDCoulombOnly:
//             std::cout << "Coulomb only\n";
//             break;
//         default:
//             std::cout << "Potential not implemented\n";
//             break;
//     }
//     return;
// }