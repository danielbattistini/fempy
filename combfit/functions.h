#ifndef COMBFIT_FUNCTIONS_H_
#define COMBFIT_FUNCTIONS_H_

#include <gsl/gsl_sf_dawson.h>

const double Pi(3.141592653589793);
const double AlphaFS(0.0072973525664);
const double RevSqrt2(0.7071067811865475);
const std::complex<double> i(0, 1);
const double EulerConst = 0.57721566490153;
const double FmToNu(5.067731237e-3);
const double NuToFm(197.3269602);
const double hbarc(197.3269602);
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

    if (QS) {
        CkValue -= 0.5 * exp(-4. * Rho2);
    } else {
        CkValue *= 2.;
    }

    return A_c * (CkValue + 1.);
}

double GeneralCoulombLednickySecond(const double &Momentum, const double &GaussR, const double &ScattLenSin,
                                    const double &EffRangeSin, const double &ScattLenTri, const double &EffRangeTri,
                                    const bool &QS, const double &RedMass, const double &Q1Q2) {
    return 1. / 3 * GeneralCoulombLednicky(Momentum, GaussR, ScattLenSin, EffRangeSin, QS, RedMass, Q1Q2) +
           2. / 3 * GeneralCoulombLednicky(Momentum, GaussR, ScattLenTri, EffRangeTri, QS, RedMass, Q1Q2);
}

/*
x[0]:  Momentum (double)

par[0]: GaussR (double)
par[1]: ScattLenSin (double)
par[2]: EffRangeSin (double)
par[3]: QS (bool)
par[4]: RedMass (double)
par[5]: Q1Q2 (double)
*/
double GeneralCoulombLednicky(double *x, double *pars) {
    return GeneralCoulombLednicky(x[0], pars[0], pars[1], pars[2], static_cast<bool>(pars[3]), pars[4], pars[5]);
}

/*
x[0]:   Momentum (double)

par[0]: GaussR (double)
par[1]: ScattLenSin (double)
par[2]: EffRangeSin (double)
par[3]: ScattLenTri (double)
par[4]: EffRangeTri (double)
par[5]: QS (bool)
par[6]: RedMass (double)
par[7]: Q1Q (double)
*/
double GeneralCoulombLednickySecond(double *x, double *pars) {
    return GeneralCoulombLednickySecond(x[0], pars[0], pars[1], pars[2], pars[3], pars[4], static_cast<bool>(pars[5]),
                                        pars[6], pars[7]);
}

/*
x[0]:   Momentum (double)
par[0]: r1 (double)
par[1]: r2 (double)
par[2]: w1 (double)
par[3]: ScattLen (double)
par[4]: EffRange (double)
par[5]: QS (bool)
par[6]: RedMass (double)
par[7]: Q1Q (double)
*/
double GeneralCoulombLednickyTwoRadii(double *x, double *pars) {
    return pars[2] *
               GeneralCoulombLednicky(x[0], pars[0], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]) +
           (1 - pars[2]) *
               GeneralCoulombLednicky(x[0], pars[1], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]);
}

/*
x[0]:   Momentum (double)
par[0]: r1 (double)
par[1]: r2 (double)
par[2]: w1 (double)
par[3]: ScattLen (double)
par[4]: EffRange (double)
par[5]: QS (bool)
par[6]: RedMass (double)
par[7]: Q1Q (double)
par[8]: overall normalization (double)
*/
double ScalableGeneralCoulombLednickyTwoRadii(double *x, double *pars) {
    double gcl1 = GeneralCoulombLednicky(x[0], pars[0], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]);
    double gcl2 = GeneralCoulombLednicky(x[0], pars[1], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]);
    return pars[8] * (pars[2] * gcl1 + (1 - pars[2]) * gcl2);
}

/*
x[0]:   Momentum (double)
par[0]: r1 (double)
par[1]: r2 (double)
par[2]: w1 (double)
par[3]: ScattLen (double)
par[4]: EffRange (double)
par[5]: QS (bool)
par[6]: RedMass (double)
par[7]: Q1Q (double)
par[8]: overall normalization (double)
*/
double ScalableShifted300GeneralCoulombLednickyTwoRadii(double *x, double *pars) {
    double gcl1 = GeneralCoulombLednicky(x[0], pars[0], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]);
    double gcl2 = GeneralCoulombLednicky(x[0], pars[1], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]);
    double cf = pars[8] * (pars[2] * gcl1 + (1 - pars[2]) * gcl2);
    // compute the CF at 300 MeV and shift
    double res1 = GeneralCoulombLednicky(300, pars[0], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]);
    double res2 = GeneralCoulombLednicky(300, pars[1], pars[3], pars[4], static_cast<bool>(pars[5]), pars[6], pars[7]);
    double residue = pars[8] * (pars[2] * res1 + (1 - pars[2]) * res2);

    return cf - (residue - 1);
}

/*
x[0]:   Momentum (double)

par[0]: r1 (double)
par[1]: r2 (double)
par[2]: w1 (double)
par[3]: ScattLenSin (double)
par[4]: EffRangeSin (double)
par[5]: ScattLenTri (double)
par[6]: EffRangeTri (double)
par[7]: QS (bool)
par[8]: RedMass (double)
par[9]: Q1Q (double)
*/
double GeneralCoulombLednickySecondTwoRadii(double *x, double *pars) {
    return pars[2] * GeneralCoulombLednickySecond(x[0], pars[0], pars[3], pars[4], pars[5], pars[6],
                                                  static_cast<bool>(pars[7]), pars[8], pars[9]) +
           (1 - pars[2]) * GeneralCoulombLednickySecond(x[0], pars[1], pars[3], pars[4], pars[5], pars[6],
                                                        static_cast<bool>(pars[7]), pars[8], pars[9]);
}

/*
x[0]:   Momentum (double)

par[0]: r1 (double)
par[1]: r2 (double)
par[2]: w1 (double)
par[3]: ScattLenSin (double)
par[4]: EffRangeSin (double)
par[5]: ScattLenTri (double)
par[6]: EffRangeTri (double)
par[7]: QS (bool)
par[8]: RedMass (double)
par[9]: Q1Q (double)
par[10]: overall normalization (double)
*/
double ScalableGeneralCoulombLednickySecondTwoRadii(double *x, double *pars) {
    double gcl1 = GeneralCoulombLednickySecond(x[0], pars[0], pars[3], pars[4], pars[5], pars[6],
                                               static_cast<bool>(pars[7]), pars[8], pars[9]);
    double gcl2 = GeneralCoulombLednickySecond(x[0], pars[1], pars[3], pars[4], pars[5], pars[6],
                                               static_cast<bool>(pars[7]), pars[8], pars[9]);
    return pars[10] * (pars[2] * gcl1 + (1 - pars[2]) * gcl2);
}
#endif  // COMBFIT_FUNCTIONS_H_
