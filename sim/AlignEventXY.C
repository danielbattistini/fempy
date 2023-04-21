#include <cmath>

int sgn(double value) { return (value > 0) - (value < 0); }

double *prod(double mat[2][2], double vec[2]) {
    static double ans[2];
    ans[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1];
    ans[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1];
    return ans;
}

void print(double mat[2][2]){
    for (int i=0; i<2; i++){
        for (int j=0; j<2; j++){
            printf("%.3f  ", mat[i][j]);
        }
        printf("\n");
    }
}

void print(double vec[2]){
    printf("%.3f  %.3f\n", vec[0], vec[1]);
}

void AlignEventXY() {
    std::vector<std::vector<double>> p = {
        {-0.162931, -0.372415},
        {-0.588803, +0.007411},
        {-0.065611, +0.495443},
        {+0.053690, +0.088543},
        {-0.465681, +0.199942},
        {-0.599639, +0.035467},
        {-0.332850, -0.267899},
        {-0.767986, +0.251907},
        {-0.303869, +0.242914},
        {+0.336314, -0.073045},
        {+0.274709, -0.006469},
        {+0.492120, -0.713790},
    };

    p = {
        {-0.4, -0.1},
        {-0.2, +0.2},
        {-0.1, +0.2},
        {+0.1, +0.4},
        {-1, -1},
        {-1.1, -1},
    };

    TGraph *gEvent = new TGraph(1);
    for (int n = 0; n < p.size(); n++) {
        gEvent->SetPoint(n, p[n][0], p[n][1]);
    }

    TCanvas *cEvent = new TCanvas("cEvent", "", 600, 600);
    cEvent->DrawFrame(-2, -2, 2, 2);
    cEvent->SetGridx(true);
    cEvent->SetGridy(true);
    gEvent->SetMarkerStyle(20);
    gEvent->Draw("p same");

    // make spheri matrix
    double spheri[2][2] = {{0, 0}, {0, 0}};
    for (int n = 0; n < p.size(); n++){
        spheri[0][0] += p[n][0] * p[n][0];
        spheri[0][1] += p[n][0] * p[n][1];
        spheri[1][0] += p[n][1] * p[n][0];
        spheri[1][1] += p[n][1] * p[n][1];
    }
    // spheri[0][0] = 6;
    // spheri[0][1] = -2;
    // spheri[1][0] = -2;
    // spheri[1][1] = 3;

    print(spheri);

    double l1 = 0.5 * (spheri[0][0] + spheri[1][1] + sgn(spheri[0][0] - spheri[1][1]) * pow((spheri[0][0] - spheri[1][1])*(spheri[0][0] - spheri[1][1]) + 4 * spheri[0][1]*spheri[0][1], 0.5));
    double l2 = 0.5 * (spheri[0][0] + spheri[1][1] - sgn(spheri[0][0] - spheri[1][1]) * pow((spheri[0][0] - spheri[1][1])*(spheri[0][0] - spheri[1][1]) + 4 * spheri[0][1]*spheri[0][1], 0.5));
    printf("l1: %f   e2: %f\n", l1, l2);
    double ev1[2];
    ev1[0] = 1;
    ev1[1] = (l1 - spheri[0][0])/spheri[0][1];
    double ev2[2];
    ev2[0] = 1;
    ev2[1] = (l2 - spheri[0][0])/spheri[0][1];

    printf("eigen\n");
    print(ev1);
    print(ev2);

    printf("\nprod\n");
    print(prod(spheri, ev1));
    print(prod(spheri, ev2));

    double bigv[2];
    bigv[0] = l1 > l2 ? ev1[0] : ev2[0];
    bigv[1] = l1 > l2 ? ev1[1] : ev2[1];
    double phi = atan(bigv[1]/bigv[0]);
    printf("phi = %.2f deg\n\n", phi * 180 / 3.14);

    // rotate
    double rot[2][2];
    rot[0][0] = cos(phi);
    rot[0][1] = -sin(phi);
    rot[1][0] = sin(phi);
    rot[1][1] = cos(phi);

    std::vector<std::vector<double>> pRot;
    for (const auto &part : p){
        double pp[2] = {part[0], part[1]};
        auto rotated = prod(rot, pp);
  
        pRot.push_back({rotated[0], rotated[1]});
    }


    TGraph *gRot = new TGraph(1);
    for (int n = 0; n < p.size(); n++) {
        gRot->SetPoint(n, pRot[n][0], pRot[n][1]);
    }
    TCanvas *cRot = new TCanvas("cRot", "", 600, 600);
    cRot->DrawFrame(-2, -2, 2, 2);
    cRot->SetGridx(true);
    cRot->SetGridy(true);
    gRot->SetMarkerStyle(20);
    gRot->Draw("p same");
}
