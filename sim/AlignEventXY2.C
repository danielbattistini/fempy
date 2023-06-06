#include <cmath>
#include "functions.hxx"

void AlignEventXY() {
    std::vector<float> p[2] = {
        {-0.162931, -0.588803, -0.065611, 0.053690, -0.465681, -0.599639, -0.332850, -0.767986, -0.303869, 0.336314, 0.274709, 0.492120},
        {-0.372415, 0.007411, 0.495443, 0.088543, 0.199942, 0.035467, -0.267899, 0.251907, 0.242914, -0.073045, -0.006469, -0.713790},
    };

    TGraph *gEvent = new TGraph(1);
    for (int n = 0; n<12; n++) {
        gEvent->SetPoint(n, p[0][n], p[1][n]);
    }
    
    TCanvas *cEvent= new TCanvas("cEvent", "", 600, 600);
    cEvent->DrawFrame(-1, -1, 1, 1);
    gEvent->SetMarkerStyle(20);
    gEvent->Draw("p same");

    double ptot = 0;
    for (int n = 0; n<12; n++) {
        ptot += pow(p[0][n] * p[0][n] + p[1][n] * p[1][n], 0.5) ;
    }

    // compute spheri matrix
    double spheri[2][2] = { {0, 0}, {0, 0}};
    for (int i = 0; i<2; i++) {
        for (int j = 0; j<2; j++) {
            for (int n = 0; n<12; n++) {
                spheri[i][j] += p[i][n] * p[j][n] / pow(p[0][n] * p[0][n] + p[1][n] * p[1][n], 0.5) ;
            }
        }
    }
    for (int i = 0; i<2; i++) {
        for (int j = 0; j<2; j++) {
            spheri[i][j] /= ptot;
        }
    }
    double l1 = 0.5 * (spheri[0][0] + spheri[1][1] + sgn(spheri[0][0] - spheri[1][1]) * pow((spheri[0][0] - spheri[1][1])*(spheri[0][0] - spheri[1][1]) + 4 * spheri[0][1]*spheri[0][1], 0.5));
    double l2 = 0.5 * (spheri[0][0] + spheri[1][1] - sgn(spheri[0][0] - spheri[1][1]) * pow((spheri[0][0] - spheri[1][1])*(spheri[0][0] - spheri[1][1]) + 4 * spheri[0][1]*spheri[0][1], 0.5));
    printf("e1: %f     e2: %f\n", l1, l2);

    // eigenv
    double phi = 0.5 * atan(2. * spheri[0][1] / (spheri[0][0] - spheri[1][1]));
    double rot[2][2];
    rot[0][0] = rot[1][1] = cos(phi);
    rot[0][1] = sin(phi);
    rot[1][0] = -sin(phi);

    // double rot[2][2];
    // rot[0][0] = l1 * cos(phi) * cos(phi) + l2 * sin(phi) * sin(phi);
    // rot[0][1] = rot[1][0] = (l1 -l2) * cos(phi) * sin(phi);
    // rot[1][1] = l1 * sin(phi) * sin(phi) + l2 * cos(phi) * cos(phi);

    printf("fiii %f\n", phi);
    printf("fiii deg %f\n", 180 / 3.14 * phi);


    

    printf("\nspheri matrix\n");
    for (int i = 0; i<2; i++) {
        for (int j = 0; j<2; j++) {
            printf("%.2f ", spheri[i][j]/ptot);
        }
        printf("\n");
    }
    printf("%f   %f    spheritrans: %f\n", l1, l2, 2*std::min(l1, l2)/(l1+l2));



    printf("\n rot matrix\n");
    for (int i = 0; i<2; i++) {
        for (int j = 0; j<2; j++) {
            printf("%.3f ", rot[i][j]);
        }
        printf("\n");
    }
    // printf("%f   %f    spheritrans: %f\n", l1, l2, 2*std::min(l1, l2)/(l1+l2));



    std::vector<float> pAlign[2] = { {}, {} };
    printf("\n");
    for (int iPart = 0; iPart < 12; iPart++){
        pAlign[0].push_back(rot[0][0] * p[0][iPart] + rot[0][1] * p[1][iPart]);
        pAlign[1].push_back(rot[1][0] * p[0][iPart] + rot[1][1] * p[1][iPart]);

        printf("%f ", pAlign[0][iPart]);
    }
    printf("rotated\n");
    for (int iPart = 0; iPart < 12; iPart++){
        printf("%f ", pAlign[1][iPart]);
        


    }

    

}
