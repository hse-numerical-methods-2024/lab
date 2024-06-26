#include <vector>
#include "main.cpp"

std::vector<std::vector<double>> diff(5, std::vector<double>(5));

void fill_diff() {
    double Stencil_3 = Differentiator<WhichD::XY, DiffMethod::Stencil3>(&F, 2.0, 2.0, h);
    double Stencil_3E = Differentiator<WhichD::XY, DiffMethod::Stencil3Extra>(&F, 2.0, 2.0, h);
    double Stencil_5 = Differentiator<WhichD::XY, DiffMethod::Stencil5>(&F, 2.0, 2.0, h);
    double Stencil_5E = Differentiator<WhichD::XY, DiffMethod::Stencil5Extra>(&F, 2.0, 2.0, h);
//    double AAD = Differentiator<WhichD::XY, DiffMethod::Stencil3>(&F, 2.0, 2.0, h);
    diff[0][0] = 0;
    diff[0][1] = Stencil_3 - Stencil_3E;
    diff[0][2] = Stencil_3 - Stencil_5;
    diff[0][3] = Stencil_3 - Stencil_5E;
//    diff[0][4] = Stencil_3 - AAD;
    diff[1][0] = Stencil_3E - Stencil_3;
    diff[1][1] = 0;
    diff[1][2] = Stencil_3E - Stencil_5;
    diff[1][3] = Stencil_3E - Stencil_5E;
//    diff[1][4] = Stencil_3E - AAD;
    diff[2][0] = Stencil_5 - Stencil_3;
    diff[2][1] = Stencil_5 - Stencil_3E;
    diff[2][2] = 0;
    diff[2][3] = Stencil_5 - Stencil_5E;
//    diff[2][4] = Stencil_5 - AAD;
    diff[3][0] = Stencil_5E - Stencil_3;
    diff[3][1] = Stencil_5E - Stencil_3E;
    diff[3][2] = Stencil_5E - Stencil_5;
    diff[3][3] = 0;
//    diff[3][4] = Stencil_5E - AAD;

//    diff[4][0] = AAD - Stencil_3;
//    diff[4][1] = AAD - Stencil_3E;
//    diff[4][2] = AAD - Stencil_5;
//    diff[4][3] = AAD  - Stencil_5E;
//    diff[4][4] = 0;

}

void test() {
    fill_diff();
    for (int i=0; i < 4; ++i) {
        for (int j=0; j < 4; ++j) {
            std::cout << fabs(diff[i][j]) << " ";
        }
        std::cout << "\n";
    }
}

int main() {
    test();
}
