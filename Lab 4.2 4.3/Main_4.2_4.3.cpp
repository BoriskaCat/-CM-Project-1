#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <Windows.h>

using namespace std;

double poly0_integral(double x) {
    return x;
}

double poly1(double x) {
    return 5.0 * x + 1;
}

double poly1_integral(double x) {
    return 2.5 * x * x + x;
}

double poly3(double x) {
    return 5.0 * x + pow(x, 3);
}

double poly3_integral(double x) {
    return pow(x, 4) / 4.0 + 2.5 * x * x;
}

double poly3_derivative(double x) {
    return 5 + 3.0 * pow(x, 2);
}

double f(double x) {
    return cos(x) + 2.0 * x;
}

double f_integral(double x) {
    return sin(x) + x * x;
}

double Runge(int l, int d, double J_hl, double J_h) {
    return (pow(l, d + 1) * J_hl - J_h) / (pow(l, d + 1) - 1);
}


int main() {
    vector<double> E_left;  // значения функции в левых узлах;
    vector<double> E_right;  // значения функции в правых узлах;
    vector<double> E_middle;  // значения функции в средних узлах;
    vector<double> E_der1;  // значения производных для вычисления теор.погрешности;
    vector<double> E_der2;  // значения производных для вычисления теор.погрешности;
    vector<double> E_der4;  // значения производных для вычисления теор.погрешности;

    string ans = "y";
    double A, B;
    double h, h2;
    int m;
    double J;
    vector<vector<double>> J_h;
    cout << "#4.2; #4.3 Approximate calculation of the integral from composite QF" << endl;
    cout << "@@ #ID: 10" << endl << endl;
    cout << "@@ Parameters:" << endl;
    cout << "   f_id (x) = cos(x) + 2x" << endl;
    cout << "   f_0 = 1" << endl;
    cout << "   f_1 = 5x + 1" << endl;
    cout << "   f_2 = 5x + x^2" << endl;
    cout << "   f_3 = 5x + x^3" << endl;

    while (!ans.compare("y")) {
        cout << endl;
        cout << "   A = "; cin >> A;
        cout << "   B = "; cin >> B;
        cout << "   m = "; cin >> m;
        h = (B - A) / m;
        h2 = h / 2;
        cout << "   h = " << h << endl;
        E_left.resize(m);
        E_right.resize(m);
        E_middle.resize(m);
        E_der1.resize(m);
        E_der2.resize(m);
        E_der4.resize(m);

        J_h.assign(4, vector<double>(5));


        cout << endl;
        cout << "@@ Results:" << endl;
        cout << "   f(x) = 1" << endl;
        J = poly0_integral(B) - poly0_integral(A);
        cout << endl << "   Exact value of the integral J(h) = " << setprecision(13) << J << endl << endl;

        cout << "   CQF of the left rectangle: ";
        J_h[0][0] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[0][0] += h;
        }
        cout << fixed << J_h[0][0] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[0][0]) << endl;
        cout << "   Theoretical error t = " << scientific << 0 << endl << endl;

        cout << "   CQF of the right rectangle: ";
        J_h[0][1] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[0][1] += h;
        }
        cout << fixed << J_h[0][1] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[0][1]) << endl;
        cout << "   Theoretical error t = " << scientific << 0 << endl << endl;

        cout << "   CQF of the middle rectangle: ";
        J_h[0][2] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[0][2] += h;
        }
        cout << fixed << J_h[0][2] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[0][2]) << endl;
        cout << "   Theoretical error t = " << 0 << endl << endl;

        cout << "   CQF trapezoid: ";
        J_h[0][3] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[0][3] += h;
        }
        cout << fixed << J_h[0][3] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[0][3]) << endl;
        cout << "   Theoretical error t = " << 0 << endl << endl;

        cout << "   Simpson's CQF: ";
        J_h[0][4] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[0][4] += h;
        }
        cout << fixed << J_h[0][4] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[0][4]) << endl;
        cout << "   Theoretical error t = " << scientific << 0 << endl << endl;


        cout << "@@ Results:" << endl;
        cout << "   f(x) = 5x + 1" << endl;
        J = poly1_integral(B) - poly1_integral(A);
        cout << endl << "   Exact value of the integral J(h) = " << fixed << J << endl << endl;

        for (int i = 0; i < m; ++i) {
            E_left[i] = poly1(A + i * h);
            E_right[i] = poly1(A + (i + 1) * h);
            E_middle[i] = poly1(A + (i + 0.5) * h);
        }

        cout << "   CQF of the left rectangle: ";
        J_h[1][0] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[1][0] += h * E_left[i];
        }
        cout << fixed << J_h[1][0] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[1][0]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 2.0 * 3 * (B - A) * h << endl << endl;

        cout << "   CQF of the right rectangle: ";
        J_h[1][1] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[1][1] += h * E_right[i];
        }
        cout << fixed << J_h[1][1] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[1][1]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 2.0 * 3 * (B - A) * h << endl << endl;

        cout << "   CQF of the middle rectangle: ";
        J_h[1][2] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[1][2] += h * E_middle[i];
        }
        cout << fixed << J_h[1][2] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[1][2]) << endl;
        cout << "   Theoretical error t = " << 0 << endl << endl;

        cout << "   CQF trapezoid: ";
        J_h[1][3] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[1][3] += h * 0.5 * (E_left[i] + E_right[i]);
        }
        cout << fixed << J_h[1][3] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[1][3]) << endl;
        cout << "   Theoretical error t = " << 0 << endl << endl;

        cout << "   Simpson's CQF: ";
        J_h[1][4] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[1][4] += h / 6.0 * (E_left[i] + 4 * E_middle[i] + E_right[i]);
        }
        cout << fixed << J_h[1][4] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[1][4]) << endl;
        cout << "   Theoretical error t = " << scientific << 0 << endl << endl;


        cout << "@@ Results:" << endl;
        cout << "   f(x) = 5x + x^3" << endl;
        J = poly3_integral(B) - poly3_integral(A);
        cout << endl << "   Exact value of the integral J(h) = " << fixed << J << endl << endl;

        for (int i = 0; i < m; ++i) {
            E_left[i] = poly3(A + i * h);
            E_right[i] = poly3(A + (i + 1) * h);
            E_middle[i] = poly3(A + (i + 0.5) * h);
            E_der1[i] = poly3_derivative(A + i * h);
            E_der2[i] = 6 * (A + i * h);
        }

        cout << "   CQF of the left rectangle: ";
        J_h[2][0] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[2][0] += h * E_left[i];
        }
        cout << fixed << J_h[2][0] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[2][0]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 2.0 * *max_element(E_der1.begin(), E_der1.end()) * (B - A) * h << endl << endl;

        cout << "   CQF of the right rectangle: ";
        J_h[2][1] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[2][1] += h * E_right[i];
        }

        cout << fixed << J_h[2][1] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[2][1]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 2.0 * *max_element(E_der1.begin(), E_der1.end()) * (B - A) * h << endl << endl;

        cout << "   CQF of the middle rectangle: ";
        J_h[2][2] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[2][2] += h * E_middle[i];
        }
        cout << fixed << J_h[2][2] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[2][2]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 24.0 * *max_element(E_der2.begin(), E_der2.end()) * (B - A) * pow(h, 2) << endl << endl;

        cout << "   CQF trapezoid: ";
        J_h[2][3] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[2][3] += h * 0.5 * (E_left[i] + E_right[i]);
        }
        cout << fixed << J_h[2][3] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[2][3]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 12.0 * *max_element(E_der2.begin(), E_der2.end()) * (B - A) * pow(h, 2) << endl << endl;

        cout << "   Simpson's CQF: ";
        J_h[2][4] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[2][4] += h / 6.0 * (E_left[i] + 4 * E_middle[i] + E_right[i]);
        }
        cout << fixed << J_h[2][4] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[2][4]) << endl;
        cout << "   Theoretical error t = " << 0 << endl << endl;


        cout << "@@ Results:" << endl;
        cout << "   f(x) = cos(x) + 5x" << endl;
        J = f_integral(B) - f_integral(A);
        cout << endl << "   Exact value of the integral J(h) = " << fixed << J << endl << endl;

        for (int i = 0; i < m; ++i) {
            E_left[i] = f(A + i * h);
            E_right[i] = f(A + (i + 1) * h);
            E_middle[i] = f(A + (i + 0.5) * h);
            E_der1[i] = -sin(A + i * h) + 2;
            E_der2[i] = -cos(A + i * h);
            E_der4[i] = -E_der2[i];
        }

        cout << "   CQF of the left rectangle: ";
        J_h[3][0] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[3][0] += h * E_left[i];
        }
        cout << fixed << J_h[3][0] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[3][0]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 2.0 * *max_element(E_der1.begin(), E_der1.end()) * (B - A) * h << endl << endl;

        cout << "   CQF of the right rectangle: ";
        J_h[3][1] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[3][1] += h * E_right[i];
        }
        cout << fixed << J_h[3][1] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[3][1]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 2.0 * *max_element(E_der1.begin(), E_der1.end()) * (B - A) * h << endl << endl;

        cout << "   CQF of the right rectangle: ";
        J_h[3][2] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[3][2] += h * E_middle[i];
        }
        cout << fixed << J_h[3][2] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[3][2]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 24.0 * *max_element(E_der2.begin(), E_der2.end()) * (B - A) * pow(h, 2) << endl << endl;

        cout << "   CQF trapezoid: ";
        J_h[3][3] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[3][3] += h * 0.5 * (E_left[i] + E_right[i]);
        }
        cout << fixed << J_h[3][3] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[3][3]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 12.0 * *max_element(E_der2.begin(), E_der2.end()) * (B - A) * pow(h, 2) << endl << endl;

        cout << "   Simpson's CQF: ";
        J_h[3][4] = 0;
        for (int i = 0; i < m; ++i) {
            J_h[3][4] += h / 6.0 * (E_left[i] + 4 * E_middle[i] + E_right[i]);
        }
        cout << fixed << J_h[3][4] << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_h[3][4]) << endl;
        cout << "   Theoretical error t = " << scientific << 1 / 2880.0 * *max_element(E_der4.begin(), E_der4.end()) * (B - A) * pow(h, 4) << endl << endl << endl;


        cout << "?? Change the value of A, B? [y / n]" << endl;
        cin >> ans;
    }
    cout << "%%%%%-----------------------------------------------------------------------------------------------------------%%%%%" << endl;
    cout << "@@ Increase the rank of the partition" << endl;
    cout << "   Hint: The rank will be increased 'l' times" << endl;
    E_left.clear();
    E_middle.clear();
    E_right.clear();
    E_der1.clear();
    E_der2.clear();
    E_der4.clear();

    int l;
    double J_hl, J_r;
    cout << "   l = "; cin >> l;
    h /= l;
    cout << "   h = " << h << endl;
    cout << endl;


    E_left.resize(m * l);
    E_right.resize(m * l);
    E_middle.resize(m * l);
    E_der1.resize(m * l);
    E_der2.resize(m * l);
    E_der4.resize(m * l);

    cout << "@@ Results:" << endl;
    cout << "   f(x) = 1" << endl;
    J = poly0_integral(B) - poly1_integral(A);
    cout << endl << "   Exact value of the integral J(h/l) = " << setprecision(13) << J << endl << endl;

    cout << "   mCQF of the left rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h;
    }
    J_r = Runge(l, 0, J_hl, J_h[0][0]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the left rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h;
    }
    J_r = Runge(l, 0, J_hl, J_h[0][1]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the right rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h;
    }
    J_r = Runge(l, 1, J_hl, J_h[0][2]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF trapezoid: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h;
    }
    J_r = Runge(l, 1, J_hl, J_h[0][3]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   Simpson's mCQF: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h;
    }
    J_r = Runge(l, 3, J_hl, J_h[0][4]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl << endl;


    cout << "@@ Results:" << endl;
    cout << "   f(x) = 5x + 1" << endl;
    J = poly1_integral(B) - poly1_integral(A);
    cout << endl << "   Exact value of the integral J(h) = " << fixed << J << endl << endl;

    for (int i = 0; i < m * l; ++i) {
        E_left[i] = poly1(A + i * h);
        E_right[i] = poly1(A + (i + 1) * h);
        E_middle[i] = poly1(A + (i + 0.5) * h);
    }

    cout << "   mCQF of the left rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_left[i];
    }
    J_r = Runge(l, 0, J_hl, J_h[1][0]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the right rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_right[i];
    }
    J_r = Runge(l, 0, J_hl, J_h[1][1]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the middle rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_middle[i];
    }
    J_r = Runge(l, 1, J_hl, J_h[1][2]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF trapezoid: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * 0.5 * (E_left[i] + E_right[i]);
    }
    J_r = Runge(l, 1, J_hl, J_h[1][3]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   Simpson's mCQF: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h / 6.0 * (E_left[i] + 4 * E_middle[i] + E_right[i]);
    }
    J_r = Runge(l, 3, J_hl, J_h[1][4]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl << endl;


    cout << "@@ Results:" << endl;
    cout << "   f(x) = 5x + x^3" << endl;
    J = poly3_integral(B) - poly3_integral(A);
    cout << endl << "   Exact value of the integral J(h) = " << fixed << J << endl << endl;

    for (int i = 0; i < m * l; ++i) {
        E_left[i] = poly3(A + i * h);
        E_right[i] = poly3(A + (i + 1) * h);
        E_middle[i] = poly3(A + (i + 0.5) * h);
        E_der1[i] = poly3_derivative(A + i * h);
        E_der2[i] = 6 * (A + i * h);
    }

    cout << "   mCQF of the left rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_left[i];
    }
    J_r = Runge(l, 0, J_hl, J_h[2][0]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the right rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_right[i];
    }
    J_r = Runge(l, 0, J_hl, J_h[2][1]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the middle rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_middle[i];
    }
    J_r = Runge(l, 1, J_hl, J_h[2][2]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF trapezoid: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * 0.5 * (E_left[i] + E_right[i]);
    }
    J_r = Runge(l, 1, J_hl, J_h[2][3]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   Simpson's mCQF: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h / 6.0 * (E_left[i] + 4 * E_middle[i] + E_right[i]);
    }
    J_r = Runge(l, 3, J_hl, J_h[2][4]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl << endl;


    cout << "@@ Results:" << endl;
    cout << "   f(x) = cos(x) + 2x" << endl;
    J = f_integral(B) - f_integral(A);
    cout << endl << "   Exact value of the integral J(h) = " << fixed << J << endl << endl;

    for (int i = 0; i < m * l; ++i) {
        E_left[i] = f(A + i * h);
        E_right[i] = f(A + (i + 1) * h);
        E_middle[i] = f(A + (i + 0.5) * h);
        E_der1[i] = -sin(A + i * h) + 2;
        E_der2[i] = -cos(A + i * h);
        E_der4[i] = -E_der2[i];
    }

    cout << "   mCQF of the left rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_left[i];
    }
    J_r = Runge(l, 0, J_hl, J_h[3][0]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the right rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_right[i];
    }
    J_r = Runge(l, 0, J_hl, J_h[3][1]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF of the middle rectangle: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * E_middle[i];
    }
    J_r = Runge(l, 1, J_hl, J_h[3][2]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   mCQF trapezoid: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h * 0.5 * (E_left[i] + E_right[i]);
    }
    J_r = Runge(l, 1, J_hl, J_h[3][3]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl;

    cout << "   Simpson's mCQF: " << endl;
    J_hl = 0;
    for (int i = 0; i < m * l; ++i) {
        J_hl += h / 6.0 * (E_left[i] + 4 * E_middle[i] + E_right[i]);
    }
    J_r = Runge(l, 3, J_hl, J_h[3][4]);
    cout << "   J(h/l) = " << fixed << J_hl << endl;
    cout << "   |J(h/l) - J| = " << scientific << abs(J - J_hl) << endl;
    cout << "   J*_r = " << scientific << J_r << endl;
    cout << "   |J*_r - J| = " << abs(J - J_r) << endl << endl << endl;
}