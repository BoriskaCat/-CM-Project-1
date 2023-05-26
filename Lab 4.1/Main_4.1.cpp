#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <Windows.h>

using namespace std;


double f_integral(double A, double B) {
    return -pow(A,2) - sin(A) + pow(B, 2) + sin(B);
}

double f(double x) {
    return cos(x) + 2 * x;
}

double poly0_integral(double A, double B) {
    return 5 * (B - A);
}

double poly0(double x) {
    return 5;
}

double poly1_integral(double A, double B) {
    return 2.5 * (B * B - A * A);
}

double poly1(double x) {
    return 5 * x;
}

double poly2_integral(double A, double B) {
    return 1.0 / 6 * (-2 * pow(A, 3) - 15 * pow(A, 2) + pow(B, 2) * (2 * B + 15));
}

double poly2(double x) {
    return 5 * x + pow(x,2);
}

double poly3_integral(double A, double B) {
    return 1.0 / 4 * (-1 * pow(A, 4) - 10 * pow(A, 2) + pow(B, 4) + 10 * pow(B, 2));
}

double poly3(double x) {
    return 5 * x + pow(x, 3);
}


int main() {
    string ans = "y";
    double A, B;
    double h;
    double J, J_c;

    cout << "@@ #4.1. Calculation of integrals by quadrature formulas" << endl;
    cout << "@@ #ID: 10" << endl << endl;
    cout << "@@ Parameters:" << endl;
    cout << "   f_id (x) = cos(x) + 2x" << endl;
    while (!ans.compare("y")) {
        cout << "   A = ";
        cin >> A;
        cout << "   B = ";
        cin >> B;
        J = f_integral(A, B);

        cout << endl;
        cout << "@@ Results:" << endl;
        cout << "   f(x) = f_id (x)" << endl;
        cout << endl << "   Exact value of the integral J = " << setprecision(13) << J << endl << endl;

        cout << "   QF of the left rectangle: ";
        J_c = (B - A) * f(A);
        cout << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the right rectangle: ";
        J_c = (B - A) * f(B);
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the middle rectangle: ";
        J_c = (B - A) * f(0.5 * (A + B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF trapezoid: ";
        J_c = 0.5 * (B - A) * (f(A) + f(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   Simpson's QF: ";
        J_c = ((B - A) / 6.0) * (f(A) + 4 * f(0.5 * (A + B)) + f(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF 3/8: ";
        h = (B - A) / 3.0;
        J_c = (B - A) * (f(A) + 3 * f(A + h) + 3 * f(A + 2 * h) + f(B)) / 8.0;
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        //deg(poly) = 0
        J = poly0_integral(A, B);

        cout << endl;
        cout << "@@ Results:" << endl;
        cout << "   f(x) = 5" << endl;
        cout << endl << "   Exact value of the integral J = " << setprecision(13) << J << endl << endl;

        cout << "   QF of the left rectangle: ";
        J_c = (B - A) * poly0(A);
        cout << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the right rectangle: ";
        J_c = (B - A) * poly0(B);
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the middle rectangle: ";
        J_c = (B - A) * poly0(0.5 * (A + B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF trapezoid: ";
        J_c = 0.5 * (B - A) * (poly0(A) + poly0(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   Simpson's QF: ";
        J_c = ((B - A) / 6) * (poly0(A) + 4 * poly0(0.5 * (A + B)) + poly0(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF 3/8: ";
        h = (B - A) / 3.0;
        J_c = (B - A) * (poly0(A) + 3 * poly0(A + h) + 3 * poly0(A + 2 * h) + poly0(B)) / 8.0;
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        //deg(poly) = 1
        J = poly1_integral(A, B);

        cout << endl;
        cout << "@@ Results:" << endl;
        cout << "   f(x) = 5x" << endl;
        cout << endl << "   Exact value of the integral J = " << setprecision(13) << J << endl << endl;

        cout << "   QF of the left rectangle: ";
        J_c = (B - A) * poly1(A);
        cout << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the right rectangle: ";
        J_c = (B - A) * poly1(B);
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the middle rectangle: ";
        J_c = (B - A) * poly1(0.5 * (A + B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF trapezoid: ";
        J_c = 0.5 * (B - A) * (poly1(A) + poly1(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   Simpson's QF: ";
        J_c = ((B - A) / 6) * (poly1(A) + 4 * poly1(0.5 * (A + B)) + poly1(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF 3/8: ";
        h = (B - A) / 3.0;
        J_c = (B - A) * (poly1(A) + 3 * poly1(A + h) + 3 * poly1(A + 2 * h) + poly1(B)) / 8.0;
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        //deg(poly) = 2
        J = poly2_integral(A, B);

        cout << endl;
        cout << "@@ Results:" << endl;
        cout << "   f(x) = 5x + x^2" << endl;
        cout << endl << "   Exact value of the integral J = " << setprecision(13) << J << endl << endl;

        cout << "   QF of the left rectangle: ";
        J_c = (B - A) * poly2(A);
        cout << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the right rectangle: ";
        J_c = (B - A) * poly2(B);
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the middle rectangle: ";
        J_c = (B - A) * poly2(0.5 * (A + B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF trapezoid: ";
        J_c = 0.5 * (B - A) * (poly2(A) + poly2(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   Simpson's QF: ";
        J_c = ((B - A) / 6) * (poly2(A) + 4 * poly2(0.5 * (A + B)) + poly2(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF 3/8: ";
        h = (B - A) / 3.0;
        J_c = (B - A) * (poly2(A) + 3 * poly2(A + h) + 3 * poly2(A + 2 * h) + poly2(B)) / 8.0;
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        //deg(poly) = 3
        J = poly3_integral(A, B);

        cout << endl;
        cout << "@@ Results:" << endl;
        cout << "   f(x) = 5x + x^3" << endl;
        cout << endl << "   Exact value of the integral J = " << setprecision(13) << J << endl << endl;

        cout << "   QF of the left rectangle: ";
        J_c = (B - A) * poly3(A);
        cout << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the right rectangle: ";
        J_c = (B - A) * poly3(B);
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF of the middle rectangle: ";
        J_c = (B - A) * poly3(0.5 * (A + B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF trapezoid: ";
        J_c = 0.5 * (B - A) * (poly3(A) + poly3(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   Simpson's QF: ";
        J_c = ((B - A) / 6) * (poly3(A) + 4 * poly3(0.5 * (A + B)) + poly3(B));
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "   QF 3/8: ";
        h = (B - A) / 3;
        J_c = (B - A) * (poly3(A) + 3 * poly3(A + h) + 3 * poly3(A + 2 * h) + poly3(B)) / 8.0;
        cout << fixed << J_c << endl;
        cout << "   |J* - J|: " << scientific << abs(J - J_c) << endl << endl;

        cout << "?? Change the value of A, B? [y / n]" << endl;
        cin >> ans;
    }
}