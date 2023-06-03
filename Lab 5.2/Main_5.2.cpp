#include <iostream>
#include "Function.h"
#include <list>
#include <iterator>
#define _USE_MATH_DEFINES
#include <math.h>


using namespace std;


void Resh(double (*f)(double x), double (*I)(double x), double a, double b, size_t n)
{
    list<double> otrezki; //для хранения отрезков с решениями
    list<double> koef; //для хранения коэффициентов для Гаусса
    list<double> uzli; //для хранения узлов для Гаусса
    double I_kf, Integr;

    Solutia(-1, 1, 1e-12, 1000, 1, otrezki, n, uzli);  //нахождение узлов
    Koef(uzli, koef, n);  //находим к коэф и печатаем коэф и узлы

    I_kf = KF_Gaus(f, a, b, uzli, koef);
    Integr = Leibniz(I, a, b);
    cout << "\n\t\tIKF = " << I_kf << endl;
    cout << "\t\tI    = " << Integr << endl;

    printf("\t\t|IKF - Sol| = %e\n", abs(I_kf - Integr));
}

void Resh_2(double (*f)(double x), double a, double b, size_t n)
{
    list<double> koef_2; //для хранения коэффициентов для Меллера
    list<double> uzli_2; //для хранения узлов Меллера
    double I_kf, Integr;

    uzli_koef_Mel(koef_2, uzli_2, n);

    cout << "\t\tKnots: ";
    print_list(uzli_2);
    cout << endl;
    cout << "\t\tCoeffs: ";
    print_list(koef_2);

    I_kf = KF_Meller(f, koef_2, uzli_2, n);
    cout << "\n\t\tIKF = " << I_kf << endl;

}


int main()
{
    cout << "#5.2. QF Gauss. QF Meller\n";
    cout << "    ID: #10" << endl << endl;

    cout << "@@ Parameters: \n";
    cout << "\tf(x) = (x + 0.8) / sqrt( x^2 + 1.2)" << endl;
    cout << "\tf_7(x) = 5x^7 - 4x^3 + x^2 - 1" << endl;
    cout << "\tf_9(x) = x^9 + 4x^3 + 8" << endl;
    cout << "\tf_13(x) = -2x^13 -1" << endl;
    cout << "\tf_15(x) = 5x^15 +4x^7 -4x^6 - 11" << endl << endl;
    cout << "\tf_mel(x) = x^2*e^2x / sqrt(1 - x^2) " << endl;


    size_t k = 0, n1, n2, n3;
    double a, b;

    cout << fixed;
    cout.precision(15);

    do
    {
        cout << "\n\t0 - Gauss" << endl;
        cout << "\t1 - Meller" << endl;
        cin >> k;

        switch (k)
        {
        case 0:
        {
            cout << "\n@@ QF Gauss:\n\n";

            cout << "\ta = ";
            cin >> a;
            cout << "\tb = ";
            cin >> b;
            cout << "Input n:\n";
            cout << "\tn1 = ";
            cin >> n1;
            cout << "\tn2 = ";
            cin >> n2;
            cout << "\tn3 = ";
            cin >> n3;

            cout << "\n\tAccuracy for polynomials:\n\n";

            cout << "\tn = 4";
            cout << "\tf_7(x) = 5x^7 - 4x^3 + x^2 - 1\n";
            Resh(f_7, I_7, a, b, 4);

            cout << endl;
            cout << endl;

            cout << "\tn = 5";
            cout << "\tf_9(x) = x^9 + 4x^3 + 8\n";
            Resh(f_9, I_9, a, b, 5);

            cout << endl;
            cout << endl;

            cout << "\tn = 7";
            cout << "\tf_13(x) = -2x^13 -1\n";
            Resh(f_13, I_13, a, b, 7);

            cout << endl;
            cout << endl;

            cout << "\tComputing n1, n2, n3:\n\n";

            cout << "\tn1 = " << n1;
            cout << "\tf(x) = (x + 0.8) / sqrt( x^2 + 1.2)\n";
            Resh(f, I, a, b, n1);

            cout << endl;
            cout << endl;

            cout << "\tn2 = " << n2;
            cout << "\tf(x) = (x + 0.8) / sqrt( x^2 + 1.2)\n";
            Resh(f, I, a, b, n2);

            cout << endl;
            cout << endl;

            cout << "\tn3 = " << n3;
            cout << "\tf(x) = (x + 0.8) / sqrt( x^2 + 1.2)\n";
            Resh(f, I, a, b, n2);
        }
        break;

        case 1:
        {
            cout << "\nQF Meller:\n\n";

            cout << "\ta = ";
            cin >> a;
            cout << "\tb = ";
            cin >> b;
            cout << "Input n:\n";
            cout << "\tn1 = ";
            cin >> n1;
            cout << "\tn2 = ";
            cin >> n2;
            cout << "\tn3 = ";
            cin >> n3;

            cout << "\n\tn1 = " << n1 << endl;
            Resh_2(f_mel, a, b, n1);
            cout << "\tn2 = " << n2 << endl;
            Resh_2(f_mel, a, b, n2);
            cout << "\tn3 = " << n3 << endl;
            Resh_2(f_mel, a, b, n3);
        }
        break;
        }

        cout << "\t0 - start again" << endl;
        cout << "\t1 - the end" << endl;
        cin >> k;

    } while (k != 1);
}