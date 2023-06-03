#include <iostream>
#include <list>
#include <iterator>
#include "SystemSolver.h"
#include "Function.h"

using namespace std;

int main()
{
    cout << "@@ #4.1. Approximate calculation of integrals using quadrature formulas (QF HADA)" << endl;
    cout << "@@ #ID: 10" << endl << endl;
    cout << "@@ Parameters:" << endl;
    cout << "   f (x) = cos(x) + 2x" << endl;
    cout << "   rho (x) = sqrt(1 - x)" << endl;
    cout << "   a = 0" << endl;
    cout << "   b = 1" << endl;

    cout << fixed;
    cout.precision(15);

    size_t k = 0, N;
    double a = 0, b = 1, eps;


    do
    {
        cout << endl;
        cout << "@@ Menu:" << endl;
        cout << "    0 - IKF" << endl;
        cout << "    1 - QF HADA" << endl;
        cout << "    2 - |IKF - Sol|" << endl;
        cout << "    3 - |QF HADA - Sol|" << endl;
        cin >> k;

        switch (k)
        {
        case 0:
        {
            list<double> uzli;

            cout << endl << "@@ IKF:" << endl;

            cout << "    N (<= 5) = ";
            cin >> N;
            cout << "    Enter values: ";
            vvod_uzli(uzli, N);

            Resh_1(a, b, uzli, N);

            uzli.clear();

        }
        break;

        case 1:
        {
            cout << endl << "@@ QF HADA:" << endl;

            cout << "    N (<= 5) = ";
            cin >> N;
            cout << "    eps = ";
            cin >> eps;

            Resh_2(a, b, N, eps);
        }
        break;

        case 2:
        {
            list<double> uzli;
            cout << "@@ |IKF - Sol|" << endl;
            cout << "    N = 5";
            cout << "    Enter values (5): ";
            vvod_uzli(uzli, 5);

            cout << endl;

            Resh_1(a, b, uzli, 5);

            uzli.clear();


        }

        case 3:
        {        
            cout << "@@ |QF HADA - Sol|" << endl;
            cout << "    N = 2";
            cout << "    eps = ";
            cin >> eps;

            cout << endl;
            Resh_2(a, b, 2, eps);
        }
        }

        cout << "@@ Menu:" << endl;
        cout << "    0 - start again" << endl;
        cout << "    1 - #the end " << endl;
        cin >> k;

    } while (k != 1);
}