#include <iostream>
#include "Function.h"
#include <list>
#include <iterator>
#include <math.h>

using namespace std;

int main()
{
    cout << "Задание №5.3\nCQF Gauss\n";
    cout << "    ID: #10" << endl << endl;

    cout << "@@ Parameters: \n";
    cout << "\tФункция f(x) = sin(x)" << endl;
    cout << "\tФункция q(x) = |0.5 - x|" << endl << endl;

    cout << fixed;
    cout.precision(15);


    size_t k = 0, n = 5, m = 40;
    double a = 0, b = 1;

    do
    {
        cout << "\ta = ";
        cin >> a;
        cout << "\tb = ";
        cin >> b;
        cout << "Input n:\n";
        cout << "\tn = ";
        cin >> n;
        cout << "\tm = ";
        cin >> m;

        Resh(f, I, a, b, n, m);

        cout << "\t0 - start again" << endl;
        cout << "\t1 - the end" << endl;
        cin >> k;

    } while (k != 1);

}