#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <Windows.h>

using namespace std;


double f(double x) {
    return exp(1.5 * x);
}


int main() {
    vector <vector <double>> table;
    double a, h;
    int m;
    string ans = "y";

    cout << "@@ #3.2. Derived functions" << endl;
    cout << "@@ #ID: 10" << endl << endl;
    cout << "@@ Parameters:" << endl;
    cout << "   f(x) = exp(1.5x)" << endl;

    while (!ans.compare("y")) {
        cout << "   m = ";
        cin >> m;
        cout << "   a = ";
        cin >> a;
        cout << "   h = ";
        cin >> h;

        table.assign(m + 1, vector<double>(8));
        for (int i = 0; i <= m; ++i) {
            table[i][0] = a + i * h; //x
            table[i][1] = f(table[i][0]); //f(x)
        }

        table[0][2] = (-3 * table[0][1] + 4 * table[1][1] - table[2][1]) / (2 * h);
        for (int i = 1; i < m; ++i) {
            table[i][2] = (table[i + 1][1] - table[i - 1][1]) / (2 * h); //f'(x)
        }
        table[m][2] = (table[m][1] - 4 * table[m - 1][1] + 3 * table[m - 2][1]) / (2 * h);

        for (int i = 0; i <= m; ++i) {
            table[i][3] = abs(1.5 * table[i][1] - table[i][2]); // val_1.1
        }

        for (int i = 0; i < m; ++i) {
            table[i][4] = abs(table[i][3] / table[i][2]); // val_1.2
        }

        table[0][5] = (table[0][1] - 2 * table[1][1] + table[2][1]) / (h * h);
        for (int i = 1; i < m; ++i) {
            table[i][5] = (table[i + 1][1] - 2 * table[i][1] + table[i - 1][1]) / (h * h); //f''(x)
        }
        table[m][5] = (table[m][1] - 2 * table[m - 1][1] + table[m - 2][1]) / (h * h);

        for (int i = 0; i <= m; ++i) {
            table[i][6] = abs(2.25 * table[i][1] - table[i][5]); // val_2.1
        }
        for (int i = 0; i <= m; ++i) {
            table[i][7] = abs(table[i][6] / table[i][5]); // val_2.2
        }

        cout << endl;
        cout << "@@ Table of values:" << endl;
        cout << "   val_1.1 = |f'_alg (x) - f'(x)|" << endl;
        cout << "   val_1.2 = |(f'_alg (x) - f'(x)) / f'_alg (x)|" << endl;
        cout << "   val_2.1 = |f''_alg (x) - f''(x)|" << endl;
        cout << "   val_2.2 = |(f''_alg (x) - f''(x)) / f''_alg (x)|" << endl << endl;

        cout.fill(' ');
        cout << setw(18) << "x";
        cout << setw(18) << "f(x)";
        cout << setw(18) << "f'_alg (x)";
        cout << setw(18) << "val_1.1";
        cout << setw(18) << "val_1.2";
        cout << setw(18) << "f''_alg (x)";
        cout << setw(18) << "val_2.1";
        cout << setw(18) << "val_2.2";

        cout << endl;
        for (int i = 0; i < m; ++i) {
            cout << fixed << setw(18) << table[i][0];
            cout << setw(18) << table[i][1];
            cout << setw(18) << table[i][2];
            cout << setw(18) << scientific << table[i][3];
            cout << setw(18) << table[i][4];
            cout << setw(18) << fixed << table[i][5];
            cout << setw(18) << scientific << table[i][6];
            cout << setw(18) << table[i][7];
            cout << endl;
        }
        cout << "?? Change the value of m, a, h? [y / n]" << endl;
        cin >> ans;
    }
}