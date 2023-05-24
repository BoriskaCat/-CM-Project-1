#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>


using namespace std;


double F;
int m, n;
double a, b;
double eps, t;

vector <vector <double>> table;


double f (double x) 
{
    return cos(x) + 2 * x;
}

double Lagrange_r() 
{
    double result = 0;

    if (n == 0) result = table[0][0];

    for (int i = 0; i < n; ++i) 
    {
        double product = 1;

        for (int j = 0; j < n; ++j) 
        {
            if (i != j)
                product *= (F - table[j][1]) / (table[i][1] - table[j][1]);
        }

        result += product * table[i][0];
    }

    return result;
}

double Lagrange(double x) 
{
    double result = 0;

    if (n == 0) result = table[0][1];

    for (int i = 0; i < n; ++i) 
    {
        double product = 1;

        for (int j = 0; j < n; ++j) 
        {
            if (i != j)
                product *= (x - table[j][0]) / (table[i][0] - table[j][0]);
        }

        result += product * table[i][1];
    }

    return result;
}

void secants(double A, double B, double eps) 
{
    double x1,
           x2 = A,
           x3 = B;

    while (abs(x3 - x2) > eps) 
    {
        x1 = x2;
        x2 = x3;
        x3 = x2 - (Lagrange(x2) - F) * (x2 - x1) / ((Lagrange(x2) - F) - (Lagrange(x1) - F));
    }

    cout << "   Found value: " << x3 << endl;
    cout << "   |f(F_found) - F|: " << scientific << abs(f(x3) - F) << endl << endl;
}

int main()
{
    string ans = "y";
    bool isAdded = false;
    cout << "@@ #2. Inverse interpolation problem" << endl << endl;
    cout << "@@ #ID: 10" << endl << endl;
    cout << "@@ Parameters:" << endl;
    cout << "   f(x) = cos(x) + 2x" << endl;
    cout << "   m = ";
    cin >> m;
    cout << "   a (> -1) = ";
    cin >> a;
    while (a <= -1) 
    { 
        cout << "   Wrong value/ Try again" << endl << "   a (> -1) = ";
        cin >> a; 
    }
    cout << "   b (> -1) = ";
    cin >> b;
    while (b <= -1) 
    { 
        cout << "   Wrong value/ Try again" << endl << "   b (> -1) = ";
        cin >> b; 
    }

    table.assign(m + 1, vector<double>(2));
    double low_bound = a;
    double h = (b - a) / m;
    for (int k = 0; k <= m; ++k) 
    {
        table[k][0] = low_bound;
        table[k][1] = f(low_bound);
        low_bound += h;
    }

    cout << endl;
    cout << "@@ Table of values:" << endl;
    cout << "     x      |       f(x)    " << endl;
    cout.precision(10);
    for (int k = 0; k <= m; ++k) 
    {
        cout << "   ";
        printf("%.3f", table[k][0]);
        if (table[k][0] / 10 < 1) cout << " ";
        printf("   |   ");
        printf("%.11f", table[k][1]);
        cout << endl;
    }

    while (!ans.compare("y")) 
    {
        cout << endl;
        cout << endl << "@@ Interpolation point:" << endl << "   x = ";
        cin >> F;

        cout << endl << "   Degree of interpolation polynomial" << endl << "   Hint: n <= " << m << endl << "   n = ";
        cin >> n;
        while (n > m) 
        {
            cout << endl << "   Wrong value. Try again" << endl << "   n = ";
            cin >> n;
        }
        cout << endl;

        cout << "@@ Method #1" << endl;

        sort(table.begin(), table.end(), [](const auto& v1, const auto& v2)
            { return abs(v1[1] - F) < abs(v2[1] - F); });

        double L_res = Lagrange_r();
        cout << "   Lagrange value: " << L_res;
        cout << endl << "   |f(L_res) - F|: " << scientific << abs(f(L_res) - F) << endl;
        cout << endl;

        cout << "@@ Method #2" << endl;
        cout << "   Measurement accuracy (10^t , t < 0) " << endl << "   t = ";
        cin >> t;
        eps = pow(10, t);
        double roots[21][2];
        double h = (b - a) / 20;
        int counter = 0;
        double x1 = a;
        double x2 = a + h;
        double y2, y1 = f(x1);

        while (x2 <= b) 
        {
            y2 = f(x2);

            if ((F - y1) * (F - y2) <= 0) 
            {
                roots[counter][0] = x1;
                roots[counter][1] = x2;
                counter++;
                isAdded = true;
            }

            x1 = x2;
            x2 += h;
            y1 = y2;
        }

        double A, B;

        if (isAdded)
        {
            for (int i = 0; i < counter; i++)
            {
                A = roots[i][0];
                B = roots[i][1];
                setprecision(8);
                secants(A, B, eps);
            }
        }
        else
            cout << "   Unable to find the value" << endl << endl;

        cout << "?? Change the value of F, n, eps? [y / n]" << endl;
        cin >> ans;
    }

}