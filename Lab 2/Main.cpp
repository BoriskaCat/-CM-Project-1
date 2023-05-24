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

vector <vector <double>> table;


double f(double x) {
    return cos(x) + 2 * x;
}

double Lagrange() 
{
    double result = 0;

    if (n == 0) result = table[0][1];

    for (int i = 0; i < n; ++i) 
    {
        double product = 1;

        for (int j = 0; j < n; ++j) 
        {
            if (i != j) 
                product *= (F - table[j][0]) / (table[i][0] - table[j][0]);
        }

        result += product * table[i][1];
    }

    return result;
}

double divDiff(int index) 
{
    double result = 0;

    for (int i = 0; i <= index; ++i) 
    {
        double tmp = 1;

        for (int j = 0; j <= index; ++j) 
        {
            if (i != j)
                tmp *= (table[i][0] - table[j][0]);
        }

        tmp = 1 / tmp;
        result += table[i][1] * tmp;  // f(x_i)/(x_i-x_j)
    }

    return result;
}


double Newton() 
{
    vector<double> coefs(n);

    for (int i = 0; i < n; ++i) 
    {
        if (i == 0)
            coefs[i] = table[i][1];
        else
            coefs[i] = divDiff(i);
    }

    double product = 1;
    double result = 0;

    for (int i = 0; i < n; i++) 
    {
        result += coefs[i] * product;
        product *= F - table[i][0];  // p(x) = a_0 + a_1(x-x_1) + a_2(x-x_1)(x-x_2)
    }

    if (n == 0) result = table[0][1];

    return result;
}

int main() 
{
    string ans = "y";
    cout << "@@ #2. Problem of algebraic interpolation." << endl << endl;
    cout << "@@ #ID: 10" << endl << endl;
    cout << "@@ Parameters:" << endl;
    cout << "   f(x) = cos(x) + 2x" << endl;
    cout << "   m = ";
    cin >> m;
    cout << "   a (> -1) = ";
    while (a <= -1) { cin >> a; }
    cout << "   b (> -1) = ";
    while (b <= -1) { cin >> b; }

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
    cout << "          x      |       f(x)    " << endl;
    cout.precision(10);
    for (int k = 0; k <= m; ++k) 
    {
        cout << "   ";
        printf("%.8f", table[k][0]);
        if (table[k][0] / 10 < 1) cout << " ";
        printf("   |   ");
        printf("%.8f", table[k][1]);
        cout << endl;
    }

    while (!ans.compare("y")) 
    {
        cout << endl << "@@ Interpolation point:" << endl << "   x = ";
        cin >> F;
        cout << endl <<  "@@ Degree of interpolation polynomial" << endl << "   Hint: n <= " << m << endl << "   n = ";
        cin >> n;
        while (n > m) 
        {
            cout << endl << "   Wrong value. Try again" << endl << "   n = ";
            cin >> n;
        }

        sort(table.begin(), table.end(), [](const auto& v1, const auto& v2)
            { return abs(v1[0] - F) < abs(v2[0] - F); });

        cout << endl;
        cout << "@@ Sorted table:" << endl;
        cout << "     x      |    f(x)    " << endl;
        for (int k = 0; k <= m; ++k) 
        {
            cout << "   ";
            printf("%.3f", table[k][0]);
            if (table[k][0] / 10 < 1) cout << " ";
            printf("   |   ");
            printf("%.3f", table[k][1]);
            cout << endl;
        }

        double L_res = Lagrange();
        double N_res = Newton();
        double f_res = f(F);
        cout << endl << "@@ Lagrange value: " << L_res;
        cout << endl << "  |L_res - f(x)|: " << abs(L_res - f_res);
        cout << endl;
        cout << endl << "@@ Newton value: " << N_res;
        cout << endl << "  |N_res - f(x)|: " << abs(N_res - f_res);
        cout << endl;
        cout << endl << "@@ True value: " << f_res << endl;
        cout << endl;

        cout << "?? Change the value of x or n? [y / n]" << endl;
        cin >> ans;
    }
    return 0;
}