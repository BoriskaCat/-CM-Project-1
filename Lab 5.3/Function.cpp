
#include <list>
#include <iterator>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

double f(double x)
{
    return (sin(x) * abs(x - 0.5));
}

double I(double x)
{
    if ((0.5 - x) >= 0) return ((-sin(x) + (x - 0.5) * cos(x) + 0.479426) - 0.479426);
    else return(-(-sin(x) + (x - 0.5) * cos(x) + 0.479426) - 0.479426);
}

double P(size_t n, double x)   ///многочлен Лежандра
{
    if (n == 0) return 1;
    else if (n == 1) return x;
    else
    {
        return(((2.0 * n - 1) / n * P(n - 1, x)) * x - double((n - 1)) / n * P(n - 2, x));
    }
}

void print_list(const list<double>& for_n)
{
    auto it = for_n.cbegin();
    while (it != for_n.end())
    {
        cout << *it << " ";
        ++it;
    }
}

size_t Tabulir(const double A, const double B, double N, list<double>& otrezki, size_t n) //табулирование (оттеление корней) возвращает кол-во отрезков
{
    list<double>::iterator it = otrezki.begin();
    double h = (B - A) / N;
    double x = A, y = x + h;

    size_t k = 0;

    while (y <= B)
    {
        if ((P(n, x) * P(n, y)) < 0)
        {
            otrezki.push_back(x);
            otrezki.push_back(y);
            ++k;
        }
        x = y;
        y += h;
    }

    if (x < B && (P(n, x) * P(n, B)) < 0)
    {
        otrezki.push_back(x);
        otrezki.push_back(B);
        ++k;
    }

    return k;
}


double Secush(double a, double b, const double eps, size_t p, size_t n)
{
    if (P(n, a) * P(n, b) > 0) exit(1);

    double x1 = a, x2 = b, x3;
    size_t k = 0;

    x3 = x2 - (P(n, x2) * (x2 - x1) * p) / P(n, x2) - P(n, x1);

    while (abs(x3 - x2) > eps)
    {
        x1 = x2;
        x2 = x3;
        x3 = x2 - (P(n, x2) * (x2 - x1) * p) / (P(n, x2) - P(n, x1));
        ++k;
    }

    return x2;
}

void Solutia(const double A, const double B, const double eps, double N, size_t p, list<double>& otrezki, size_t n, list<double>& korni)
{
    size_t k = Tabulir(A, B, N, otrezki, n);

    auto it1 = otrezki.cbegin(), it2 = otrezki.cbegin();
    it2++;

    for (size_t i = 0; i < k - 1; i++)
    {
        korni.push_back(Secush(*it1, *it2, eps, p, n));
        it1++; it1++; it2++; it2++;
    }
    korni.push_back(Secush(*it1, *it2, eps, p, n));
}


void Koef(list<double>& uzli, list<double>& koef, size_t n)    //находит и печатает узлы кф (корень пробел узел)
{
    auto it = uzli.cbegin();
    double u = 0;

    cout << "\t\tKnots: " << "\t\t\t" << " Coeffs: " << endl;
    while (it != uzli.end())
    {
        u = (2 * (1 - pow(*it, 2))) / (pow(n, 2) * pow(P(n - 1, *it), 2));
        koef.push_back(u);
        cout << "\t\t" << *it << "\t" << u << endl;
        ++it;
    }
}

double Leibniz(double (*I)(double x), double a, double b)     //вычисление собственного интеграла по фомуле Ньютона - Лейблица
{
    return(I(b) - I(a));
}

void Resh(double (*f)(double x), double (*I)(double x), double a, double b, size_t n, size_t m)
{
    list<double> otrezki; //для хранения отрезков с решениями
    list<double> koef; //для хранения коэффициентов для Гаусса
    list<double> uzli; //для хранения узлов для Гаусса
    double I_kf, Integr, sum = 0, x, z_1, z_2;
    double h = (b - a) / m;

    Solutia(-1, 1, 1e-12, 1000, 1, otrezki, n, uzli);  //нахождение узлов
    Koef(uzli, koef, n);  //находим к коэф и печатаем коэф и узлы

    auto it_ko = koef.cbegin();
    auto it_uz = uzli.cbegin();

    for (size_t j = 0; j < m; ++j)
    {
        it_ko = koef.cbegin();
        it_uz = uzli.cbegin();
        z_1 = a + j * h;
        z_2 = a + (j + 1) * h;

        for (size_t i = 1; i <= n; ++i)
        {
            x = h / 2 * (*it_uz) + (z_1 + z_2) / 2;
            sum += (*it_ko) * f(x);
            ++it_ko;
            ++it_uz;
        }
    }

    sum *= h / 2;

    I_kf = sum;
    Integr = Leibniz(I, a, b);

    cout << "\n\t\tIKF = " << I_kf << endl;
    cout << "\t\tSol    = " << Integr << endl;

    printf("\t\t|IKF - Sol| = %e\n", abs(I_kf - Integr));

}