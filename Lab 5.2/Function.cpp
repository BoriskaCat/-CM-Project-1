#include <list>
#include <iterator>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

double f(double x)
{
    return ((x + 0.8) / sqrt(x * x + 1.2));
}

double I(double x)
{
    return (4.0 / 5 * log(abs(sqrt(5 * x * x + 6) + sqrt(5.0) * x)) + sqrt((5 * x * x + 6) / 5.0));
}


double f_mel(double x)
{
    return (exp(2 * x) * x * x);
}

double f_7(double x)
{
    return (5 * pow(x, 7) - 4 * pow(x, 3) + pow(x, 2) - 1);
}

double I_7(double x)
{
    return(5.0 / 8 * pow(x, 8) - pow(x, 4) + 1.0 / 3 * pow(x, 3) - x);
}

double f_9(double x)
{
    return (pow(x, 9) + 4 * pow(x, 3) + pow(x, 1) + 8);
}

double I_9(double x)
{
    return(1.0 / 10 * pow(x, 10) + pow(x, 4) + 0.5 * pow(x, 2) + 8 * x);
}

double f_13(double x)
{
    return (-2 * pow(x, 13) - 1);
}

double I_13(double x)
{
    return(-1.0 / 7 * pow(x, 14) - x);
}

double f_15(double x)
{
    return (5 * pow(x, 15) + 4 * pow(x, 7) - 4 * pow(x, 6) - 11);
}

double I_15(double x)
{
    return(5.0 / 16 * pow(x, 16) + 0.5 * pow(x, 8) - 4.0 / 7 * pow(x, 7) - 11 * x);
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

void print_list(const list<size_t>& for_n)     //для того чтобы напечатать список из n
{
    auto it = for_n.cbegin();
    while (it != for_n.end())
    {
        cout << *it << " ";
        ++it;
    }
}

void print_list(const list<double>& for_n)     //для того чтобы напечатать список из n
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

void print_otrezki(const list<double>& otrezki) //печатает отрезки на которых есть нечетные корни
{
    auto it = otrezki.cbegin();
    while (it != otrezki.end())
    {
        cout << "[" << *it;
        ++it;
        cout << "; " << *it << "]" << endl;
        ++it;
    }

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

    /* cout << "\t\tКорень x = " << x2 << endl;
     cout << "\t\tКол-во шагов " << k << endl;
     cout << "\t\tНевязка " << abs(P(n, x2) - 0) << endl;*/

    return x2;
}

void Solutia(const double A, const double B, const double eps, double N, size_t p, list<double>& otrezki, size_t n, list<double>& korni)
{
    size_t k = Tabulir(A, B, N, otrezki, n);

    /*cout << endl;

    cout << "Отделение корней:" << endl;

    print_otrezki(otrezki);

    cout << endl;


    cout << "Количество отрезков, на которых есть корни нечетной кратности: " << k << endl;

    cout << endl;*/


    auto it1 = otrezki.cbegin(), it2 = otrezki.cbegin();
    it2++;

    for (size_t i = 0; i < k - 1; i++)
    {
        /*cout << "Отрезок: ";
        cout << "[" << *it1 << "; " << *it2 << "]" << endl;

        cout << "\tМетод секущих:\n";*/
        korni.push_back(Secush(*it1, *it2, eps, p, n));

        /*cout << endl;*/

        it1++; it1++; it2++; it2++;
    }

    /*cout << "Отрезок: ";
    cout << "[" << *it1 << "; " << *it2 << "]" << endl;

    cout << "\tМетод секущих:\n";*/
    korni.push_back(Secush(*it1, *it2, eps, p, n));

    /*cout << endl;*/
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

double KF_Gaus(double (*f)(double x), double a, double b, list<double>& uzli, list<double>& koef)
{
    double sum = 0, t;
    auto it_ko = koef.cbegin(), it_uz = uzli.cbegin();

    for (size_t i = 0; i < koef.size(); ++i)
    {
        t = (b - a) / 2 * (*it_uz) + (b + a) / 2;
        sum += (*it_ko) * f(t);

        ++it_uz;
        ++it_ko;
    }

    sum *= (b - a) / 2;

    return sum;
}

void uzli_koef_Mel(list<double>& koef_2, list<double>& uzli_2, size_t n)
{
    koef_2.push_back(M_PI / n);

    for (size_t i = 1; i <= n; ++i)
    {
        uzli_2.push_back(cos((2 * i - 1) * M_PI / (2 * n)));
    }
}


double KF_Meller(double (*f)(double x), list<double>& koef_2, list<double>& uzli_2, size_t n)
{
    double sum = 0;
    auto it_ko = koef_2.cbegin(), it_uz = uzli_2.cbegin();

    for (size_t i = 1; i <= n; ++i)
    {
        sum += f(*it_uz);
        ++it_uz;
    }

    sum *= *it_ko;

    return sum;
}