#include <iostream>
#include <list>
#include <iterator>
#include <complex>
#include <math.h>
#include "SystemSolver.h"


using namespace std;

std::complex<double> i(0, 1.0);



double fq(double x)
{
    return (sin(x) * sqrt(1 - x));
}

double f(double x)
{
    return (sin(x));
}

double q(double x)
{
    return (sqrt(1 - x));
}

double Leibniz(double (*I)(double x), double a, double b)     //вычисление собственного интеграла по фомуле Ньютона - Лейблица
{
    return(I(b) - I(a));
}

void vvod_uzli(list<double>& uzli, size_t N)   //ввод узлов
{
    double x;

    for (size_t i = 1; i <= N; ++i)
    {
        cin >> x;
        uzli.push_back(x);
    }
}

void print_list(list<double>& list)   //вывод узлов
{
    auto it = list.cbegin();

    for (size_t i = 0; i < list.size(); ++i)
    {
        cout << *it << " ";
        ++it;
    }
}

double moment(double (*q)(double), double a, double b, size_t m, size_t k)   //подсчет момента с помощью скф средних прямоуг, где m - кол-во разбиений
{
    double sum = 0, h = (b - a) / m;
    double x;

    for (size_t j = 0; j < m; ++j)
    {
        x = a + (j + 0.5) * h;
        sum += q(x) * pow(x, k);
    }

    sum *= h;

    return(sum);
}


void list_moment(list<double>& mom, size_t N, double a, double b, size_t m, double (*q)(double))  //заносит значения мооментов в список
{
    double y;

    for (size_t i = 0; i < N; ++i)
    {
        y = moment(q, a, b, m, i);
        mom.push_back(y);
    }
}


void Koef(size_t N, list<double>& uzli, list<double>& mom, list<double>& ko)
{
    matrix* X, * M, * A;     //х - матр узлов, М - матр моменто, А - искомая матрица коэффициентов

    X = matrix_alloc(N, N);
    M = matrix_alloc(N, 1);
    A = matrix_alloc(N, 1);

    cout << "Матрица узлов:" << endl;
    X = matrix_vvod_uzli(X, uzli);
    matrix_vivod(X);

    cout << "\nСтолбец моментов:" << endl;
    M = matrix_vvod_moment(M, mom);
    matrix_vivod(M);

    cout << "\nСтолбец коэффициентов:" << endl;
    A = Gauss_urav(X, M, A);
    matrix_vivod(A);

    from_matrix_to_list(A, ko);

    matrix_free(X);
    matrix_free(M);
    matrix_free(A);
}


double IKF(double (*f)(double), double a, double b, size_t N, list<double> ko, list<double> uzli)
{
    double sum = 0;
    auto it_ko = ko.cbegin();
    auto it_uz = uzli.cbegin();

    for (size_t i = 1; i <= N; ++i)
    {
        sum += (*it_ko) * f(*it_uz);

        ++it_ko;
        ++it_uz;
    }

    return sum;
}

void Resh_1(double a, double b, list<double> uzli, size_t N)
{
    list<double> mom;
    list<double> ko;
    size_t m = 1000;
    double I_kf, Integr;

    list_moment(mom, N, a, b, m, q);
    cout << endl;

    cout << "Moments: ";
    print_list(mom);
    cout << endl;
    cout << endl;

    Koef(N, uzli, mom, ko);

    auto it_ko = ko.cbegin();
    auto it_uz = uzli.cbegin();
    cout << endl;
    cout << " Knots:" << "\t\t\t" << "Coeffs:" << endl;
    for (size_t i = 1; i <= N; ++i)
    {
        cout << *it_uz << "\t" << *it_ko << endl;
        ++it_uz; ++it_ko;
    }

    I_kf = IKF(f, a, b, N, ko, uzli);
    Integr = 0.25020169514301414;
    cout << "\n\tI_kf = " << I_kf << endl;
    cout << "\tI    = " << Integr << endl;

    printf("\t|I_kf - I| = %e\n", abs(I_kf - Integr));
}

void Koef_2(size_t N, list<double>& mom_2, list<double>& a_i)
{
    matrix* X, * M, * A;     //х - матр системы (3), М - столб свобод членов, А - искомый столб

    X = matrix_alloc(N, N);
    M = matrix_alloc(N, 1);
    A = matrix_alloc(N, 1);

    cout << "Матрица системы (3):" << endl;
    X = matrix_vvod_system_2(X, mom_2);
    matrix_vivod(X);

    cout << endl;

    cout << "Столбец свободных членов матрицы (3):" << endl;
    M = matrix_vvod_stolb_2(M, mom_2, N);
    matrix_vivod(M);

    cout << endl;

    cout << "Столбец решение матрицы (3):" << endl;
    A = Gauss_urav(X, M, A);
    matrix_vivod(A);

    from_matrix_to_list(A, a_i);


    matrix_free(X);
    matrix_free(M);
    matrix_free(A);
}

double Omega_n(double x, list<double>& a_i, size_t N)
{
    double sum = 0;
    auto it = a_i.cbegin();

    for (size_t i = 0; i < N; ++i)
    {
        sum += pow(x, i) * (*it);
        ++it;
    }

    sum += pow(x, N);

    return sum;
}



size_t Tabulir(const double A, const double B, double N, list<double>& otrezki, list<double>& a_i, size_t n) //табулирование (отделение корней) возвращает кол-во отрезков
{                                                                                                            //n - из задачи N число отрезочков
    list<double>::iterator it = otrezki.begin();
    double h = (B - A) / N;
    double x = A, y = x + h;

    size_t k = 0;

    while (y <= B)
    {
        if ((Omega_n(x, a_i, n) * Omega_n(y, a_i, n)) < 0)
        {
            otrezki.push_back(x);
            otrezki.push_back(y);
            ++k;
        }
        x = y;
        y += h;
    }

    if (x < B && (Omega_n(x, a_i, n) * Omega_n(B, a_i, n)) < 0)
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

double Secush(double a, double b, const double eps, size_t p, list<double>& a_i, size_t n)
{
    if (Omega_n(a, a_i, n) * Omega_n(b, a_i, n) > 0) exit(1);

    double x1 = a, x2 = b, x3;
    size_t k = 0;

    x3 = x2 - (Omega_n(x2, a_i, n) * (x2 - x1) * p) / (Omega_n(x2, a_i, n) - Omega_n(x1, a_i, n));

    while (abs(x3 - x2) > eps)
    {
        x1 = x2;
        x2 = x3;
        x3 = x2 - (Omega_n(x2, a_i, n) * (x2 - x1) * p) / (Omega_n(x2, a_i, n) - Omega_n(x1, a_i, n));
        ++k;
    }

    return x2;
}

void Solutia(const double A, const double B, const double eps, double N, size_t p, list<double>& otrezki, list<double>& a_i, size_t n, list<double>& uzli)
{
    size_t k = Tabulir(A, B, N, otrezki, a_i, n);

    auto it1 = otrezki.cbegin(), it2 = otrezki.cbegin();
    it2++;

    for (size_t i = 0; i < k - 1; i++)
    {
        uzli.push_back(Secush(*it1, *it2, eps, p, a_i, n));

        it1++; it1++; it2++; it2++;
    }

    uzli.push_back(Secush(*it1, *it2, eps, p, a_i, n));

    cout << endl;
}



void Resh_2(double a, double b, size_t N, double eps)
{
    list<double> mom_2;
    list<double> otrezki;
    size_t m = 1000;
    list<double> a_i;
    list<double> uzli;
    list<double> ko;
    double I_kfnst, Integr;

    list_moment(mom_2, 2 * N, a, b, m, q);

    cout << endl;
    cout << "    Moments: ";
    print_list(mom_2);

    cout << endl; cout << endl;

    Koef_2(N, mom_2, a_i);

    cout << endl;

    cout << "a_i: ";
    print_list(a_i);

    cout << endl;

    Solutia(a, b, eps, 1000, 1, otrezki, a_i, N, uzli);

    cout << "    knots (w_n(x) = 0):  ";
    print_list(uzli);

    cout << endl; cout << endl;

    Koef(N, uzli, mom_2, ko);
    cout << endl;

    auto it_ko = ko.cbegin();
    auto it_uz = uzli.cbegin();
    cout << " Knots:" << "\t\t\t" << "    Coeffs:" << endl;
    for (size_t i = 1; i <= N; ++i)
    {
        cout << *it_uz << "\t" << *it_ko << endl;
        ++it_uz; ++it_ko;
    }

    I_kfnst = IKF(f, a, b, N, ko, uzli);
    Integr = 0.25020169514301414;
    cout << "\n\t\QF HADA = " << I_kfnst << endl;
    cout << "\tI       = " << Integr << endl;

    printf("\t|QF HADA - Sol| = %e\n", abs(I_kfnst - Integr));

}