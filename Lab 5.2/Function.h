#pragma once
#include <list>
#include <iostream>
using namespace std;

double f(double x);

double I(double x);

double f_mel(double x);

double f_7(double x);

double I_7(double x);

double f_9(double x);

double I_9(double x);

double f_13(double x);

double I_13(double x);

double f_15(double x);

double I_15(double x);

double P(size_t n, double x);   ///многочлен Лежандра

void print_list(const list<size_t>& for_n);     //для того чтобы напечатать список из n

size_t Tabulir(const double A, const double B, double N, list<double>& otrezki, size_t n); //табулирование (оттеление корней) возвращает кол-во отрезков

void print_otrezki(const list<double>& otrezki); //печатает отрезки на которых есть нечетные корни

double Secush(double a, double b, const double eps, size_t p, size_t n);

void Solutia(const double A, const double B, const double eps, double N, size_t p, list<double>& otrezki, size_t n, list<double>& korni);

void print_list(const list<double>& for_n);

void Koef(list<double>& uzli, list<double>& koef, size_t n);    //находит и печатает узлы кф (корень пробел узел)

double Leibniz(double (*I)(double x), double a, double b);    //вычисление собственного интеграла по фомуле Ньютона - Лейблица

double KF_Gaus(double (*f)(double x), double a, double b, list<double>& uzli, list<double>& koef);

void uzli_koef_Mel(list<double>& koef_2, list<double>& uzli_2, size_t n);

double KF_Meller(double (*f)(double x), list<double>& koef_2, list<double>& uzli_2, size_t n);