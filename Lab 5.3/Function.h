#pragma once
#include <list>
#include <iostream>
using namespace std;

double f(double x);

double I(double x);

double P(size_t n, double x);   ///многочлен Лежандра

void print_list(const list<double>& for_n);

size_t Tabulir(const double A, const double B, double N, list<double>& otrezki, size_t n); //табулирование (оттеление корней) возвращает кол-во отрезков

double Secush(double a, double b, const double eps, size_t p, size_t n);

void Solutia(const double A, const double B, const double eps, double N, size_t p, list<double>& otrezki, size_t n, list<double>& korni);

void Koef(list<double>& uzli, list<double>& koef, size_t n);   //находит и печатает узлы кф (корень пробел узел)

double Leibniz(double (*I)(double x), double a, double b);    //вычисление собственного интеграла по фомуле Ньютона - Лейблица

void Resh(double (*f)(double x), double (*I)(double x), double a, double b, size_t n, size_t m);
