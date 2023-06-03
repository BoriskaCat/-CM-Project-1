#pragma once
#include <list>

using namespace std;

double fq(double x);

double f(double x);

double q(double x);

double I(double x);

double Leibniz(double (*I)(double x), double a, double b);     //вычисление собственного интеграла по фомуле Ньютона - Лейблица

void vvod_uzli(list<double>& uzli, size_t N);   //ввод узлов

void print_list(list<double>& uzli);   //вывод узлов

double moment(double (*q)(double), double a, double b, size_t m, size_t k);   //подсчет момента с помощью скф средних прямоуг, где m - кол-во разбиений

void list_moment(list<double>& mom, size_t N, double a, double b, size_t m, double (*q)(double));  //заносит значения мооментов в список

void Koef(size_t N, list<double>& uzli, list<double>& mom, list<double>& ko);

double IKF(double (*f)(double), double a, double b, size_t N, list<double> ko, list<double> uzli);

void Resh_1(double a, double b, list<double> uzli, size_t N);

void Koef_2(size_t N, list<double>& mom_2, list<double>& a_i);

double Omega_n(double x, list<double>& a_i, size_t N);

size_t Tabulir(const double A, const double B, double N, list<double>& otrezki, list<double>& a_i, size_t n); //табулирование (оттеление корней) возвращает кол-во отрезков

void print_otrezki(const list<double>& otrezki); //печатает отрезки на которых есть нечетные корни

double Secush(double a, double b, const double eps, size_t p, list<double>& a_i, size_t n);

void Solutia(const double A, const double B, const double eps, double N, size_t p, list<double>& otrezki, list<double>& a_i, size_t n, list<double>& uzli);

void Resh_2(double a, double b, size_t N, double eps);