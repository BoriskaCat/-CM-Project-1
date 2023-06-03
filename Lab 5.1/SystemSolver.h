#pragma once
#include <stdio.h>
#include <list>

using namespace std;

typedef struct matrix matrix;


double* matrix_constget(const matrix* m, size_t i, size_t j);
double* matrix_get(matrix* m, size_t i, size_t j);
matrix* matrix_str_swap(matrix* m, size_t i1, size_t i2);
matrix* matrix_str_multi_cislo(matrix* m, double a, size_t i);
matrix* matrix_str_mainsum(matrix* m, size_t i1, size_t i2, double a);
matrix* matrix_alloc(size_t h, size_t w);
void matrix_free(matrix* m);
matrix* matrix_copy(const matrix* m);
matrix* matrix_assign(matrix* m1, const matrix* m2);
matrix* Gauss_urav(matrix* m, matrix* b, matrix* x);
matrix* matrix_vvod(matrix* m);
void matrix_vivod(matrix* m);
matrix* matrix_vvod_moment(matrix* m, list<double>& mom);
matrix* matrix_vvod_uzli(matrix* m, list<double>& uzli);
void from_matrix_to_list(matrix* m, list<double>& ko);
matrix* matrix_vvod_system_2(matrix* m, list<double>& mom_2);
matrix* matrix_vvod_stolb_2(matrix* m, list<double>& mom_2, size_t N);
