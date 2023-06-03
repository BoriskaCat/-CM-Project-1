#include <iostream>
#include <list>
#include <iterator>

using namespace std;

typedef struct matrix
{
    double* data;
    size_t h, w;
} matrix;

double* matrix_constget(const matrix* m, size_t i, size_t j)
{
    if (!m->data || i > m->h || j > m->w) return NULL;
    return m->data + m->w * i + j;
}

double* matrix_get(matrix* m, size_t i, size_t j)
{
    if (!m->data || i >= m->h || j >= m->w) return NULL;
    return m->data + m->w * i + j;
}

matrix* matrix_str_swap(matrix* m, size_t i1, size_t i2)
{
    if (!m || i1 >= m->h || i2 >= m->h) return NULL;
    double k;
    for (size_t j = 0; j < m->w; ++j)
    {
        k = *matrix_get(m, i1, j);
        *matrix_get(m, i1, j) = *matrix_get(m, i2, j);
        *matrix_get(m, i2, j) = k;
    }
    return m;

}

matrix* matrix_str_multi_cislo(matrix* m, double a, size_t i)
{
    if (!m || i >= m->h) return NULL;
    for (size_t j = 0; j < m->w; ++j)
    {
        *matrix_get(m, i, j) *= a;
    }
    return m;
}

matrix* matrix_str_mainsum(matrix* m, size_t i1, size_t i2, double a)
{
    if (!m || i1 >= m->h || i2 >= m->h) return NULL;
    for (size_t j = 0; j < m->w; ++j)
    {
        *matrix_get(m, i1, j) += *matrix_get(m, i2, j) * a;
    }
    return m;
}

matrix* matrix_alloc(size_t h, size_t w)
{
    matrix* m = new (matrix);
    if (!m) return NULL;
    m->data = new double[h * w];

    if (!m->data)
    {
        delete m;
        return NULL;
    }
    m->h = h;
    m->w = w;
    return m;
}

void matrix_free(matrix* m)
{
    delete m->data;
    delete m;
}

matrix* matrix_vvod(matrix* m)
{
    if (!m) return NULL;
    double* x;
    for (size_t i = 0; i < m->h; ++i)
    {
        for (size_t j = 0; j < m->w; ++j)
        {
            x = matrix_get(m, i, j);
            cin >> *x;
        }
    }

    return m;
}

matrix* matrix_vvod_moment(matrix* m, list<double>& mom)
{
    if (!m) return NULL;
    double* x;
    auto it = mom.cbegin();
    for (size_t i = 0; i < m->h; ++i)
    {
        for (size_t j = 0; j < m->w; ++j)
        {
            x = matrix_get(m, i, j);
            *x = *it;
        }
        ++it;
    }

    return m;
}

matrix* matrix_vvod_stolb_2(matrix* m, list<double>& mom_2, size_t N)
{
    if (!m) return NULL;
    double* x;
    auto it = mom_2.cbegin();

    for (size_t i = 0; i < N; ++i)
    {
        ++it;
    }

    for (size_t i = 0; i < m->h; ++i)
    {
        for (size_t j = 0; j < m->w; ++j)
        {
            x = matrix_get(m, i, j);
            *x = -(*it);
        }
        ++it;
    }

    return m;
}

matrix* matrix_vvod_uzli(matrix* m, list<double>& uzli)
{
    if (!m) return NULL;
    double* x;
    auto it = uzli.cbegin();
    for (size_t i = 0; i < m->h; ++i)
    {
        it = uzli.cbegin();
        for (size_t j = 0; j < m->w; ++j)
        {
            x = matrix_get(m, i, j);
            *x = pow(*it, i);
            ++it;
        }
    }

    return m;
}

matrix* matrix_vvod_system_2(matrix* m, list<double>& mom_2)
{
    if (!m) return NULL;
    double* x;
    auto it = mom_2.cbegin(), it_2 = mom_2.cbegin();;
    for (size_t i = 0; i < m->h; ++i)
    {
        it = it_2;

        for (size_t j = 0; j < m->w; ++j)
        {
            x = matrix_get(m, i, j);
            *x = *it;
            ++it;
        }
        ++it_2;
    }

    return m;
}

void from_matrix_to_list(matrix* m, list<double>& ko)
{
    double* x;
    for (size_t i = 0; i < m->h; ++i)
    {
        for (size_t j = 0; j < m->w; ++j)
        {
            x = matrix_get(m, i, j);
            ko.push_back(*x);
        }
    }
}

void matrix_vivod(matrix* m)
{
    double* x;
    for (size_t i = 0; i < m->h; ++i)
    {
        for (size_t j = 0; j < m->w; ++j)
        {
            x = matrix_get(m, i, j);
            printf("%lf ", *x);
        }
        printf("\n");
    }
}

matrix* matrix_copy(const matrix* m)
{
    if (!m) return NULL;
    matrix* m_copy;
    m_copy = matrix_alloc(m->h, m->w);
    if (!m_copy) return NULL;
    memcpy(m_copy->data, m->data, sizeof(double) * m->h * m->w);
    return m_copy;
}

matrix* matrix_assign(matrix* m1, const matrix* m2)
{
    if (!m1 || !m2 || m1->h != m2->h || m1->w != m2->w) return NULL;
    m1 = matrix_copy(m2);
    return m1;
}

matrix* Gauss_urav(matrix* m, matrix* b, matrix* x)
{
    if (!x || !m || !b) return NULL;
    size_t j = 0, i = 0, k = 0;                         
    double p;                                                  

    for (; k < m->w; ++k)
    {
        while (*matrix_constget(m, i, j) == 0 && i < m->h)     
        {
            ++i;

            if (*matrix_constget(m, i, j) == 0 && i == m->h - 1)
            {
                ++j;
                i = 0;
            }
        }
        if (i != m->h)                        
        {
            m = matrix_str_swap(m, k, i);     b = matrix_str_swap(b, k, i);
            for (i = k + 1; i < m->h; ++i)
            {
                if (*matrix_constget(m, i, j) != 0)
                {
                    p = *matrix_constget(m, i, j);
                    matrix_str_multi_cislo(m, -*matrix_constget(m, k, j), i);   matrix_str_multi_cislo(b, -*matrix_constget(m, k, j), i);
                    matrix_str_mainsum(m, i, k, p);         matrix_str_mainsum(b, i, k, p);
                }
            }
        }
        j = k + 1; i = k + 1;
    }


    if (*matrix_constget(m, m->h - 1, m->w - 1) == 0)
    {
        return NULL;
    }

    for (j = m->w - 1; j > 0; --j)
    {
        for (i = 0; i < j; ++i)          
        {
            p = *matrix_constget(m, i, j);
            matrix_str_multi_cislo(m, -*matrix_constget(m, j, j), i);  matrix_str_multi_cislo(b, -*matrix_constget(m, j, j), i);
            matrix_str_mainsum(m, i, j, p);   matrix_str_mainsum(b, i, j, p);
        }
    }


    for (i = 0; i < x->h; ++i)
    {
        p = *matrix_get(m, i, i);
        *matrix_get(m, i, i) /= p; *matrix_get(b, i, 0) /= p;
    }

    x = matrix_assign(x, b);

    return x;
}
