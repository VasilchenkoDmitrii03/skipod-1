/*
 * Последовательная / OpenMP-программа
 * Один и тот же файл: без флага OpenMP (-qsmp=omp / -fopenmp) собирается
 * последовательная версия, с флагом — многопоточная
 *
 * Запуск: ./prog M N [delta] [файл_решения]
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "problem.h"
#ifdef _OPENMP
#include <omp.h>
#endif

// Текущее время в секундах (источник разный для seq и OpenMP)
static double wtime(void)
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
#endif
}

int main(int argc, char **argv)
{
    int M, N, i, j, iter, max_iter, nthreads = 1;
    double h1, h2, h, eps, inv_eps2, hh, delta, t0, t1;
    double *w, *r, *z, *p, *Ap, *Q, *B, *Dinv;
    double zr, zr_new, pAp, pp, alpha, beta, diff = 0.0;
    size_t sz;
    const char *out = NULL;

    if (argc < 3) {
        fprintf(stderr, "usage: %s M N [delta] [out_file]\n", argv[0]);
        return 1;
    }
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    delta = argc > 3 ? atof(argv[3]) : 1e-6;
    if (argc > 4) {
        out = argv[4];
    }
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    // Шаги сетки и параметр фиктивной области eps = h^2
    h1 = (B1 - A1) / M;
    h2 = (B2 - A2) / N;
    h = h1 > h2 ? h1 : h2;
    eps = h * h;
    inv_eps2 = 1.0 / (eps * eps);
    hh = h1 * h2;
    max_iter = (M - 1) * (N - 1);

    // Все сеточные функции — одномерные массивы размера (M+1)*(N+1)
    sz = (size_t)(M + 1) * (N + 1);
    w    = calloc(sz, sizeof(double));
    r    = calloc(sz, sizeof(double));
    z    = calloc(sz, sizeof(double));
    p    = calloc(sz, sizeof(double));
    Ap   = calloc(sz, sizeof(double));
    Q    = calloc(sz, sizeof(double));
    B    = calloc(sz, sizeof(double));
    Dinv = calloc(sz, sizeof(double));
    if (!w || !r || !z || !p || !Ap || !Q || !B || !Dinv) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }

#define ID(i, j) ((size_t)(i) * (N + 1) + (j))

    t0 = wtime();

    /*
     * Коэффициенты Q_ij, F_ij и правая часть B = F + Q*v0
     * В граничные узлы кладём само v0 = phi, чтобы ниже применять оператор
     * ко всем внутренним узлам единой пятиточечной формулой
     */
#pragma omp parallel for private(j)
    for (i = 0; i <= M; i++) {
        for (j = 0; j <= N; j++) {
            double x = A1 + i * h1, y = A2 + j * h2;
            if (i == 0 || i == M || j == 0 || j == N) {
                w[ID(i, j)] = phi(x, y);
            } else {
                double S = cell_area(x - 0.5 * h1, x + 0.5 * h1,
                                     y - 0.5 * h2, y + 0.5 * h2);
                double q = (1.0 - S / hh) * inv_eps2;
                Q[ID(i, j)] = q;
                B[ID(i, j)] = S / hh * f_rhs(x, y) + q * phi(x, y);
                Dinv[ID(i, j)] = 1.0 / (2.0 / (h1 * h1) + 2.0 / (h2 * h2) + q);
            }
        }
    }

    // Старт метода: r0 = B - A*w0, z0 = D^{-1}*r0, p1 = z0
    zr = 0.0;
#pragma omp parallel for private(j) reduction(+:zr)
    for (i = 1; i < M; i++) {
        for (j = 1; j < N; j++) {
            size_t k = ID(i, j);
            double Aw = -(w[k + N + 1] - 2.0 * w[k] + w[k - N - 1]) / (h1 * h1)
                        - (w[k + 1] - 2.0 * w[k] + w[k - 1]) / (h2 * h2)
                        + Q[k] * w[k];
            r[k] = B[k] - Aw;
            z[k] = Dinv[k] * r[k];
            p[k] = z[k];
            zr += z[k] * r[k];
        }
    }
    zr *= hh;

    for (iter = 1; iter <= max_iter; iter++) {
        // Ap = A*p и знаменатель шага (Ap, p)_E
        pAp = 0.0;
#pragma omp parallel for private(j) reduction(+:pAp)
        for (i = 1; i < M; i++) {
            for (j = 1; j < N; j++) {
                size_t k = ID(i, j);
                Ap[k] = -(p[k + N + 1] - 2.0 * p[k] + p[k - N - 1]) / (h1 * h1)
                        - (p[k + 1] - 2.0 * p[k] + p[k - 1]) / (h2 * h2)
                        + Q[k] * p[k];
                pAp += Ap[k] * p[k];
            }
        }
        alpha = zr / (pAp * hh);

        // Шаг вдоль p: обновляем w, r, z, попутно считаем (p,p) и (z,r)
        pp = 0.0;
        zr_new = 0.0;
#pragma omp parallel for private(j) reduction(+:pp,zr_new)
        for (i = 1; i < M; i++) {
            for (j = 1; j < N; j++) {
                size_t k = ID(i, j);
                w[k] += alpha * p[k];
                r[k] -= alpha * Ap[k];
                z[k] = Dinv[k] * r[k];
                pp += p[k] * p[k];
                zr_new += z[k] * r[k];
            }
        }
        zr_new *= hh;

        // Остановка по норме приращения: ||w^{k+1}-w^k||_E = |alpha|*||p||_E
        diff = fabs(alpha) * sqrt(pp * hh);
        if (diff < delta) {
            break;
        }

        // Новое сопряжённое направление p = z + beta*p
        beta = zr_new / zr;
        zr = zr_new;
#pragma omp parallel for private(j)
        for (i = 1; i < M; i++) {
            for (j = 1; j < N; j++) {
                size_t k = ID(i, j);
                p[k] = z[k] + beta * p[k];
            }
        }
    }
    if (iter > max_iter) {
        iter = max_iter;
    }

    t1 = wtime();

    {
        // Контрольные величины — чтобы сверять версии между собой
        double sum = 0.0, wmax = -1e300;
#ifdef CHECK_EXACT
        double err = 0.0;
#endif
        for (i = 1; i < M; i++) {
            for (j = 1; j < N; j++) {
                sum += w[ID(i, j)] * hh;
                if (w[ID(i, j)] > wmax) {
                    wmax = w[ID(i, j)];
                }
#ifdef CHECK_EXACT
                {
                    double x = A1 + i * h1, y = A2 + j * h2;
                    if (in_domain(x, y) && fabs(w[ID(i, j)] - phi(x, y)) > err) {
                        err = fabs(w[ID(i, j)] - phi(x, y));
                    }
                }
#endif
            }
        }
        printf("M=%d N=%d threads=%d delta=%g iters=%d time=%.4f s "
               "last_diff=%.3e sum=%.10f max=%.10f w(0,0)=%.10f",
               M, N, nthreads, delta, iter, t1 - t0, diff, sum, wmax,
               w[ID(M / 2, N / 2)]);
#ifdef CHECK_EXACT
        printf(" err_C=%.3e", err);
#endif
        printf("\n");
    }

    // По запросу выгружаем решение: строки "x y w(x,y)" по всем узлам
    if (out) {
        FILE *fp = fopen(out, "w");
        if (fp) {
            for (i = 0; i <= M; i++) {
                for (j = 0; j <= N; j++) {
                    fprintf(fp, "%.6f %.6f %.10f\n",
                            A1 + i * h1, A2 + j * h2, w[ID(i, j)]);
                }
            }
            fclose(fp);
        }
    }

    free(w);
    free(r);
    free(z);
    free(p);
    free(Ap);
    free(Q);
    free(B);
    free(Dinv);
    return 0;
}
