//
// Created by vistart on 2022/1/25.
//

#include <iostream>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <vector>
#include <map>

#define CONE_TOL (1e-9)
#define CONE_THRESH (1e-8)
#define EXP_CONE_MAX_ITERS (100)
#define BOX_CONE_MAX_ITERS (25)
#define POW_CONE_MAX_ITERS (20)

/* Box cone limits (+ or -) taken to be INF */
#define MAX_BOX_VAL (1e15)

#ifndef SCS
#define SCS(x) _scs_##x
#endif

typedef int    scs_int;
typedef double  scs_float;
#define scs_free free
#define scs_malloc malloc
#define scs_calloc calloc
#define scs_realloc realloc

#define scs_printf printf

#define POWF pow
#define SQRTF sqrt

#define SCS_NULL 0 /* NULL type */

/** This defines the data matrices which should be supplied in compressed
 *  sparse column format with zero based indexing.
 */
typedef struct {
    /** Matrix values, size: number of non-zeros. */
    scs_float* x;
    /** Matrix row indices, size: number of non-zeros. */
    scs_int* i;
    /** Matrix column pointers, size: `n+1`. */
    scs_int* p;
    /** Number of rows. */
    scs_int m;
    /** Number of columns. */
    scs_int n;
} ScsMatrix;

/** Struct containing all settings. */
typedef struct {
    /** Whether to heuristically rescale the data before solve. */
    scs_int normalize;
    /** Initial dual scaling factor (may be updated if adaptive_scale is on). */
    scs_float scale;
    /** Whether to adaptively update `scale`. */
    scs_int adaptive_scale;
    /** Primal constraint scaling factor. */
    scs_float rho_x;
    /** Maximum iterations to take. */
    scs_int max_iters;
    /** Absolute convergence tolerance. */
    scs_float eps_abs;
    /** Relative convergence tolerance. */
    scs_float eps_rel;
    /** Infeasible convergence tolerance. */
    scs_float eps_infeas;
    /** Douglas-Rachford relaxation parameter. */
    scs_float alpha;
    /** Time limit in secs (can be fractional). */
    scs_float time_limit_secs;
    /** Whether to log progress to stdout. */
    scs_int verbose;
    /** Whether to use warm start (put initial guess in ScsSolution struct). */
    scs_int warm_start;
    /** Memory for acceleration. */
    scs_int acceleration_lookback;
    /** Interval to apply acceleration. */
    scs_int acceleration_interval;
    /** String, if set will dump raw prob data to this file. */
    const char* write_data_filename;
    /** String, if set will log data to this csv file (makes SCS very slow). */
    const char* log_csv_filename;
} ScsSettings;

/** Struct containing problem data. */
typedef struct {
    /** A has `m` rows. */
    scs_int m;
    /** A has `n` cols, P has `n` cols and `n` rows. */
    scs_int n;
    /** A is supplied in CSC format (size `m` x `n`). */
    ScsMatrix* A;
    /** P is supplied in CSC format (size `n` x `n`), must be symmetric positive
     * semidefinite. Only pass in the upper triangular entries. If `P = 0` then
     * set `P = SCS_NULL`. */
    ScsMatrix* P;
    /** Dense array for b (size `m`). */
    scs_float* b;
    /** Dense array for c (size `n`). */
    scs_float* c;
} ScsData;

/** Cone data. Rows of data matrix `A` must be specified in this exact order. */
typedef struct {
    /** Number of linear equality constraints (primal zero, dual free). */
    scs_int z;
    /** Number of positive orthant cones. */
    scs_int l;
    /** Upper box values, `len(bu) = len(bl) = max(bsize-1, 0)`. */
    scs_float* bu;
    /** Lower box values, `len(bu) = len(bl) = max(bsize-1, 0)`. */
    scs_float* bl;
    /** Total length of box cone (includes scale `t`). */
    scs_int bsize;
    /** Array of second-order cone constraints, `len(q) = qsize`. */
    scs_int* q;
    /** Length of second-order cone array `q`. */
    scs_int qsize;
    /** Array of semidefinite cone constraints, `len(s) = ssize`. */
    scs_int* s;
    /** Length of semidefinite constraints array `s`. */
    scs_int ssize;
    /** Number of primal exponential cone triples. */
    scs_int ep;
    /** Number of dual exponential cone triples. */
    scs_int ed;
    /** Array of power cone params, must be in `[-1, 1]`, negative values are
     * interpreted as specifying the dual cone, `len(p) = psize ` */
    scs_float* p;
    /** Number of (primal and dual) power cone triples. */
    scs_int psize;
} ScsCone;

/** Contains primal-dual solution arrays or a certificate of infeasibility.
 *  Check the exit flag to determine whether this contains a solution or a
 *  certificate. If when passed into SCS the members `x`, `y`, `s` are
 *  NULL then SCS will allocate memory for them which should be managed
 *  by the user to prevent memory leaks.
 */
typedef struct {
    /** Primal variable. */
    scs_float* x;
    /** Dual variable. */
    scs_float* y;
    /** Slack variable. */
    scs_float* s;
} ScsSolution;

/** Contains information about the solve run at termination. */
typedef struct {
    /** Number of iterations taken. */
    scs_int iter;
    /** Status string, e.g. 'solved'. */
    char status[128];
    /** Linear system solver used. */
    char lin_sys_solver[128];
    /** Status as scs_int, defined in glbopts.h. */
    scs_int status_val;
    /** Number of updates to scale. */
    scs_int scale_updates;
    /** Primal objective. */
    scs_float pobj;
    /** Dual objective. */
    scs_float dobj;
    /** Primal equality residual. */
    scs_float res_pri;
    /** Dual equality residual. */
    scs_float res_dual;
    /** Duality gap. */
    scs_float gap;
    /** Infeasibility cert residual. */
    scs_float res_infeas;
    /** Unbounded cert residual. */
    scs_float res_unbdd_a;
    /** Unbounded cert residual. */
    scs_float res_unbdd_p;
    /** Time taken for setup phase (milliseconds). */
    scs_float setup_time;
    /** Time taken for solve phase (milliseconds). */
    scs_float solve_time;
    /** Final scale parameter. */
    scs_float scale;
    /** Complementary slackness. */
    scs_float comp_slack;
    /** Number of rejected AA steps. */
    scs_int rejected_accel_steps;
    /** Number of accepted AA steps. */
    scs_int accepted_accel_steps;
    /** Total time (milliseconds) spent in the linear system solver. */
    scs_float lin_sys_time;
    /** Total time (milliseconds) spent in the cone projection. */
    scs_float cone_time;
    /** Total time (milliseconds) spent in the acceleration routine. */
    scs_float accel_time;
} ScsInfo;

/** Contains normalization variables. */
typedef struct {
    scs_float* D, * E; /* for normalization */
    scs_int m;        /* Length of D */
    scs_int n;        /* Length of E */
    scs_float primal_scale, dual_scale;
} ScsScaling;

/* private data to help cone projection step */
struct SCS_CONE_WORK {
    /*
     * cone_boundaries will contain array of indices of rows of A corresponding to
     * cone boundaries, boundaries[0] is starting index for cones of size larger
     * than 1
     */
    ScsCone* k; /* original cone information */
    scs_int* cone_boundaries;
    scs_int cone_boundaries_len;
    scs_int scaled_cones; /* boolean, whether the cones have been scaled */
    scs_float* s;         /* used for Moreau decomposition in projection */
    scs_int m;            /* total length of cone */
    /* box cone quantities */
    scs_float box_t_warm_start;
#ifdef USE_LAPACK
    /* workspace for eigenvector decompositions: */
    scs_float* Xs, * Z, * e, * work;
    blas_int lwork;
#endif
};

/** Struct containing cone projection workspace. Implemented by cones. */
typedef struct SCS_CONE_WORK ScsConeWork;

#define _MAX_RAND_VAL (1073741823) /* 2^30 - 1 */
/************ see the book for explanations and caveats! *******************/
/************ in particular, you need two's complement arithmetic **********/

#define KK 100                                  /* the long lag */
#define LL 37                                   /* the short lag */
#define MM (1L << 30)                           /* the modulus */
#define mod_diff(x, y) (((x) - (y)) & (MM - 1)) /* subtraction mod MM */

long ran_x[KK]; /* the generator state */

#ifdef __STDC__
void ran_array(long aa[], int n)
#else
void ran_array(aa, n) /* put n new random numbers in aa */
long* aa;         /* destination */
int n;                /* array length (must be at least KK) */
#endif
{
    register int i, j;
    for (j = 0; j < KK; j++)
        aa[j] = ran_x[j];
    for (; j < n; j++)
        aa[j] = mod_diff(aa[j - KK], aa[j - LL]);
    for (i = 0; i < LL; i++, j++)
        ran_x[i] = mod_diff(aa[j - KK], aa[j - LL]);
    for (; i < KK; i++, j++)
        ran_x[i] = mod_diff(aa[j - KK], ran_x[i - LL]);
}

/* the following routines are from exercise 3.6--15 */
/* after calling ran_start, get new randoms by, e.g., "x=ran_arr_next()" */

#define QUALITY 1009 /* recommended quality level for high-res use */
long ran_arr_buf[QUALITY];
long ran_arr_dummy = -1, ran_arr_started = -1;
long* ran_arr_ptr = &ran_arr_dummy; /* the next random number, or -1 */

#define TT 70             /* guaranteed separation between streams */
#define is_odd(x) ((x)&1) /* units bit of x */

#ifdef __STDC__
void ran_start(long seed)
#else
void ran_start(seed)  /* do this before using ran_array */
long seed;        /* selector for different streams */
#endif
{
    register int t, j;
    long x[KK + KK - 1]; /* the preparation buffer */
    register long ss = (seed + 2) & (MM - 2);
    for (j = 0; j < KK; j++) {
        x[j] = ss; /* bootstrap the buffer */
        ss <<= 1;
        if (ss >= MM)
            ss -= MM - 2; /* cyclic shift 29 bits */
    }
    x[1]++; /* make x[1] (and only x[1]) odd */
    for (ss = seed & (MM - 1), t = TT - 1; t;) {
        for (j = KK - 1; j > 0; j--)
            x[j + j] = x[j], x[j + j - 1] = 0; /* "square" */
        for (j = KK + KK - 2; j >= KK; j--)
            x[j - (KK - LL)] = mod_diff(x[j - (KK - LL)], x[j]),
            x[j - KK] = mod_diff(x[j - KK], x[j]);
        if (is_odd(ss)) { /* "multiply by z" */
            for (j = KK; j > 0; j--)
                x[j] = x[j - 1];
            x[0] = x[KK]; /* shift the buffer cyclically */
            x[LL] = mod_diff(x[LL], x[KK]);
        }
        if (ss)
            ss >>= 1;
        else
            t--;
    }
    for (j = 0; j < LL; j++)
        ran_x[j + KK - LL] = x[j];
    for (; j < KK; j++)
        ran_x[j - LL] = x[j];
    for (j = 0; j < 10; j++)
        ran_array(x, KK + KK - 1); /* warm things up */
    ran_arr_ptr = &ran_arr_started;
}

#define ran_arr_next() (*ran_arr_ptr >= 0 ? *ran_arr_ptr++ : ran_arr_cycle())
long ran_arr_cycle(void) {
    if (ran_arr_ptr == &ran_arr_dummy)
        ran_start(314159L); /* the user forgot to initialize */
    ran_array(ran_arr_buf, QUALITY);
    ran_arr_buf[KK] = -1;
    ran_arr_ptr = ran_arr_buf + 1;
    return ran_arr_buf[0];
}


/* uniform random number in [-1,1] */
static scs_float rand_scs_float(void) {
    return 2 * (((scs_float)ran_arr_next()) / _MAX_RAND_VAL) - 1; /* in [-1, 1] */
}

void SCS(finish_cone)(ScsConeWork* c) {
#ifdef USE_LAPACK
    if (c->Xs) {
        scs_free(c->Xs);
    }
    if (c->Z) {
        scs_free(c->Z);
    }
    if (c->e) {
        scs_free(c->e);
    }
    if (c->work) {
        scs_free(c->work);
    }
#endif
    if (c->cone_boundaries) {
        scs_free(c->cone_boundaries);
    }
    if (c->s) {
        scs_free(c->s);
    }
    if (c) {
        scs_free(c);
    }
}

static inline scs_int get_sd_cone_size(scs_int s) {
    return (s * (s + 1)) / 2;
}

/*
 * boundaries will contain array of indices of rows of A corresponding to
 * cone boundaries, boundaries[0] is starting index for cones of size strictly
 * larger than 1, boundaries malloc-ed here so should be freed.
 */
void set_cone_boundaries(const ScsCone* k, ScsConeWork* c) {
    scs_int i, s_cone_sz, count = 0;
    scs_int cone_boundaries_len =
        1 + k->qsize + k->ssize + k->ed + k->ep + k->psize;
    scs_int* b = (scs_int*)scs_calloc(cone_boundaries_len, sizeof(scs_int));
    /* cones that can be scaled independently */
    b[count] = k->z + k->l + k->bsize;
    count += 1; /* started at 0 now move to first entry */
    for (i = 0; i < k->qsize; ++i) {
        b[count + i] = k->q[i];
    }
    count += k->qsize;
    for (i = 0; i < k->ssize; ++i) {
        s_cone_sz = get_sd_cone_size(k->s[i]);
        b[count + i] = s_cone_sz;
    }
    count += k->ssize; /* add ssize here not ssize * (ssize + 1) / 2 */
    /* exp cones */
    for (i = 0; i < k->ep + k->ed; ++i) {
        b[count + i] = 3;
    }
    count += k->ep + k->ed;
    /* power cones */
    for (i = 0; i < k->psize; ++i) {
        b[count + i] = 3;
    }
    count += k->psize;
    /* other cones */
    c->cone_boundaries = b;
    c->cone_boundaries_len = cone_boundaries_len;
}

static scs_int set_up_sd_cone_work_space(ScsConeWork* c, const ScsCone* k) {
    scs_int i;
#ifdef USE_LAPACK
    blas_int n_max = 0;
    blas_int neg_one = -1;
    blas_int info = 0;
    scs_float wkopt = 0.0;
#if VERBOSITY > 0
#define _STR_EXPAND(tok) #tok
#define _STR(tok) _STR_EXPAND(tok)
    scs_printf("BLAS(func) = '%s'\n", _STR(BLAS(func)));
#endif
    /* eigenvector decomp workspace */
    for (i = 0; i < k->ssize; ++i) {
        if (k->s[i] > n_max) {
            n_max = (blas_int)k->s[i];
        }
    }
    c->Xs = (scs_float*)scs_calloc(n_max * n_max, sizeof(scs_float));
    c->Z = (scs_float*)scs_calloc(n_max * n_max, sizeof(scs_float));
    c->e = (scs_float*)scs_calloc(n_max, sizeof(scs_float));

    /* workspace query */
    BLAS(syev)
        ("Vectors", "Lower", &n_max, c->Xs, &n_max, SCS_NULL, &wkopt, &neg_one,
            &info);

    if (info != 0) {
        scs_printf("FATAL: syev failure, info = %li\n", (long)info);
        return -1;
    }
    c->lwork = (blas_int)(wkopt + 1); /* +1 for int casting safety */
    c->work = (scs_float*)scs_calloc(c->lwork, sizeof(scs_float));

    if (!c->Xs || !c->Z || !c->e || !c->work) {
        return -1;
    }
    return 0;
#else
    for (i = 0; i < k->ssize; i++) {
        if (k->s[i] > 1) {
            scs_printf(
                "FATAL: Cannot solve SDPs without linked blas+lapack libraries\n");
            scs_printf(
                "Install blas+lapack and re-compile SCS with blas+lapack library "
                "locations\n");
            return -1;
        }
    }
    return 0;
#endif
}
ScsConeWork* SCS(init_cone)(ScsCone* k, scs_int m) {
    ScsConeWork* c = (ScsConeWork*)scs_calloc(1, sizeof(ScsConeWork));
    c->k = k;
    c->m = m;
    c->scaled_cones = 0;
    set_cone_boundaries(k, c);
    c->s = (scs_float*)scs_calloc(m, sizeof(scs_float));
    if (k->ssize && k->s) {
        if (set_up_sd_cone_work_space(c, k) < 0) {
            SCS(finish_cone)(c);
            return SCS_NULL;
        }
    }
    return c;
}

/*
 * Routine to scale the limits of the box cone by the scaling diagonal mat D > 0
 *
 *  want (t, s) \in K <==> (t', s') \in K'
 *
 *  (t', s') = (d0 * t, D s) (overloading D to mean D[1:])
 *    (up to scalar scaling factor which we can ignore due to conic prooperty)
 *
 *   K = { (t, s) | t * l <= s <= t * u, t >= 0 } =>
 *       { (t, s) | d0 * t * D l / d0 <= D s <= d0 * t D u / d0, t >= 0 } =>
 *       { (t', s') | t' * l' <= s' <= t' u', t >= 0 } = K'
 *  where l' = D l  / d0, u' = D u / d0.
 */
static void normalize_box_cone(ScsCone* k, scs_float* D, scs_int bsize) {
    scs_int j;
    for (j = 0; j < bsize - 1; j++) {
        if (k->bu[j] >= MAX_BOX_VAL) {
            k->bu[j] = INFINITY;
        }
        else {
            k->bu[j] = D ? D[j + 1] * k->bu[j] / D[0] : k->bu[j];
        }
        if (k->bl[j] <= -MAX_BOX_VAL) {
            k->bl[j] = -INFINITY;
        }
        else {
            k->bl[j] = D ? D[j + 1] * k->bl[j] / D[0] : k->bl[j];
        }
    }
}

void scale_box_cone(ScsCone* k, ScsConeWork* c, ScsScaling* scal) {
    if (k->bsize && k->bu && k->bl) {
        c->box_t_warm_start = 1.;
        if (scal) {
            /* also does some sanitizing */
            normalize_box_cone(k, &(scal->D[k->z + k->l]), k->bsize);
        }
    }
}

/* Project onto { (t, s) | t * l <= s <= t * u, t >= 0 }, Newton's method on t
   tx = [t; s], total length = bsize, under Euclidean metric 1/r_box.
*/
static scs_float proj_box_cone(scs_float* tx, const scs_float* bl,
    const scs_float* bu, scs_int bsize,
    scs_float t_warm_start, scs_float* r_box) {
    scs_float* x, gt, ht, t_prev, t = t_warm_start;
    scs_float rho_t = 1, * rho = SCS_NULL, r;
    scs_int iter, j;

    if (bsize == 1) { /* special case */
        tx[0] = MAX(tx[0], 0.0);
        return tx[0];
    }
    x = &(tx[1]);

    if (r_box) {
        rho_t = 1.0 / r_box[0];
        rho = &(r_box[1]);
    }

    /* should only require about 5 or so iterations, 1 or 2 if warm-started */
    for (iter = 0; iter < BOX_CONE_MAX_ITERS; iter++) {
        t_prev = t;
        gt = rho_t * (t - tx[0]); /* gradient */
        ht = rho_t;               /* hessian */
        for (j = 0; j < bsize - 1; j++) {
            r = rho ? 1.0 / rho[j] : 1.;
            if (x[j] > t * bu[j]) {
                gt += r * (t * bu[j] - x[j]) * bu[j]; /* gradient */
                ht += r * bu[j] * bu[j];              /* hessian */
            }
            else if (x[j] < t * bl[j]) {
                gt += r * (t * bl[j] - x[j]) * bl[j]; /* gradient */
                ht += r * bl[j] * bl[j];              /* hessian */
            }
        }
        t = MAX(t - gt / MAX(ht, 1e-8), 0.); /* newton step */
#if VERBOSITY > 3
        scs_printf("iter %i, t_new %1.3e, t_prev %1.3e, gt %1.3e, ht %1.3e\n", iter,
            t, t_prev, gt, ht);
        scs_printf("ABS(gt / (ht + 1e-6)) %.4e, ABS(t - t_prev) %.4e\n",
            ABS(gt / (ht + 1e-6)), ABS(t - t_prev));
#endif
        /* TODO: sometimes this check can fail (ie, declare convergence before it
         * should) if ht is very large, which can happen with some pathological
         * problems.
         */
        if (ABS(gt / MAX(ht, 1e-6)) < 1e-12 * MAX(t, 1.) ||
            ABS(t - t_prev) < 1e-11 * MAX(t, 1.)) {
            break;
        }
    }
    if (iter == BOX_CONE_MAX_ITERS) {
        scs_printf("warning: box cone proj hit maximum %i iters\n", (int)iter);
    }
    for (j = 0; j < bsize - 1; j++) {
        if (x[j] > t * bu[j]) {
            x[j] = t * bu[j];
        }
        else if (x[j] < t * bl[j]) {
            x[j] = t * bl[j];
        }
        /* x[j] unchanged otherwise */
    }
    tx[0] = t;

#if VERBOSITY > 3
    scs_printf("box cone iters %i\n", (int)iter + 1);
#endif
    return t;
}
/* Self-rolled basic linear algebra routines */

/* a *= b */
void SCS(scale_array)(scs_float* a, const scs_float b, scs_int len) {
    scs_int i;
    for (i = 0; i < len; ++i)
        a[i] *= b;
}

/* x'*y */
scs_float SCS(dot)(const scs_float* x, const scs_float* y, scs_int len) {
    scs_int i;
    scs_float ip = 0.0;
    for (i = 0; i < len; ++i) {
        ip += x[i] * y[i];
    }
    return ip;
}

/* ||v||_2^2 */
scs_float SCS(norm_sq)(const scs_float* v, scs_int len) {
    scs_int i;
    scs_float nmsq = 0.0;
    for (i = 0; i < len; ++i) {
        nmsq += v[i] * v[i];
    }
    return nmsq;
}

/* ||v||_2 */
scs_float SCS(norm_2)(const scs_float* v, scs_int len) {
    return SQRTF(SCS(norm_sq)(v, len));
}

scs_float SCS(norm_inf)(const scs_float* a, scs_int len) {
    scs_float tmp, max = 0.0;
    scs_int i;
    for (i = 0; i < len; ++i) {
        tmp = ABS(a[i]);
        if (tmp > max) {
            max = tmp;
        }
    }
    return max;
}

/* axpy a += sc*b */
void SCS(add_scaled_array)(scs_float* a, const scs_float* b, scs_int n,
    const scs_float sc) {
    scs_int i;
    for (i = 0; i < n; ++i) {
        a[i] += sc * b[i];
    }
}

scs_float SCS(mean)(const scs_float* x, scs_int n) {
    scs_int i;
    scs_float mean = 0.;
    for (i = 0; i < n; ++i) {
        mean += x[i];
    }
    return mean / n;
}

/* project onto SOC of size q*/
static void proj_soc(scs_float* x, scs_int q) {
    if (q == 0) {
        return;
    }
    if (q == 1) {
        x[0] = MAX(x[0], 0.);
        return;
    }
    scs_float v1 = x[0];
    scs_float s = SCS(norm_2)(&(x[1]), q - 1);
    scs_float alpha = (s + v1) / 2.0;

    if (s <= v1) {
        return;
    }
    else if (s <= -v1) {
        memset(&(x[0]), 0, q * sizeof(scs_float));
    }
    else {
        x[0] = alpha;
        SCS(scale_array)(&(x[1]), alpha / s, q - 1);
    }
}

/* size of X is get_sd_cone_size(n) */
static scs_int proj_semi_definite_cone(scs_float* X, const scs_int n,
    ScsConeWork* c) {
    /* project onto the positive semi-definite cone */
#ifdef USE_LAPACK
    scs_int i, first_idx;
    blas_int nb = (blas_int)n;
    blas_int ncols_z;
    blas_int nb_plus_one = (blas_int)(n + 1);
    blas_int one_int = 1;
    scs_float zero = 0., one = 1.;
    scs_float sqrt2 = SQRTF(2.0);
    scs_float sqrt2_inv = 1.0 / sqrt2;
    scs_float* Xs = c->Xs;
    scs_float* Z = c->Z;
    scs_float* e = c->e;
    scs_float* work = c->work;
    blas_int lwork = c->lwork;
    blas_int info = 0;
    scs_float sq_eig_pos;

#endif

    if (n == 0) {
        return 0;
    }
    if (n == 1) {
        X[0] = MAX(X[0], 0.);
        return 0;
    }

#ifdef USE_LAPACK

    /* copy lower triangular matrix into full matrix */
    for (i = 0; i < n; ++i) {
        memcpy(&(Xs[i * (n + 1)]), &(X[i * n - ((i - 1) * i) / 2]),
            (n - i) * sizeof(scs_float));
    }
    /*
       rescale so projection works, and matrix norm preserved
       see http://www.seas.ucla.edu/~vandenbe/publications/mlbook.pdf pg 3
     */
     /* scale diags by sqrt(2) */
    BLAS(scal)(&nb, &sqrt2, Xs, &nb_plus_one); /* not n_squared */

    /* Solve eigenproblem, reuse workspaces */
    BLAS(syev)("Vectors", "Lower", &nb, Xs, &nb, e, work, &lwork, &info);
    if (info != 0) {
        scs_printf("WARN: LAPACK syev error, info = %i\n", (int)info);
        if (info < 0) {
            return info;
        }
    }

    first_idx = -1;
    /* e is eigvals in ascending order, find first entry > 0 */
    for (i = 0; i < n; ++i) {
        if (e[i] > 0) {
            first_idx = i;
            break;
        }
    }

    if (first_idx == -1) {
        /* there are no positive eigenvalues, set X to 0 and return */
        memset(X, 0, sizeof(scs_float) * get_sd_cone_size(n));
        return 0;
    }

    /* Z is matrix of eigenvectors with positive eigenvalues */
    memcpy(Z, &Xs[first_idx * n], sizeof(scs_float) * n * (n - first_idx));

    /* scale Z by sqrt(eig) */
    for (i = first_idx; i < n; ++i) {
        sq_eig_pos = SQRTF(e[i]);
        BLAS(scal)(&nb, &sq_eig_pos, &Z[(i - first_idx) * n], &one_int);
    }

    /* Xs = Z Z' = V E V' */
    ncols_z = (blas_int)(n - first_idx);
    BLAS(syrk)("Lower", "NoTrans", &nb, &ncols_z, &one, Z, &nb, &zero, Xs, &nb);

    /* undo rescaling: scale diags by 1/sqrt(2) */
    BLAS(scal)(&nb, &sqrt2_inv, Xs, &nb_plus_one); /* not n_squared */

    /* extract just lower triangular matrix */
    for (i = 0; i < n; ++i) {
        memcpy(&(X[i * n - ((i - 1) * i) / 2]), &(Xs[i * (n + 1)]),
            (n - i) * sizeof(scs_float));
    }
    return 0;

#else
    scs_printf("FAILURE: solving SDP but no blas/lapack libraries were found!\n");
    scs_printf("SCS will return nonsense!\n");
    SCS(scale_array)(X, NAN, n);
    return -1;
#endif
}

static scs_float exp_newton_one_d(scs_float rho, scs_float y_hat,
    scs_float z_hat, scs_float w) {
    scs_float t_prev, t = MAX(w - z_hat, MAX(-z_hat, 1e-9));
    scs_float f = 1., fp = 1.;
    scs_int i;
    for (i = 0; i < EXP_CONE_MAX_ITERS; ++i) {
        t_prev = t;
        f = t * (t + z_hat) / rho / rho - y_hat / rho + log(t / rho) + 1;
        fp = (2 * t + z_hat) / rho / rho + 1 / t;

        t = t - f / fp;

        if (t <= -z_hat) {
            t = -z_hat;
            break;
        }
        else if (t <= 0) {
            t = 0;
            break;
        }
        else if (ABS(t - t_prev) < CONE_TOL) {
            break;
        }
        else if (SQRTF(f * f / fp) < CONE_TOL) {
            break;
        }
    }
    if (i == EXP_CONE_MAX_ITERS) {
        scs_printf("warning: exp cone newton step hit maximum %i iters\n", (int)i);
        scs_printf("rho=%1.5e; y_hat=%1.5e; z_hat=%1.5e; w=%1.5e; f=%1.5e, "
            "fp=%1.5e, t=%1.5e, t_prev= %1.5e\n",
            rho, y_hat, z_hat, w, f, fp, t, t_prev);
    }
    return t + z_hat;
}

static void exp_solve_for_x_with_rho(const scs_float* v, scs_float* x,
    scs_float rho, scs_float w) {
    x[2] = exp_newton_one_d(rho, v[1], v[2], w);
    x[1] = (x[2] - v[2]) * x[2] / rho;
    x[0] = v[0] - rho;
}

static scs_float exp_calc_grad(const scs_float* v, scs_float* x, scs_float rho,
    scs_float w) {
    exp_solve_for_x_with_rho(v, x, rho, w);
    if (x[1] <= 1e-12) {
        return x[0];
    }
    return x[0] + x[1] * log(x[1] / x[2]);
}

static void exp_get_rho_ub(const scs_float* v, scs_float* x, scs_float* ub,
    scs_float* lb) {
    *lb = 0;
    *ub = 0.125;
    while (exp_calc_grad(v, x, *ub, v[1]) > 0) {
        *lb = *ub;
        (*ub) *= 2;
    }
}

/* project onto the exponential cone, v has dimension *exactly* 3 */
static scs_int proj_exp_cone(scs_float* v) {
    scs_int i;
    scs_float ub, lb, rho, g, x[3];
    scs_float r = v[0], s = v[1], t = v[2];

    /* v in cl(Kexp) */
    if ((s * exp(r / s) - t <= CONE_THRESH && s > 0) ||
        (r <= 0 && s == 0 && t >= 0)) {
        return 0;
    }

    /* -v in Kexp^* */
    if ((r > 0 && r * exp(s / r) + exp(1) * t <= CONE_THRESH) ||
        (r == 0 && s <= 0 && t <= 0)) {
        memset(v, 0, 3 * sizeof(scs_float));
        return 0;
    }

    /* special case with analytical solution */
    if (r < 0 && s < 0) {
        v[1] = 0.0;
        v[2] = MAX(v[2], 0);
        return 0;
    }

    /* iterative procedure to find projection, bisects on dual variable: */
    exp_get_rho_ub(v, x, &ub, &lb); /* get starting upper and lower bounds */
    for (i = 0; i < EXP_CONE_MAX_ITERS; ++i) {
        rho = (ub + lb) / 2; /* halfway between upper and lower bounds */
        g = exp_calc_grad(v, x, rho, x[1]); /* calculates gradient wrt dual var */
        if (g > 0) {
            lb = rho;
        }
        else {
            ub = rho;
        }
        if (ub - lb < CONE_TOL) {
            break;
        }
    }
#if VERBOSITY > 10
    scs_printf("exponential cone proj iters %i\n", (int)i);
#endif
    if (i == EXP_CONE_MAX_ITERS) {
        scs_printf("warning: exp cone outer step hit maximum %i iters\n", (int)i);
        scs_printf("r=%1.5e; s=%1.5e; t=%1.5e\n", r, s, t);
    }
    v[0] = x[0];
    v[1] = x[1];
    v[2] = x[2];
    return 0;
}

static scs_float pow_calc_x(scs_float r, scs_float xh, scs_float rh,
    scs_float a) {
    scs_float x = 0.5 * (xh + SQRTF(xh * xh + 4 * a * (rh - r) * r));
    return MAX(x, 1e-12);
}

static scs_float pow_calcdxdr(scs_float x, scs_float xh, scs_float rh,
    scs_float r, scs_float a) {
    return a * (rh - 2 * r) / (2 * x - xh);
}

static scs_float pow_calc_f(scs_float x, scs_float y, scs_float r,
    scs_float a) {
    return POWF(x, a) * POWF(y, (1 - a)) - r;
}

static scs_float pow_calc_fp(scs_float x, scs_float y, scs_float dxdr,
    scs_float dydr, scs_float a) {
    return POWF(x, a) * POWF(y, (1 - a)) * (a * dxdr / x + (1 - a) * dydr / y) -
        1;
}

static void proj_power_cone(scs_float* v, scs_float a) {
    scs_float xh = v[0], yh = v[1], rh = ABS(v[2]);
    scs_float x = 0.0, y = 0.0, r;
    scs_int i;
    /* v in K_a */
    if (xh >= 0 && yh >= 0 &&
        CONE_THRESH + POWF(xh, a) * POWF(yh, (1 - a)) >= rh) {
        return;
    }

    /* -v in K_a^* */
    if (xh <= 0 && yh <= 0 &&
        CONE_THRESH + POWF(-xh, a) * POWF(-yh, 1 - a) >=
        rh * POWF(a, a) * POWF(1 - a, 1 - a)) {
        v[0] = v[1] = v[2] = 0;
        return;
    }

    r = rh / 2;
    for (i = 0; i < POW_CONE_MAX_ITERS; ++i) {
        scs_float f, fp, dxdr, dydr;
        x = pow_calc_x(r, xh, rh, a);
        y = pow_calc_x(r, yh, rh, 1 - a);

        f = pow_calc_f(x, y, r, a);
        if (ABS(f) < CONE_TOL) {
            break;
        }

        dxdr = pow_calcdxdr(x, xh, rh, r, a);
        dydr = pow_calcdxdr(y, yh, rh, r, (1 - a));
        fp = pow_calc_fp(x, y, dxdr, dydr, a);

        r = MAX(r - f / fp, 0);
        r = MIN(r, rh);
    }
    v[0] = x;
    v[1] = y;
    v[2] = (v[2] < 0) ? -(r) : (r);
}

/* project onto the primal K cone in the paper */
/* the r_y vector determines the INVERSE metric, ie, project under the
 * diag(r_y)^-1 norm.
 */
static scs_int proj_cone(scs_float* x, const ScsCone* k, ScsConeWork* c,
    scs_int normalize, scs_float* r_y) {
    scs_int i, status;
    scs_int count = 0;
    scs_float* r_box = SCS_NULL;

    if (k->z) { /* doesn't use r_y */
      /* project onto primal zero / dual free cone */
        memset(x, 0, k->z * sizeof(scs_float));
        count += k->z;
    }

    if (k->l) { /* doesn't use r_y */
      /* project onto positive orthant */
        for (i = count; i < count + k->l; ++i) {
            x[i] = MAX(x[i], 0.0);
        }
        count += k->l;
    }

    if (k->bsize) { /* DOES use r_y */
        if (r_y) {
            r_box = &(r_y[count]);
        }
        /* project onto box cone */
        c->box_t_warm_start = proj_box_cone(&(x[count]), k->bl, k->bu, k->bsize,
            c->box_t_warm_start, r_box);
        count += k->bsize; /* since b = (t,s), len(s) = bsize - 1 */
    }

    if (k->qsize && k->q) { /* doesn't use r_y */
      /* project onto second-order cones */
        for (i = 0; i < k->qsize; ++i) {
            proj_soc(&(x[count]), k->q[i]);
            count += k->q[i];
        }
    }

    if (k->ssize && k->s) { /* doesn't use r_y */
      /* project onto PSD cones */
        for (i = 0; i < k->ssize; ++i) {
            status = proj_semi_definite_cone(&(x[count]), k->s[i], c);
            if (status < 0) {
                return status;
            }
            count += get_sd_cone_size(k->s[i]);
        }
    }

    if (k->ep) { /* doesn't use r_y */
                 /*
                  * exponential cone is not self dual, if s \in K
                  * then y \in K^* and so if K is the primal cone
                  * here we project onto K^*, via Moreau
                  * \Pi_C^*(y) = y + \Pi_C(-y)
                  */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < k->ep; ++i) {
            proj_exp_cone(&(x[count + 3 * i]));
        }
        count += 3 * k->ep;
    }

    /* dual exponential cone */
    if (k->ed) { /* doesn't use r_y */
      /*
       * exponential cone is not self dual, if s \in K
       * then y \in K^* and so if K is the primal cone
       * here we project onto K^*, via Moreau
       * \Pi_C^*(y) = y + \Pi_C(-y)
       */
        scs_int idx;
        scs_float r, s, t;
        SCS(scale_array)(&(x[count]), -1, 3 * k->ed); /* x = -x; */
#ifdef _OPENMP
#pragma omp parallel for private(r, s, t, idx)
#endif
        for (i = 0; i < k->ed; ++i) {
            idx = count + 3 * i;
            r = x[idx];
            s = x[idx + 1];
            t = x[idx + 2];

            proj_exp_cone(&(x[idx]));

            x[idx] -= r;
            x[idx + 1] -= s;
            x[idx + 2] -= t;
        }
        count += 3 * k->ed;
    }

    if (k->psize && k->p) { /* doesn't use r_y */
        scs_float v[3];
        scs_int idx;
        /* don't use openmp for power cone
        ifdef _OPENMP
        pragma omp parallel for private(v, idx)
        endif
        */
        for (i = 0; i < k->psize; ++i) { /* doesn't use r_y */
            idx = count + 3 * i;
            if (k->p[i] >= 0) {
                /* primal power cone */
                proj_power_cone(&(x[idx]), k->p[i]);
            }
            else {
                /* dual power cone, using Moreau */
                v[0] = -x[idx];
                v[1] = -x[idx + 1];
                v[2] = -x[idx + 2];

                proj_power_cone(v, -k->p[i]);

                x[idx] += v[0];
                x[idx + 1] += v[1];
                x[idx + 2] += v[2];
            }
        }
        count += 3 * k->psize;
    }
    /* project onto OTHER cones */
    return 0;
}

scs_int SCS(proj_dual_cone)(scs_float* x, ScsConeWork* c, ScsScaling* scal,
    scs_float* r_y) {
    scs_int status, i;
    ScsCone* k = c->k;

    if (!c->scaled_cones) {
        scale_box_cone(k, c, scal);
        c->scaled_cones = 1;
    }

    /* copy s = x */
    memcpy(c->s, x, c->m * sizeof(scs_float));

    /* x -> - Rx */
    for (i = 0; i < c->m; ++i) {
        x[i] *= r_y ? -r_y[i] : -1;
    }

    /* project -x onto cone, x -> \Pi_{C^*}^{R^{-1}}(-x) under r_y metric */
    status = proj_cone(x, k, c, scal ? 1 : 0, r_y);

    /* return x + R^{-1} \Pi_{C^*}^{R^{-1}} ( -x )  */
    for (i = 0; i < c->m; ++i) {
        if (r_y) {
            x[i] = x[i] / r_y[i] + c->s[i];
        }
        else {
            x[i] += c->s[i];
        }
    }

    return status;
}

void gen_random_prob_data(scs_int nnz, scs_int col_nnz, ScsData* d, ScsCone* k,
    ScsSolution* opt_sol, scs_int seed) {
    scs_int n = d->n;
    scs_int m = d->m;
    ScsMatrix* A = d->A = (ScsMatrix*)scs_calloc(1, sizeof(ScsMatrix));
    scs_float* b = d->b = (scs_float*)scs_calloc(m, sizeof(scs_float));
    scs_float* c = d->c = (scs_float*)scs_calloc(n, sizeof(scs_float));
    scs_float* x = opt_sol->x = (scs_float*)scs_calloc(n, sizeof(scs_float));
    scs_float* y = opt_sol->y = (scs_float*)scs_calloc(m, sizeof(scs_float));
    scs_float* s = opt_sol->s = (scs_float*)scs_calloc(m, sizeof(scs_float));
    /* temporary variables */
    scs_float* z = (scs_float*)scs_calloc(m, sizeof(scs_float));
    ScsConeWork* tmp_cone_work;
    scs_int i, j, r, rn, rm;

    A->i = (scs_int*)scs_calloc(nnz, sizeof(scs_int));
    A->p = (scs_int*)scs_calloc((n + 1), sizeof(scs_int));
    A->x = (scs_float*)scs_calloc(nnz, sizeof(scs_float));
    A->n = d->n;
    A->m = d->m;
    /* y, s >= 0 and y'*s = 0 */
    for (i = 0; i < m; i++) {
        y[i] = z[i] = rand_scs_float();
    }
    
    tmp_cone_work = SCS(init_cone)(k, m);
    SCS(proj_dual_cone)(y, tmp_cone_work, SCS_NULL, SCS_NULL);
    SCS(finish_cone)(tmp_cone_work);

    for (i = 0; i < m; i++) {
        b[i] = s[i] = y[i] - z[i];
    }

    for (i = 0; i < n; i++) {
        x[i] = rand_scs_float();
    }

    /*
     c = -A'*y
     b = A*x + s
     */
    ran_start(seed);
    A->p[0] = 0;
    for (j = 0; j < n; j++) { /* column */
        r = 0;
        for (i = 0; i < m && r < col_nnz; ++i) {
            /* generate a unique sorted array via Knuths alg */
            rn = m - i;
            rm = col_nnz - r;
            if ((ran_arr_next() % rn) < rm) {
                A->x[r + j * col_nnz] = rand_scs_float();
                A->i[r + j * col_nnz] = i;
                b[i] += A->x[r + j * col_nnz] * x[j];
                c[j] -= A->x[r + j * col_nnz] * y[i];
                r++;
            }
        }
        A->p[j + 1] = (j + 1) * col_nnz;
    }
    scs_free(z);
}

static void print_d_A(const ScsMatrix* A, int nnonz, int n) {
    scs_printf("matrix A:\n");
    scs_printf("x: ");
    for (int i = 0; i < nnonz; i++) {
        scs_printf("%6.2f ", A->x[i]);
    }
    scs_printf("\n");
    scs_printf("i: ");
    for (int i = 0; i < nnonz; i++) {
        scs_printf("%6d ", A->i[i]);
    }
    scs_printf("\n");
    scs_printf("p: ");
    for (int i = 0; i <= n; i++) {
        scs_printf("%4d ", A->p[i]);
    }
    scs_printf("\n");
}

static void print_d_b(const scs_float* b, int m) {
    scs_printf("vector b:\n");
    for (int i = 0; i < m; i++) {
        scs_printf("%8.2f ", b[i]);
    }
    scs_printf("\n");
}

static void print_d_c(const scs_float* c, int n) {
    scs_printf("vector c:\n");
    for (int i = 0; i < n; i++) {
        scs_printf("%8.2f ", c[i]);
    }
    scs_printf("\n");
}

static void print_d(const ScsData* d, int nnonz) {
    scs_printf("m: %d, n: %d\n", d->m, d->n);
    print_d_A(d->A, nnonz, d->n);
    print_d_b(d->b, d->m);
    print_d_c(d->c, d->n);
}

static void print_sol_prim(scs_float* prim_sol, int n) {
    scs_printf("Primal Solution(s):\n");
    for (int i = 0; i < n; i++) {
        scs_printf("x[%d]: %8.4f\n", i, prim_sol[i]);
    }
}

static void print_sol_dual(scs_float* dual_sol, int m) {
    scs_printf("Dual Solution(s):\n");
    for (int i = 0; i < m; i++) {
        scs_printf("y[%d]: %8.4f\n", i, dual_sol[i]);
    }
}

SCIP_RETCODE execmain(int argc, const char** argv) {
    ScsCone* k = (ScsCone*)scs_calloc(1, sizeof(ScsCone));
    ScsData* d = (ScsData*)scs_calloc(1, sizeof(ScsData));
    ScsSettings* stgs = (ScsSettings*)scs_calloc(1, sizeof(ScsSettings));
    ScsSolution* sol = (ScsSolution*)scs_calloc(1, sizeof(ScsSolution));
    ScsSolution* opt_sol = (ScsSolution*)scs_calloc(1, sizeof(ScsSolution));
    ScsInfo info = { 0 };
    const scs_float p_f = 0.1;
    int seed = 1234;
    const scs_int n = 10;
    const scs_int m = 20;
    const scs_int col_nnz = (scs_int)ceil(sqrt(n));
    const scs_int nnz = n * col_nnz;
    scs_int exitflag;
    scs_float perr, derr;
    scs_int success;
    const char* fail;

    k->z = m; // (scs_int)floor(m * p_f);
    k->l = m - k->z;

    d->m = m;
    d->n = n;
    gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, seed);
    print_d(d, nnz);

    SCIP_LPI* lpi;
    SCIP_CALL(SCIPlpiCreate(&lpi, NULL, "prob", SCIP_OBJSEN_MINIMIZE));
    double lb = -SCIPlpiInfinity(lpi);
    double ub = SCIPlpiInfinity(lpi);
    for (int i = 0; i < n; i++) {
        SCIP_CALL(SCIPlpiAddCols(lpi, 1, &d->c[i], &lb, &ub, NULL, 0, NULL, NULL, NULL));
    }
    
    std::map<int, std::vector<std::pair<int, double>>> cons;
    for (int i = 0; i < d->A->p[d->n]; i += col_nnz) {
        int v = i / col_nnz;
        for (int j = i; j < i + col_nnz; j++) {
            auto it = cons.find(d->A->i[j]);
            if (it == cons.end()) {
                cons.insert({ d->A->i[j], {{ v, d->A->x[j] }} });
                continue;
            }
            it->second.push_back({v, d->A->x[j]});
        }
    }
    for (auto i = cons.begin(); i != cons.end(); i++) {
        int beg = 0;
        int ind[i->second.size()];
        double val[i->second.size()];
        int p = 0;
        std::cout << i->first << ":";
        for (auto j : i->second) {
            std::cout << "(" << j.first << "," << j.second << ") ";
            ind[p] = j.first;
            val[p] = j.second;
            p++;
        }
        std::cout << std::endl;
        SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lb, &d->b[i->first], NULL, i->second.size(), &beg, ind, val));
    }
    /**
    for (int i = 0; i < d->A->p[d->n]; i++) {
        auto it = cons.find(d->A->i[i]);
        if (it != cons.end()) {
            //std::cout << it->first << " size:" << it->second.size() << std::endl;
            it->second.push_back(d->A->x[i]);
        }
        else {
            cons.insert({ d->A->i[i], {d->A->x[i]} });
        }
    }
    /**
    std::cout << cons.size() << std::endl;
    for (auto it = cons.begin(); it != cons.end(); it++) {
        std::cout << it->first << ":";
        for (auto i : it->second) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }
    /**
    for (int i = 0; i < n; i += col_nnz) {
        int beg = 0;
        int ind[col_nnz];
        double val[col_nnz];
        for (int j = 0; j < col_nnz; j++) {
            ind[j] = d->A->i[i + j];
            val[j] = d->A->x[i + j];
        }
        SCIP_CALL(SCIPlpiAddRows(lpi, 1, &lb, &d->b[i], NULL, col_nnz, &beg, ind, val));
    }
    */
    SCIP_CALL(SCIPlpiSolvePrimal(lpi));
    double objval[1];
    double primsol[n];
    double dualsol[m];
    SCIP_CALL(SCIPlpiGetSol(lpi, objval, primsol, dualsol, NULL, NULL));
    scs_printf("Objective: %8.4f\n", *objval);
    print_sol_prim(primsol, n);
    print_sol_dual(dualsol, m);
    /**
    scs_set_default_settings(stgs);
    stgs->eps_abs = 1e-5;
    stgs->eps_rel = 1e-5;

    exitflag = scs(d, k, stgs, sol, &info);

    perr = SCS(dot)(d->c, sol->x, d->n) - SCS(dot)(d->c, opt_sol->x, d->n);
    derr = -SCS(dot)(d->b, sol->y, d->m) + SCS(dot)(d->b, opt_sol->y, d->m);
    scs_printf("true obj %4e\n", SCS(dot)(d->c, opt_sol->x, d->n));
    scs_printf("primal obj error %4e\n", perr);
    scs_printf("dual obj error %4e\n", derr);

    success = ABS(perr) < 1e-4 && ABS(derr) < 1e-4 && exitflag == SCS_SOLVED;

    //mu_assert("small_lp: SCS failed to produce outputflag SCS_SOLVED", success);
    fail = verify_solution_correct(d, k, stgs, &info, sol, exitflag);
    SCS(free_data)(d);
    SCS(free_cone)(k);
    SCS(free_sol)(sol);
    SCS(free_sol)(opt_sol);
    scs_free(stgs);*/
    return SCIP_OKAY;
}

int main(int argc, const char* argv[]) {
    printf("Hello, SCIP! This problem would be solved by using SCIP integrated with SCS.\n");
    return execmain(argc, argv) != SCIP_OKAY ? 1 : 0;
}