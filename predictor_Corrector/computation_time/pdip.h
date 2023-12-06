#include "blasfeo_path/blasfeo.h"

#ifndef QP_REAL
#define QP_REAL double
#endif

typedef struct pdip_options
{
    int maxiter;
    QP_REAL corr_tol;
    QP_REAL pred_tol;
    QP_REAL eta;

    int iter_count;
} pdip_options;

typedef struct pdip_data
{
    int n;
    int m;
    int nx;
    int ni;

    struct blasfeo_dmat sA_big;
    struct blasfeo_dvec sb_big;
    struct blasfeo_dvec sb_big2;

    struct blasfeo_dmat sA_tmp;
    struct blasfeo_dvec sb_tmp;

    struct blasfeo_dvec sg;

    struct blasfeo_dvec sb;

    struct blasfeo_dvec sx;
    struct blasfeo_dvec slam;
    struct blasfeo_dvec sy;

    struct blasfeo_dvec syp;
    struct blasfeo_dvec slamp;

    struct blasfeo_dvec srplam_update1;
    struct blasfeo_dvec srplam_update2;

    struct blasfeo_dvec s_sol;
    struct blasfeo_dvec s_dely;

    QP_REAL *Afinal;
    QP_REAL *bfinal;

} pdip_data;

void pdip_initialize(pdip_data *pd);

void pdip_update(
    pdip_data *pd,
    QP_REAL *H,
    QP_REAL *g,
    QP_REAL *A,
    QP_REAL *lba,
    QP_REAL *uba,
    QP_REAL *lbx,
    QP_REAL *ubx,
    QP_REAL *x,
    QP_REAL *lam,
    QP_REAL *y);

int pdip_solve(pdip_data *pd, pdip_options *options);

void get_solx(pdip_data *pd, QP_REAL *x);

void get_solgamma(pdip_data *pd, QP_REAL *gamma);

void get_sollam(pdip_data *pd, QP_REAL *lam);

void get_soly(pdip_data *pd, QP_REAL *y);

void print_pdip_data(pdip_data *pd);
