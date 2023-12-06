#include <math.h>
#include <string.h>
#include "pdip.h"

void pdip_allocate_mem(pdip_data *pd)
{
    blasfeo_allocate_dmat(
        pd->nx + pd->ni,
        pd->nx + pd->ni,
        &pd->sA_big);

    blasfeo_allocate_dvec(
        pd->nx + pd->ni,
        &pd->sb_big);

    blasfeo_allocate_dvec(
        pd->nx + pd->ni,
        &pd->sb_big2);

    blasfeo_allocate_dmat(
        pd->nx + pd->ni,
        pd->nx + pd->ni,
        &pd->sA_tmp);

    blasfeo_allocate_dvec(
        pd->nx + pd->ni,
        &pd->sb_tmp);

    blasfeo_allocate_dvec(
        pd->nx,
        &pd->sg);

    blasfeo_allocate_dvec(
        pd->ni,
        &pd->sb);

    blasfeo_allocate_dvec(
        pd->nx,
        &pd->sx);
    blasfeo_allocate_dvec(
        pd->ni,
        &pd->slam);
    blasfeo_allocate_dvec(
        pd->ni,
        &pd->sy);

    blasfeo_allocate_dvec(
        pd->ni,
        &pd->syp);

    blasfeo_allocate_dvec(
        pd->ni,
        &pd->slamp);

    blasfeo_allocate_dvec(
        pd->ni,
        &pd->srplam_update1);
    blasfeo_allocate_dvec(
        pd->ni,
        &pd->srplam_update2);

    blasfeo_allocate_dvec(
        pd->nx + pd->ni,
        &pd->s_sol);
    blasfeo_allocate_dvec(
        pd->ni,
        &pd->s_dely);

    pd->Afinal = (QP_REAL *) calloc(pd->ni * pd->nx, sizeof(QP_REAL));
    pd->bfinal = (QP_REAL *) calloc(pd->ni, sizeof(QP_REAL));
   
}

void pdip_free_mem(pdip_data *pd)
{
    blasfeo_free(&pd->sA_big);
    blasfeo_free(&pd->sb_big);
    blasfeo_free(&pd->sb_big2);
    blasfeo_free(&pd->sA_tmp);
    blasfeo_free(&pd->sb_tmp);
    blasfeo_free(&pd->sg);
    blasfeo_free(&pd->sb);
    blasfeo_free(&pd->sx);
    blasfeo_free(&pd->slam);
    blasfeo_free(&pd->sy);
    blasfeo_free(&pd->syp);
    blasfeo_free(&pd->slamp);
    blasfeo_free(&pd->srplam_update1);
    blasfeo_free(&pd->srplam_update2);
    blasfeo_free(&pd->s_sol);
    blasfeo_free(&pd->s_dely);

    free(pd->Afinal);
    free(pd->bfinal);
}

void pdip_initialize(pdip_data *pd)
{
    pd->nx = pd->n;
    pd->ni = 2*pd->m + 2*pd->n;
    pdip_allocate_mem(pd);
    blasfeo_dgecpsc(pd->nx + pd->ni, pd->nx + pd->ni, 0.0, &pd->sA_big, 0, 0, &pd->sA_big, 0, 0);
    blasfeo_dveccpsc(pd->nx + pd->ni, 0.0, &pd->sb_big, 0, &pd->sb_big, 0);
}

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
    QP_REAL *y)
{
    int ii, jj;

    for (ii = 0; ii < pd->m; ii++)
    {
        for (jj = 0; jj < pd->n; jj++)
        {
            pd->Afinal[ii*pd->n+jj] = A[ii*pd->n+jj];
        }
        pd->bfinal[ii] = lba[ii];
    }

    for (ii = 0; ii < pd->m; ii++)
    {
        for (jj = 0; jj < pd->n; jj++)
        {
            pd->Afinal[pd->m*pd->n + ii*pd->n+jj] = -A[ii*pd->n+jj];
        }
        pd->bfinal[pd->m+ii] = -uba[ii];
    }

    for (ii = 0; ii < pd->n; ii++)
    {
        for (jj = 0; jj < pd->n; jj++)
        {
            pd->Afinal[2*pd->m*pd->n+ii*pd->n+jj] = ii==jj? 1.0: 0.0;
        }
        pd->bfinal[2*pd->m+ii] = lbx[ii];
    }

    for (ii = 0; ii < pd->n; ii++)
    {
        for (jj = 0; jj < pd->n; jj++)
        {
            pd->Afinal[(2*pd->m+pd->n)*pd->n + ii*pd->n+jj] = ii==jj? -1.0: 0.0;;
        }
        pd->bfinal[2*pd->m+pd->n+ii] = -ubx[ii];
    }
    

    if (H != 0)
    {
        blasfeo_pack_tran_dmat(pd->nx, pd->nx, H, pd->nx, &pd->sA_big, 0, 0);
    }

    blasfeo_pack_tran_dmat(pd->nx, pd->ni, pd->Afinal, pd->nx, &pd->sA_big, pd->nx, 0);
    blasfeo_pack_dmat(pd->nx, pd->ni, pd->Afinal, pd->nx, &pd->sA_big, 0, pd->nx);
    blasfeo_dgecpsc(pd->nx, pd->ni, -1.0, &pd->sA_big, 0, pd->nx, &pd->sA_big, 0, pd->nx);
    

    if (g != 0)
        blasfeo_pack_dvec(pd->nx, g, 1, &pd->sg, 0);

    blasfeo_pack_dvec(pd->ni, pd->bfinal, 1, &pd->sb, 0);

    if (x != 0)
        blasfeo_pack_dvec(pd->nx, x, 1, &pd->sx, 0);
    if (lam != 0)
        blasfeo_pack_dvec(pd->ni, lam, 1, &pd->slam, 0);
    if (y != 0)
        blasfeo_pack_dvec(pd->ni, y, 1, &pd->sy, 0);
}

void pdip_compute_predictor_rhs(pdip_data *pd)
{
    // Create -rd
    blasfeo_dgemv_n(pd->nx, pd->nx, -1.0, &pd->sA_big, 0, 0, &pd->sx, 0, -1.0, &pd->sg, 0, &pd->sb_big, 0);
    blasfeo_dgemv_t(pd->ni, pd->nx, 1.0, &pd->sA_big, pd->nx, 0, &pd->slam, 0, 1.0, &pd->sb_big, 0, &pd->sb_big, 0);

    // Create -rplamda - y
    blasfeo_dgemv_n(pd->ni, pd->nx, -1.0, &pd->sA_big, pd->nx, 0, &pd->sx, 0, 1.0, &pd->sb, 0, &pd->sb_big, pd->nx);
}

void pdip_compute_corrector_rhs(pdip_data *pd, QP_REAL sigma, QP_REAL mu)
{
    int ii;
    blasfeo_dveccp(pd->nx + pd->ni, &pd->sb_big, 0, &pd->sb_big2, 0);

    blasfeo_dvecmul(pd->ni, &pd->s_sol, pd->nx, &pd->s_dely, 0, &pd->srplam_update2, 0);
    for (ii = 0; ii < pd->ni; ii++)
    {
        blasfeo_dvecin1(sigma * mu / blasfeo_dvecex1(&pd->slam, ii),
                        &pd->srplam_update1, ii);
        blasfeo_dvecin1(1 / blasfeo_dvecex1(&pd->slam, ii) * blasfeo_dvecex1(&pd->srplam_update2, ii),
                        &pd->srplam_update2, ii);
    }

    blasfeo_dvecad(pd->ni, 1.0, &pd->srplam_update1, 0, &pd->sb_big2, pd->nx);
    blasfeo_dvecad(pd->ni, -1.0, &pd->srplam_update2, 0, &pd->sb_big2, pd->nx);
}

void blasfeo_linsolve(int m, struct blasfeo_dmat *sA, struct blasfeo_dvec *sb, struct blasfeo_dvec *sx, struct blasfeo_dmat *sA_tmp)
{
    blasfeo_dtrsv_lnu(m, sA_tmp, 0, 0, sb, 0, sx, 0);
    blasfeo_dtrsv_unn(m, sA_tmp, 0, 0, sx, 0, sx, 0);
}

QP_REAL compute_alpha(int m, struct blasfeo_dvec *sy, int yi, struct blasfeo_dvec *sdely, int delyi, QP_REAL eta)
{
    int ii;
    QP_REAL alpha = 1.0, alpha_tmp;
    for (ii = 0; ii < m; ii++)
    {
        if (sdely->pa[delyi + ii] < 0)
        {
            alpha_tmp = -eta * sy->pa[yi + ii] / sdely->pa[delyi + ii];
            alpha = alpha_tmp < alpha ? alpha_tmp : alpha;
        }
    }
    return alpha;
}

QP_REAL compute_norm_inf(int m, struct blasfeo_dvec *sx, int yi)
{
    QP_REAL norm_inf = 0.0, norm_inf_tmp;
    int ii;
    for (ii = 0; ii < m; ii++)
    {
        norm_inf_tmp = fabs(blasfeo_dvecex1(sx, ii));
        norm_inf = norm_inf_tmp > norm_inf ? norm_inf_tmp : norm_inf;
    }
    return norm_inf;
}

int pdip_solve(pdip_data *pd, pdip_options *options)
{
    int iter_count = 0, ii, break_stat = 0;
    QP_REAL alpha_p, alpha_d, mu, mu_p, sigma;
    if (pd->ni == 0)
    {
        /* code for no constraint. Not implemented now */
    }
    else
    {
        while (1)
        {
            for (ii = 0; ii < pd->ni; ii++)
            {
                blasfeo_dgein1(blasfeo_dvecex1(&pd->sy, ii) / blasfeo_dvecex1(&pd->slam, ii), &pd->sA_big, pd->nx + ii, pd->nx + ii);
            }

            mu = blasfeo_dvecmuldot(pd->ni, &pd->sy, 0, &pd->slam, 0, &pd->slamp, 0) / pd->ni;

            pdip_compute_predictor_rhs(pd);

            // LU decomposition
            blasfeo_dgetrf_np(pd->nx + pd->ni, pd->nx + pd->ni, &pd->sA_big, 0, 0, &pd->sA_tmp, 0, 0);
            blasfeo_linsolve(pd->nx + pd->ni, &pd->sA_big, &pd->sb_big, &pd->s_sol, &pd->sA_tmp);

            blasfeo_dgemv_n(pd->ni, pd->nx, 1.0, &pd->sA_big, pd->nx, 0, &pd->s_sol, 0, -1.0, &pd->sb_big, pd->nx, &pd->s_dely, 0);
            blasfeo_dvecad(pd->ni, -1.0, &pd->sy, 0, &pd->s_dely, 0);

            alpha_p = compute_alpha(pd->ni, &pd->sy, 0, &pd->s_dely, 0, options->eta);
            alpha_d = compute_alpha(pd->ni, &pd->slam, 0, &pd->s_sol, pd->nx, options->eta);

            blasfeo_dveccp(pd->ni, &pd->sy, 0, &pd->syp, 0);
            blasfeo_dvecad(pd->ni, alpha_p, &pd->s_dely, 0, &pd->syp, 0);

            blasfeo_dveccp(pd->ni, &pd->slam, 0, &pd->slamp, 0);
            blasfeo_dvecad(pd->ni, alpha_d, &pd->s_sol, pd->nx, &pd->slamp, 0);

            mu_p = blasfeo_dvecmuldot(pd->ni, &pd->syp, 0, &pd->slamp, 0, &pd->slamp, 0) / pd->ni;

            sigma = pow(mu_p / mu, 3);

            pdip_compute_corrector_rhs(pd, sigma, mu);

            blasfeo_linsolve(pd->nx + pd->ni, &pd->sA_big, &pd->sb_big2, &pd->s_sol, &pd->sA_tmp);

            blasfeo_dgemv_n(pd->ni, pd->nx, 1.0, &pd->sA_big, pd->nx, 0, &pd->s_sol, 0, -1.0, &pd->sb_big, pd->nx, &pd->s_dely, 0);
            blasfeo_dvecad(pd->ni, -1.0, &pd->sy, 0, &pd->s_dely, 0);

            alpha_p = compute_alpha(pd->ni, &pd->sy, 0, &pd->s_dely, 0, options->eta);
            alpha_d = compute_alpha(pd->ni, &pd->slam, 0, &pd->s_sol, pd->nx, options->eta);

            blasfeo_dvecad(pd->nx, alpha_p, &pd->s_sol, 0, &pd->sx, 0);
            blasfeo_dvecad(pd->ni, alpha_p, &pd->s_dely, 0, &pd->sy, 0);
            blasfeo_dvecad(pd->ni, alpha_d, &pd->s_sol, pd->nx, &pd->slam, 0);

            iter_count++;

            if (compute_norm_inf(pd->nx, &pd->s_sol, 0) <= options->corr_tol)
            {
                // printf("iter_count: %d\n", iter_count);
                options->iter_count = iter_count;
                return 0;
            }

            if (iter_count >= options->maxiter)
            {
                options->iter_count = iter_count;
                return 1;
            }

            
        }
    }
}

void get_solx(pdip_data *pd, QP_REAL *x)
{
    blasfeo_unpack_dvec(pd->nx, &pd->sx, 0, x, 1);
}

void get_sollam(pdip_data *pd, QP_REAL *lam)
{
    blasfeo_unpack_dvec(pd->ni, &pd->slam, 0, lam, 1);
}

void get_soly(pdip_data *pd, QP_REAL *y)
{
    blasfeo_unpack_dvec(pd->ni, &pd->sy, 0, y, 1);
}

void print_pdip_data(pdip_data *pd)
{
    printf("Abig:\n");
    blasfeo_print_dmat(
        pd->nx + pd->ni,
        pd->nx + pd->ni,
        &pd->sA_big,
        0,
        0);

    printf("bbig:\n");
    blasfeo_print_dvec(
        pd->nx + pd->ni,
        &pd->sb_big,
        0);
}