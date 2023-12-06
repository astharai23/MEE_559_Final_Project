#include "MPC_simplified_qpoases_prep.h"
#include "qp_pcond.h"
#include "qp_qpoases.h"
#include "qp_solver.h"
#include "pdip.h"
#define __STATIC__ static


void MPC_simplified_call(double *z_op, double *x0, double *pd, double *Q11, double *Q22, double *R, double *umax, double *z_opt, int *exit_flag, double *u0)
{
	__STATIC__ int _init_stat=0;

	int ii, jj;

	__STATIC__ qp_real_t __w__[21];


	__STATIC__ qp_dim dim = {2, 1, 3, 50, 0, 0, 0};

	__STATIC__ qp_elem qp;
	__STATIC__ qp_real_t __z_tmp__[54];
	__STATIC__ qp_dim dim_pcond = {2, 50, 52, 1, 0, 0, 50};

	__STATIC__ qp_elem qp_pcond;
	__STATIC__ qp_pcond_workmem pcond_workmem;
	__STATIC__ qp_generator qp_gen;
	__STATIC__ qp_input_params qp_params;
	__STATIC__ qp_options options;
	__STATIC__ qp_real_t *arg[7];
	__STATIC__ qp_real_t *tmp_arg[8];
	__STATIC__ int p_shift[7] = {3, 0, 0, 0, 0, 0, 0};
	__STATIC__ qpoases_qp_dim dim2 = {54, 2};
	__STATIC__ qpoases_qp_elem qp2;

	__STATIC__ pdip_data pds;

	__STATIC__ pdip_options options_pd;

	__STATIC__ qp_real_t lam[2 * 56], y[2 * 56];
	for (ii = 0; ii < 2 * 56; ii++)
	{
		lam[ii] = 1.0;
		y[ii] = 1.0;
	}


	if (_init_stat == 0)
	{
		_init_stat = 1;

		allocate_qp_mem(&dim, &qp);
		allocate_qp_mem(&dim_pcond, &qp_pcond);
		allocate_pcond_workmem(&dim, &dim_pcond, &pcond_workmem);

		qp_params.n_params = 7;
		qp_params.arg = arg;
		qp_params.p_shift = p_shift;

		arg[0] = z_op;
		arg[1] = x0;
		arg[2] = pd;
		arg[3] = Q11;
		arg[4] = Q22;
		arg[5] = R;
		arg[6] = umax;

		qp_gen.form_H = MPC_simplified_form_H;
		qp_gen.form_g = MPC_simplified_form_g;
		qp_gen.form_HN = MPC_simplified_form_HN;
		qp_gen.form_gN = MPC_simplified_form_gN;
		qp_gen.form_dynamics = MPC_simplified_form_dynamics;
		qp_gen.form_bounds = MPC_simplified_form_bounds;
		qp_gen.form_boundsN = MPC_simplified_form_boundsN;
		qp_gen.form_constraints = 0;
		qp_gen.form_constraintsN = 0;

		// options.max_iter = 20;
		// options.tol = 1e-06;

		compute_qpelements(&dim, &qp, &qp_gen, &qp_params, tmp_arg, __w__);
		pcondense_qp(&dim, &qp, &dim_pcond, &pcond_workmem, &qp_pcond);
		allocate_qpoases_qp_mem(&dim2, &qp2);
		qp_to_qpoases(&dim_pcond, &qp_pcond, &dim2, &qp2);
		// initialize_qpOASES(&qpProb, &dim2, &qp2, &options, z_op);

		pds.n = dim2.nz;
		pds.m = dim2.na;

		options_pd.maxiter = 20;
		options_pd.corr_tol = 1e-3;
		// options.pred_tol = 1e-3;
		options_pd.eta = 0.99;

		pdip_initialize(&pds);
	}
	arg[0] = z_op;
	arg[1] = x0;
	arg[2] = pd;
	arg[3] = Q11;
	arg[4] = Q22;
	arg[5] = R;
	arg[6] = umax;
	compute_qpelements(&dim, &qp, &qp_gen, &qp_params, tmp_arg, __w__);
	pcondense_qp(&dim, &qp, &dim_pcond, &pcond_workmem, &qp_pcond);
	qp_to_qpoases(&dim_pcond, &qp_pcond, &dim2, &qp2);
	// exit_flag[0] = call_qpOASES(&qpProb, &dim2, &qp2, &options, __z_tmp__);

	pdip_update(&pds, qp2.H, qp2.g, qp2.A, qp2.lba, qp2.uba, qp2.zlow, qp2.zupp, __z_tmp__, lam, y);
	exit_flag[0] = pdip_solve(&pds, &options_pd);
	get_solx(&pds, __z_tmp__);
	// get_sollam(&pds, lam);
	// get_soly(&pds, y);
	for (jj = 0; jj < dim_pcond.N; jj++)
	{
		simulate_linear_dynamics(&dim, dim_pcond.block_size, qp.ABt + jj * dim_pcond.block_size * dim.nx * dim.nz, qp.b + jj * dim_pcond.block_size * dim.nx, __z_tmp__ + jj * dim_pcond.nz, z_opt + jj * dim_pcond.block_size * dim.nz);
	}

	u0[0] = z_opt[2];

}
