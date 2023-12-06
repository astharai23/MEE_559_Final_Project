/* This file was automatically generated by CasADi 3.6.3.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int int
#endif

int MPC_simplified_form_H(casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int MPC_simplified_form_H_alloc_mem(void);
int MPC_simplified_form_H_init_mem(int mem);
void MPC_simplified_form_H_free_mem(int mem);
int MPC_simplified_form_H_checkout(void);
void MPC_simplified_form_H_release(int mem);
void MPC_simplified_form_H_incref(void);
void MPC_simplified_form_H_decref(void);
casadi_int MPC_simplified_form_H_n_in(void);
casadi_int MPC_simplified_form_H_n_out(void);
casadi_real MPC_simplified_form_H_default_in(casadi_int i);
const char* MPC_simplified_form_H_name_in(casadi_int i);
const char* MPC_simplified_form_H_name_out(casadi_int i);
const casadi_int* MPC_simplified_form_H_sparsity_in(casadi_int i);
const casadi_int* MPC_simplified_form_H_sparsity_out(casadi_int i);
int MPC_simplified_form_H_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define MPC_simplified_form_H_SZ_ARG 15
#define MPC_simplified_form_H_SZ_RES 2
#define MPC_simplified_form_H_SZ_IW 0
#define MPC_simplified_form_H_SZ_W 2
int MPC_simplified_form_g(casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int MPC_simplified_form_g_alloc_mem(void);
int MPC_simplified_form_g_init_mem(int mem);
void MPC_simplified_form_g_free_mem(int mem);
int MPC_simplified_form_g_checkout(void);
void MPC_simplified_form_g_release(int mem);
void MPC_simplified_form_g_incref(void);
void MPC_simplified_form_g_decref(void);
casadi_int MPC_simplified_form_g_n_in(void);
casadi_int MPC_simplified_form_g_n_out(void);
casadi_real MPC_simplified_form_g_default_in(casadi_int i);
const char* MPC_simplified_form_g_name_in(casadi_int i);
const char* MPC_simplified_form_g_name_out(casadi_int i);
const casadi_int* MPC_simplified_form_g_sparsity_in(casadi_int i);
const casadi_int* MPC_simplified_form_g_sparsity_out(casadi_int i);
int MPC_simplified_form_g_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define MPC_simplified_form_g_SZ_ARG 15
#define MPC_simplified_form_g_SZ_RES 2
#define MPC_simplified_form_g_SZ_IW 0
#define MPC_simplified_form_g_SZ_W 5
int MPC_simplified_form_HN(casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int MPC_simplified_form_HN_alloc_mem(void);
int MPC_simplified_form_HN_init_mem(int mem);
void MPC_simplified_form_HN_free_mem(int mem);
int MPC_simplified_form_HN_checkout(void);
void MPC_simplified_form_HN_release(int mem);
void MPC_simplified_form_HN_incref(void);
void MPC_simplified_form_HN_decref(void);
casadi_int MPC_simplified_form_HN_n_in(void);
casadi_int MPC_simplified_form_HN_n_out(void);
casadi_real MPC_simplified_form_HN_default_in(casadi_int i);
const char* MPC_simplified_form_HN_name_in(casadi_int i);
const char* MPC_simplified_form_HN_name_out(casadi_int i);
const casadi_int* MPC_simplified_form_HN_sparsity_in(casadi_int i);
const casadi_int* MPC_simplified_form_HN_sparsity_out(casadi_int i);
int MPC_simplified_form_HN_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define MPC_simplified_form_HN_SZ_ARG 7
#define MPC_simplified_form_HN_SZ_RES 1
#define MPC_simplified_form_HN_SZ_IW 0
#define MPC_simplified_form_HN_SZ_W 1
int MPC_simplified_form_gN(casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int MPC_simplified_form_gN_alloc_mem(void);
int MPC_simplified_form_gN_init_mem(int mem);
void MPC_simplified_form_gN_free_mem(int mem);
int MPC_simplified_form_gN_checkout(void);
void MPC_simplified_form_gN_release(int mem);
void MPC_simplified_form_gN_incref(void);
void MPC_simplified_form_gN_decref(void);
casadi_int MPC_simplified_form_gN_n_in(void);
casadi_int MPC_simplified_form_gN_n_out(void);
casadi_real MPC_simplified_form_gN_default_in(casadi_int i);
const char* MPC_simplified_form_gN_name_in(casadi_int i);
const char* MPC_simplified_form_gN_name_out(casadi_int i);
const casadi_int* MPC_simplified_form_gN_sparsity_in(casadi_int i);
const casadi_int* MPC_simplified_form_gN_sparsity_out(casadi_int i);
int MPC_simplified_form_gN_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define MPC_simplified_form_gN_SZ_ARG 7
#define MPC_simplified_form_gN_SZ_RES 1
#define MPC_simplified_form_gN_SZ_IW 0
#define MPC_simplified_form_gN_SZ_W 5
int MPC_simplified_form_dynamics(casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int MPC_simplified_form_dynamics_alloc_mem(void);
int MPC_simplified_form_dynamics_init_mem(int mem);
void MPC_simplified_form_dynamics_free_mem(int mem);
int MPC_simplified_form_dynamics_checkout(void);
void MPC_simplified_form_dynamics_release(int mem);
void MPC_simplified_form_dynamics_incref(void);
void MPC_simplified_form_dynamics_decref(void);
casadi_int MPC_simplified_form_dynamics_n_in(void);
casadi_int MPC_simplified_form_dynamics_n_out(void);
casadi_real MPC_simplified_form_dynamics_default_in(casadi_int i);
const char* MPC_simplified_form_dynamics_name_in(casadi_int i);
const char* MPC_simplified_form_dynamics_name_out(casadi_int i);
const casadi_int* MPC_simplified_form_dynamics_sparsity_in(casadi_int i);
const casadi_int* MPC_simplified_form_dynamics_sparsity_out(casadi_int i);
int MPC_simplified_form_dynamics_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define MPC_simplified_form_dynamics_SZ_ARG 15
#define MPC_simplified_form_dynamics_SZ_RES 4
#define MPC_simplified_form_dynamics_SZ_IW 0
#define MPC_simplified_form_dynamics_SZ_W 21
int MPC_simplified_form_bounds(casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int MPC_simplified_form_bounds_alloc_mem(void);
int MPC_simplified_form_bounds_init_mem(int mem);
void MPC_simplified_form_bounds_free_mem(int mem);
int MPC_simplified_form_bounds_checkout(void);
void MPC_simplified_form_bounds_release(int mem);
void MPC_simplified_form_bounds_incref(void);
void MPC_simplified_form_bounds_decref(void);
casadi_int MPC_simplified_form_bounds_n_in(void);
casadi_int MPC_simplified_form_bounds_n_out(void);
casadi_real MPC_simplified_form_bounds_default_in(casadi_int i);
const char* MPC_simplified_form_bounds_name_in(casadi_int i);
const char* MPC_simplified_form_bounds_name_out(casadi_int i);
const casadi_int* MPC_simplified_form_bounds_sparsity_in(casadi_int i);
const casadi_int* MPC_simplified_form_bounds_sparsity_out(casadi_int i);
int MPC_simplified_form_bounds_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define MPC_simplified_form_bounds_SZ_ARG 15
#define MPC_simplified_form_bounds_SZ_RES 4
#define MPC_simplified_form_bounds_SZ_IW 0
#define MPC_simplified_form_bounds_SZ_W 4
int MPC_simplified_form_boundsN(casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int MPC_simplified_form_boundsN_alloc_mem(void);
int MPC_simplified_form_boundsN_init_mem(int mem);
void MPC_simplified_form_boundsN_free_mem(int mem);
int MPC_simplified_form_boundsN_checkout(void);
void MPC_simplified_form_boundsN_release(int mem);
void MPC_simplified_form_boundsN_incref(void);
void MPC_simplified_form_boundsN_decref(void);
casadi_int MPC_simplified_form_boundsN_n_in(void);
casadi_int MPC_simplified_form_boundsN_n_out(void);
casadi_real MPC_simplified_form_boundsN_default_in(casadi_int i);
const char* MPC_simplified_form_boundsN_name_in(casadi_int i);
const char* MPC_simplified_form_boundsN_name_out(casadi_int i);
const casadi_int* MPC_simplified_form_boundsN_sparsity_in(casadi_int i);
const casadi_int* MPC_simplified_form_boundsN_sparsity_out(casadi_int i);
int MPC_simplified_form_boundsN_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#define MPC_simplified_form_boundsN_SZ_ARG 7
#define MPC_simplified_form_boundsN_SZ_RES 2
#define MPC_simplified_form_boundsN_SZ_IW 0
#define MPC_simplified_form_boundsN_SZ_W 1
#ifdef __cplusplus
} /* extern "C" */
#endif
