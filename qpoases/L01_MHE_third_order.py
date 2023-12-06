import casadi as ca
from dynopt import Constraints, Solver, Parameter, Problem
from dynopt.mhe import MHE


N = 14
Ts = 0.02

x = ca.SX.sym('x', 3)
pe = ca.SX.sym('pe')
M_inv = ca.SX.sym('M_inv')
D = ca.SX.sym('D')
Ki = ca.SX.sym("Ki")
Rp = ca.SX.sym("Rp")

x1 = x[0]
x2 = x[1]
x3 = x[2]

# Ki = Parameter('Ki')
Tg = Parameter('Tg')
# Rp = Parameter('Rp')

fx = ca.vertcat(
    -D*M_inv * x1 + M_inv * x2 + M_inv * pe,
    -1/(Rp*Tg)*x1 - 1/Tg * x2 + 1/Tg * x3,
    -Ki*x1,
)

hx = x[0]


mhe_prob = MHE(
    name='mhe_prob',
    x=x,
    u=pe,
    p=ca.vertcat(M_inv, D, Rp, Ki),
    f=fx,
    h=hx,
    p_low=ca.SX([0.1, 0, 0.1, 0.1]),
    p_upp=ca.SX([10, 5, 10, 10]),
    arrival_cost="qr",
    solver=Solver(
        N=N, 
        qp_solver='qpoases', 
        qp_solver_options={'tol': 1e-6, 'max_iter': 20},
        ts=Ts, 
        integrator='rk4', 
        nlp_method='ggn', 
        cond_blocksize=N,
        ),
    params=[Tg],
)

mhe_prob.generate(
    filename='mhe_test',
    interface='sfun',
)

# def fc(x, u):
#     ki = 2.0
#     Tg = 0.2
#     Rp = 0.05

#     x1 = x[0]
#     x2 = x[1]
#     x3 = x[2]

#     M = x[3]
#     D = x[4]
#     # M = 4.0
#     # D = 1.5

#     return ca.vertcat(
#         -D/M * x1 + 1/M * x2 + 1/M * u,
#         -1/(Rp*Tg)*x1 - 1/Tg *x2 + 1/Tg * x3,
#         -ki*x1,
#         0,
#         0
#     )


# def Hc(x):
#     return x[0]


# N = 50

# x = StateVars("x", NX)


# w = InputVars("w", 3)

# V = Parameters("V")
# W = Parameters("W", (3, 3))

# y = Parameters("y", 1, N+1)
# u = Parameters("u", 1, N)

# delta = Parameters('delta')


# initial_cost = 0
# stage_cost = ca.vertcat(ca.sqrt(V) @ (y-Hc(x)), ca.sqrt(delta) * x[1:], ca.sqrt(W) @ w)
# terminal_cost = ca.vertcat(ca.sqrt(V) @ (y-Hc(x)), ca.sqrt(delta) * x[1:])


# mhe = DynamicOpt(
#     name='mhe',
#     x=x,
#     u=w,
#     c0=initial_cost,
#     ck=stage_cost,
#     cN=terminal_cost,
#     f=ode_integrate(fc, "rk4", 0.02, x, u)+ca.vertcat(w, 0, 0),
#     x_low=ca.vertcat(-ca.inf, -ca.inf, -ca.inf, 0.1, 0),
#     x_upp=ca.vertcat(ca.inf, ca.inf, ca.inf, 10, 5),
#     x0_low=ca.vertcat(-ca.inf, -ca.inf, -ca.inf, 0.1, 0),
#     x0_upp=ca.vertcat(ca.inf, ca.inf, ca.inf, 10, 5),
#     ext_params=[V, W, u, y, delta],
#     params_priority=[0, 0, 1, 1, 0],
#     n_horizon=N,
#     use_ggn=True
# )


# mhe.export_code("qpoases",
#            outdir='gen_code',
#            print_help=True,
#            output_indices={
#                   "x_hat": slice(-NX, -2),
#                   "p_hat": slice(-2, 0)},
#            mex=True,
#            options={
#                'maxIter': 5,
#                }
#            )
