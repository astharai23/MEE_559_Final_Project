import casadi as ca
from dynopt import Constraints, Solver, Parameter, Problem
from dynopt.mhe import MHE


N = 50
Ts = 0.02

x = ca.SX.sym('x', 2)
u = ca.SX.sym('u', 1)

M = 4
D = 1.5
Rp = 0.05
Tg = 0.2
Ki = 2


w = x[0]
wd = x[1]


pinv = u[0]

x0 = Parameter("x0", 2)
Pd = Parameter("pd")
Q11 = Parameter("Q11")
Q22 = Parameter("Q22")
R = Parameter("R")
umax = Parameter("umax")

g1 = Constraints(x0=x0, lbu=-umax, ubu=umax)

solver_s1 = Solver(N=N, qp_solver="qpoases", ts = Ts, integrator="rk4", nlp_method="sqp", cond_blocksize=N)


cost_func = w.T @  Q11 @ w + wd @  Q22@  wd + u @ R@  u


terminal_cost = w* Q11* w + wd * Q22* wd

fx = ca.vertcat(
    wd, 
    - (D/(M*Tg)+ 1/(Rp*M*Tg))*w - ((D/M) +(1/Tg))*wd - 1/(M*Tg) * (Pd - u), 
)



mpc = Problem(name = "MPC_simplified", x=x, u=u, z_alg=None, c0=None, ck=cost_func, cN=terminal_cost, 
               f=fx, constraints=g1, params=[x0,  Pd, Q11, Q22, R, umax], solver=solver_s1)

mpc.generate("MPC_simplified_qpoases", output_indices={"u0": 2}, interface="sfun")