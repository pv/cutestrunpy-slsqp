from scipy.optimize import (minimize, NonlinearConstraint, 
                            LinearConstraint, Bounds)
import pycutestmgr as cute
import time
import numpy as np
from numpy.linalg import norm
import scipy.sparse as spc
import warnings
cute.clearCache('HS71')

def print_header():
    print("|{0:^10}|{1:^5}|{2:^5}|{3:^6}|{4:^10}|{5:^10}|{6:^10}|{7:^10}|{8:^10}|{9:^10}|"
          .format("name", "n", "m", "nnz", "niters", "f evals", "CG iters", "time", "opt", "c viol"))
    s = "-"*9 + ":"
    s1 = "-"*4 + ":"
    s2 =  ":" + "-"*23 + ":"
    s3 = "-"*5 + ":"
    print("|{0:^10}|{1:^5}|{2:^5}|{3:^6}|{4:^10}|{5:^10}|{6:^10}|{7:^10}|{8:^10}|{9:^10}|"
          .format(s, s1, s1, s3, s, s, s, s, s, s, s2, s))

def print_problem_sol(name, n, m, nnz, niters, nfev, cg_niters, time, opt, c_viol, method):
    print("|{0:^10}|{1:^5}|{2:^5}|{3:^6}|{4:^10}|{5:^10}|{6:^10}|{7:^10}|{8:^1.4e}|{9:^1.4e}|"
          .format(name, n, m, nnz, niters, nfev, cg_niters, round(time, 2), opt, c_viol))

default_options = {'sparse_jacobian':True, 'maxiter':1000, 'xtol':1e-7, 'gtol':1e-7}
list_feasible_box_constr = ["HS13", "HS105", "BROYDNBD"]

def solve_problem(prob):
    name = prob[0]
    sifParams = prob[1]
    options = prob[2]        
    options = prob[2].copy()
    options.update(default_options)
    cute.clearCache(name)
    cute.prepareProblem(name, sifParams=sifParams)
    problem = cute.importProblem(name)
    info = problem.getinfo()
    x0 = info["x"]
    n = info["n"]
    m = info["m"]
    nnz = info["nnzj"]
    v0 = info["v"]
    eq_list = info["equatn"]

    def convert_to_inf(bounds):
        new_bounds = np.zeros_like(bounds)
        i = 0
        for b in bounds:
            if b >= 1e20:
                new_bounds[i] = np.inf
            elif b <= -1e20:
                new_bounds[i] = -np.inf
            else:
                new_bounds[i] = b
            i += 1
        return new_bounds

    # Define upper and lower bound
    c_lb = convert_to_inf(info["cl"])
    c_ub = convert_to_inf(info["cu"])
    lb = convert_to_inf(info["bl"])
    ub = convert_to_inf(info["bu"])
    if name in list_feasible_box_constr:
        x0 = np.maximum(lb+1e-10, np.minimum(ub-1e-10, x0))

    # Define function and constraints
    def fun(x):
        return problem.obj(x)[0]

    def grad(x):
        _, g1 = problem.obj(x, True)
        return g1

    def lagr_hess(x, v):
        def matvec(p):
            return problem.hprod(p, x, v)
        return spc.linalg.LinearOperator((n, n), matvec=matvec)


    def constr(x):
        return problem.cons(x)

    def jac(x):
        _, A1 = problem.cons(x, True)
        return spc.csc_matrix(A1)

    def hess(x):
        return spc.csc_matrix((n, n))


    # Constraints
    constr = NonlinearConstraint(constr,  c_lb, c_ub, jac, lagr_hess)
    box = Bounds(lb, ub, name in list_feasible_box_constr)

    start_time = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = minimize(fun, x0, method='trust-constr', jac=grad, hess=hess, constraints=constr, 
                          bounds=box, options=options)
    total_time = time.time() - start_time
    # Print Results
    print_problem_sol(name, n, m, nnz, result.niter, result.nfev,
                      result.cg_niter, total_time, result.optimality, 
                      result.constr_violation, result.method)
    return result

ip_problems = [("CORKSCRW", {"T": 50}, {}),  
               ("COSHFUN", {"M": 20}, {"initial_barrier_parameter": 0.1, "initial_tr_radius": 5,"initial_barrier_tolerance":1, "initial_constr_penalty":0.01}),
               ("DIXCHLNV", {"N":100}, {}),
               ("HAGER4", {"N": 1000}, {}),
               ("HIMMELBK", {}, {}),
               ("NGONE", {"HNS": 49}, {"initial_tr_radius": 100}), #  many iteractions
               ("OPTCNTRL", {}, {}), 
               ("OPTMASS", {"N": 200}, {}),
               ("ORTHREGF", {"NPTS": 20}, {}),
               ("SVANBERG", {"N": 500}, {}), 
               ("READING1", {"N": 100}, {})]

problems_si = [("HS7", {}, {}),
            ("HS10", {}, {}),
            ("HS11", {}, {}),
            ("HS13", {}, {}),
            ("HS14", {}, {}),
            ("HS16", {}, {}),
            ("HS17", {}, {}),
            ("HS19", {}, {}),
            ("HS20", {}, {}),
            ("HS22", {}, {}),
            ("HS24", {}, {}),
            ("HS26", {}, {}),
            ("HS28", {}, {}),
            ("HS31", {}, {}),
            ("HS32", {}, {}),
            ("HS33", {}, {}),
            ("HS39", {}, {}),
            ("HS46", {}, {}),
            ("HS51", {}, {"initial_barrier_parameter": 0.1, "initial_tr_radius": 1,"initial_barrier_tolerance":0.1, "initial_constr_penalty":1}),
            ("HS52", {}, {}),
            ("HS53", {}, {}),
            ("HS63", {}, {}),
            ("HS64", {}, {}),
            ("HS65", {}, {}),
            ("HS70", {}, {}),
            ("HS71", {}, {}),
            ("HS72", {}, {}),
            ("HS73", {}, {}),
            ("HS74", {}, {}),
            ("HS75", {}, {}),
            ("HS77", {}, {}),
            ("HS78", {}, {}),
            ("HS79", {}, {}),
            ("HS80", {}, {}),
            ("HS81", {}, {}),
            ("HS83", {}, {}),
            ("HS84", {}, {}),
            ("HS85", {}, {"initial_barrier_parameter": 0.001}),
            ("HS86", {}, {}),
            ("HS93", {}, {}),
            ("HS95", {}, {}),
            ("HS96", {}, {}),
            ("HS97", {}, {"initial_barrier_parameter": 100}),
            ("HS98", {}, {"initial_barrier_parameter": 100}),
            ("HS99", {}, {}), # Fails (As in the original  paper)
            ("HS100", {}, {}),
            ("HS104", {}, {"initial_barrier_parameter": 100}),
            ("HS105", {}, {}),
            ("HS106", {}, {}),
            ("HS107", {}, {}),
            ("HS108", {}, {}),
            ("HS109", {}, {}),
            ("HS111", {}, {}),
            ("HS112", {}, {"initial_barrier_parameter": 10}),
            ("HS113", {}, {}),
            ("HS114", {}, {"initial_barrier_parameter": 10}),
            ("HS116", {}, {}),
            ("HS117", {}, {}),
            ("HS118", {}, {}),
            ("HS119", {}, {})]
sqp_problems = [("HAGER2", {"N": 5000}, {}),
                ("HAGER3", {"N": 5000}, {}),
                ("ORTHREGA", {'LEVELS': 5}, {}),
                ("ORTHREGC", {'NPTS': 500}, {}),
                ("ORTHREGD", {'NPTS': 5000}, {"initial_tr_radius": 100,"initial_constr_penalty":1}),
                ("DTOC1ND", {'N': 500, 'NX':2, 'NY':4}, {"initial_tr_radius": 1000, "initial_constr_penalty":1}),
                ("DTOC2", {'N': 500, 'NX':2, 'NY':4}, {}),
                ("DTOC3", {}, {}),
                ("DTOC4", {'N': 5000}, {}),
                ("DTOC5", {'N': 5000}, {}),
                ("DTOC6", {'N': 1001}, {}),
                ("EIGENA2", {'N': 50}, {}),
                ("EIGENC2", {'M': 10}, {}),
                ("ARTIF", {'N': 1000}, {}),
                ("BRATU3D", {'P':17}, {}),
                ]

other = [("HS109", {}, {}),]

# Print Header
print_header()
# Print Table
for prob in problems_si:
    solve_problem(prob)

# Print Header
print_header()
# Print Table
for prob in ip_problems:
    solve_problem(prob)

# Print Header
print_header()
# Print Table
for prob in sqp_problems:
    solve_problem(prob)
