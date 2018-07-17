from scipy.optimize import (minimize, NonlinearConstraint, 
                            LinearConstraint, Bounds)
import pycutestmgr as cute
import time
import numpy as np
from numpy.linalg import norm
import scipy.sparse as spc
import warnings
import sys
import os
import json

import scipy

print(scipy.__version__)

check_tol = 1e-3

import psutil
import resource
max_mem = int(0.7 * psutil.virtual_memory().available)
resource.setrlimit(resource.RLIMIT_AS, (max_mem, max_mem))

def print_header():
    hdr = ("|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|{:^8}|"
           .format("name", "nit", "nfev", "success", "cons_ok", "ok_trust_constr"))
    print("="*len(hdr))
    print(hdr)
    print("-"*len(hdr))
    sys.stdout.flush()

def print_problem_sol(*args):
    fmt = "|{:^8}"*len(args) + "|"
    print(fmt.format(*args))
    sys.stdout.flush()

default_options = {'sparse_jacobian':True, 'maxiter':1000, 'xtol':1e-7, 'gtol':1e-7}
list_feasible_box_constr = ["HS13", "HS105", "BROYDNBD"]

def solve_problem(prob):
    name = prob[0]
    cache_file = os.path.join(os.environ['PYCUTEST_CACHE'],
                              'result-' + name + '.json')

    try:
        with open(cache_file, 'r') as f:
            result = json.load(f)
        if len(result) != 6+1:
            raise IOError()
    except (IOError, ValueError):
        result = None

    hash_data = [check_tol, prob[1], prob[2]]
    hash_data = json.loads(json.dumps(hash_data))

    if result is None or result[-1] != hash_data:
        result = list(_solve_problem(prob)) + [hash_data]
        with open(cache_file, 'w') as f:
            json.dump(result, f)

    print_problem_sol(*result[:-1])    
    
def _solve_problem(prob):
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
    cons_eq_mask = (c_ub == c_lb) & (~np.isinf(c_ub))
    
    cons_ub_mask = (~np.isinf(c_ub)) & (~cons_eq_mask)
    cons_lb_mask = (~np.isinf(c_lb)) & (~cons_eq_mask)

    def constr_eq(x):
        c = problem.cons(x)
        a = (c - c_lb)[cons_eq_mask]
        return a

    def constr_onesided(x):
        c = problem.cons(x)
        a = (c - c_lb)[cons_lb_mask]
        b = (c_ub - c)[cons_ub_mask]
        return np.r_[a, b]

    def jac_eq(x):
        _, A1 = problem.cons(x, True)
        return A1[cons_eq_mask,:]

    def jac_onesided(x):
        _, A1 = problem.cons(x, True)
        return np.vstack([A1[cons_lb_mask,:],
                          -A1[cons_ub_mask,:]])

    constr_dict = []

    if cons_eq_mask.any():
        constr_dict.append(dict(type='eq',
                                fun=constr_eq,
                                jac=jac_eq))

    if cons_ub_mask.any() or cons_lb_mask.any():
        constr_dict.append(dict(type='ineq',
                                fun=constr_onesided,
                                jac=jac_onesided))

    box = [(a, b) for a, b in zip(lb, ub)]

    constr_t = NonlinearConstraint(constr,  c_lb, c_ub, jac, lagr_hess)
    box_t = Bounds(lb, ub, name in list_feasible_box_constr)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = minimize(fun, x0, method='SLSQP', jac=grad, hess=hess, constraints=constr_dict, 
                          bounds=box, options=options)

        result0 = minimize(fun, x0, method='trust-constr', jac=grad, hess=hess, constraints=constr_t, 
                          bounds=box_t, options=options)

    # Print Results
    if result.success:
        c = constr_dict[0]['fun'](result.x)
        c0 = constr_dict[0]['fun'](result0.x)
        cons_ok = int((c > -check_tol * (1 + np.linalg.norm(c0))).all())

        if result0.status:
            agree_trust_constr = int(np.linalg.norm(result.x - result0.x) < check_tol * (1 + np.linalg.norm(result0.x)))
        else:
            agree_trust_constr = "--"
    else:
        cons_ok = "--"
        agree_trust_constr = "--"

    return (name, result.nit, result.nfev, int(result.success),
            cons_ok, agree_trust_constr)


ip_problems = [("CORKSCRW", {"T": 50}, {}),  
               ("COSHFUN", {"M": 20}, {"initial_barrier_parameter": 0.1, "initial_tr_radius": 5,"initial_barrier_tolerance":1, "initial_constr_penalty":0.01}),
               ("DIXCHLNV", {"N":100}, {}),
               ("HAGER4", {"N": 100}, {}),
               ("HIMMELBK", {}, {}),
               ("NGONE", {"HNS": 49}, {"initial_tr_radius": 100}), #  many iteractions
               ("OPTCNTRL", {}, {}), 
               ("OPTMASS", {"N": 200}, {}),
               ("ORTHREGF", {"NPTS": 10}, {}),
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
sqp_problems = [("HAGER2", {"N": 100}, {}),
                ("HAGER3", {"N": 100}, {}),
                ("ORTHREGA", {'LEVELS': 3}, {}),
                ("ORTHREGC", {'NPTS': 50}, {}),
                ("ORTHREGD", {'NPTS': 50}, {"initial_tr_radius": 100,"initial_constr_penalty":1}),
                ("DTOC1ND", {'N': 100, 'NX':2, 'NY':4}, {"initial_tr_radius": 1000, "initial_constr_penalty":1}),
                ("DTOC2", {'N': 100, 'NX':2, 'NY':4}, {}),
                ("DTOC3", {'N': 100}, {}),
                ("DTOC4", {'N': 100}, {}),
                ("DTOC5", {'N': 100}, {}),
                ("DTOC6", {'N': 101}, {}),
                ("EIGENA2", {'N': 10}, {}),
                ("EIGENC2", {'M': 2}, {}),
                ("ARTIF", {'N': 100}, {}),
                ("BRATU3D", {'P':10}, {}),
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
