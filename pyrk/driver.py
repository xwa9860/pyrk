# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is an example driver for the simulation. It should soon be refactored to
result in an input file, an input parser, a solver interface, and output
scripts.
"""

import numpy as np
from scipy.integrate import ode
import importlib
import argparse
from utils import logger
from utils.logger import pyrklog
from inp import sim_info
from ur import units
from utils import plotter


np.set_printoptions(precision=5, threshold=np.inf)

infile = importlib.import_module("input")

si = sim_info.SimInfo(timer=infile.ti,
                      components=infile.components,
                      iso=infile.fission_iso,
                      e=infile.spectrum,
                      n_precursors=infile.n_pg,
                      n_decay=infile.n_dg,
                      kappa=infile.kappa,
                      feedback=infile.feedback,
                      #rho_ext=infile.rho_ext
                      )

n_components = len(si.components)

# _y is the matrix of dimension timesteps*nb of equations(unknowns)
_y = np.zeros(shape=(si.timer.timesteps(), si.n_entries()), dtype=float)


def update_n(t, y_n):
    """This function updates the neutronics block.

    :param t: the time [s] at which the update is occuring.
    :type t: float.
    :param y_n: The array that solves the neutronics block at time t
    :type y_n: np.ndarray.
    """
    t_idx = si.timer.t_idx(t*units.seconds)
    n_n = len(y_n)
    si.y[t_idx][:n_n] = y_n


def update_th(t, y_n, y_th, si):
    """This function updates the thermal hydraulics block.

    :param t: the time [s] at which the update is occuring.
    :type t: float.

    :param y_th: The array that solves thermal hydraulics block at time t
    :type y_th: np.ndarray.
    """
    t_idx = si.timer.t_idx(t*units.seconds)
    for idx, comp in enumerate(si.components):
        comp.update_temp(t_idx, y_th[idx]*units.kelvin)
    n_n = len(y_n)
    si.y[t_idx][n_n:] = y_th


def update_f(t, y):
    """ update f by updating n(neutronics) and th(thermal-hydraulics) arrays
    """
    idx = 1+si.n_pg+si.n_dg
    y_n = y[:idx]
    y_th = y[idx:]
    update_n(t, y_n)
    update_th(t, y_n, y_th)


def f_n(t, y):
    """Returns the neutronics block solution at time t

    :param t: the time [s] at which the update is occuring.
    :type t: float.
    :param y: TODO
    :type y: np.ndarray
    """
    n_n = 1 + si.n_pg + si.n_dg
    end_pg = 1 + si.n_pg
    f = np.zeros(shape=(n_n,), dtype=float)
    i = 0
    f[i] = si.ne.dpdt(si.timer.ts, si.components, y[0], y[1:end_pg])
    for j in range(0, si.n_pg):
        i += 1
        f[i] = si.ne.dzetadt(t, y[0], y[i], j)
    assert(i == end_pg-1)
    for k in range(0, si.n_dg):
        i += 1
        f[i] = si.ne.dwdt(y[0], y[i], k)
    return f


def f_th(t, y_th, si):
    """Returns the thermal hydraulics solution at time t

    :param t: the time [s] at which the update is occuring.
    :type t: float.
    :param y: TODO
    :type y: np.ndarray
    """
    t_idx = si.timer.t_idx(t*units.seconds)
    f = units.Quantity(np.zeros(shape=(si.n_components(),), dtype=float),
                       'kelvin / second')
    power = si.y[t_idx][0]
    #power = _y[t_idx][0]
    o_i = 1+si.n_pg
    o_f = 1+si.n_pg+si.n_dg
    omegas = si.y[t_idx][o_i:o_f]
    for idx, comp in enumerate(si.components):
        f[idx] = si.th.dtempdt(component=comp,
                               power=power,
                               omegas=omegas,
                               t_idx=t_idx)
    return f


def f(t, y):
    i_th = 1+si.n_pg+si.n_dg
    y_th = y[i_th:]
    to_ret = np.concatenate((f_n(t, y), f_th(t, y_th)))
    return to_ret


def y0():
    """The initial conditions for y"""
    i = 0
    f = np.zeros(shape=(si.n_entries(),), dtype=float)
    f[i] = 1 # todo change the code to have this from input
    # real power is 236 MWth, but normalized is 1
    for j in range(0, si.n_pg):
        i += 1
        f[i] = f[0]*si.ne._pd.betas()[j]/(si.ne._pd.lambdas()[j]*si.ne._pd.Lambda())
    for k in range(0, si.n_dg):
        i += 1
        f[i] = 0
    for idx, comp in enumerate(si.components):
        f[i+idx+1] = comp.T0.magnitude
    assert len(f) == si.n_entries()
    si.y[0] = f
    return f


def y0_n(si):
    """The initial conditions for y_n, the neutronics sub-block of y"""
    idx = si.n_pg+si.n_dg + 1
    y = y0(si)[:idx]
    return y


def y0_th(si):
    """The initial conditions for y_th, the thermal hydraulics sub-block of
    y"""
    thidx = si.n_pg+si.n_dg + 1
    y = y0(si)[thidx:]
    return y


def solve(si, y, infile):
    """Conducts the solution step, based on the dopri5 integrator in scipy"""
    eqn = ode(f).set_integrator('dopri5', nsteps=infile.nsteps)
    eqn.set_initial_value(y0(), si.timer.t0.magnitude)
    while (eqn.successful() and eqn.t < si.timer.tf.magnitude):
        si.timer.advance_one_timestep()
        eqn.integrate(si.timer.current_time().magnitude)
        update_f(eqn.t, eqn.y)
    return _y

def post_proc():
    '''solution from the equation systems may not be temperature, for example in
    spherical system, the equations are solved for U=rT, where r and T are
    radius and temperature of the point where the equation is solved for.
    '''
    pass

def log_results(si):
    pyrklog.info("\nReactivity : \n"+str(si.ne._rho))
    pyrklog.info("\nFinal Result : \n"+np.array_str(si.y))
    for comp in si.components:
        pyrklog.info("\n" + comp.name + ":\n" + np.array_str(comp.T.magnitude))
    pyrklog.info("\nPrecursor lambdas: \n"+str(si.ne._pd.lambdas()))
    pyrklog.info("\nDelayed neutron frac: \n"+str(si.ne._pd.beta()))
    pyrklog.info("\nPrecursor betas: \n"+str(si.ne._pd.betas()))
    pyrklog.info("\nDecay kappas: \n"+str(si.ne._dd.kappas()))
    pyrklog.info("\nDecay lambdas: \n"+str(si.ne._dd.lambdas()))


def print_logo():
    with open('logo.txt', 'r') as logo:
        pyrklog.critical("\nWelcome to PyRK.\n" +
                         "(c) Kathryn D. Huff\n" +
                         "Your simulation is starting.\n" +
                         "Perhaps it's time for a coffee.\n" +
                         logo.read())


def main(args):
    np.set_printoptions(precision=5, threshold=np.inf)
    logger.set_up_pyrklog(args.logfile)
    infile = importlib.import_module(args.infile)

    si = sim_info.SimInfo(timer=infile.ti,
                          components=infile.components,
                          iso=infile.fission_iso,
                          e=infile.spectrum,
                          n_precursors=infile.n_pg,
                          n_decay=infile.n_dg,
                          kappa=infile.kappa,
                          feedback=infile.feedback,
                          rho_ext=infile.rho_ext,
                          plotdir=args.plotdir)
    print_logo()
    # n_components = len(si.components)
    sol = solve(si=si, y=si.y, infile=infile)
    log_results(si)
    plotter.plot(sol, si)
    pyrklog.critical("\nSimulation succeeded.\n")


"""Run it as a script"""
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description='PyRK parameters')
    ap.add_argument('--infile', help='the name of the input file',
                    default='input')
    ap.add_argument('--logfile', help='the name of the log file',
                    default='pyrk.log')
    ap.add_argument('--plotdir', help='the name of the directory of output plots',
                    default='images')
    ap.add_argument('--outfile', help='the name of the output database',
                    default='pyrk.h5')
    args = ap.parse_args()
    main(args)
