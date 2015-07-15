# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
This is an example driver for the simulation. It should soon be refactored to
result in an input file, an input parser, a solver interface, and output
scripts.
"""

import numpy as np
from scipy.integrate import ode
import importlib

from scipy.integrate._ode import vode
class my_vode(vode):
    def step(self, *args):
        itask = self.call_args[2]
        self.rwork[0] = args[4]
        self.call_args[2] = 5
        r = self.run(*args)
        self.call_args[2] = itask
        return r


from utils.logger import logger
from inp import sim_info
from ur import units
from utils import plotter

from th_component import THSuperComponent

np.set_printoptions(precision=5, threshold=np.inf)

infile = importlib.import_module("input")

si = sim_info.SimInfo(timer=infile.ti,
                      components=infile.components,
                      iso=infile.fission_iso,
                      e=infile.spectrum,
                      n_precursors=infile.n_pg,
                      n_decay=infile.n_dg,
                      kappa=infile.kappa,
                      rho_ext=infile.rho_ext,
                      feedback=infile.feedback,
<<<<<<< HEAD
                      rho_ext=infile.rho_ext,
                      plot_dir=infile.plot_dir
                      )
=======
                      output_plot_dir=infile.output_plot_dir)
>>>>>>> f8c62ad56461d52f7449f41ce3dde677dde15739

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
    _y[t_idx][:n_n] = y_n


def update_th(t, y_n, y_th):
    """This function updates the thermal hydraulics block.

    :param t: the time [s] at which the update is occuring.
    :type t: float.
    :param y_th: The array that solves thermal hydraulics block at time t
    :type y_th: np.ndarray.
    """
    t_idx = si.timer.t_idx(t*units.seconds)
    for idx, comp in enumerate(si.components):
        #if isinstance(comp, THSuperComponent):
        #    comp.update_temp_R(t_idx, y_th[idx+1]*units.kelvin, y_th[idx-1]*units.kelvin)
        #else:
        comp.update_temp(t_idx, y_th[idx]*units.kelvin)
    n_n = len(y_n)
    _y[t_idx][n_n:] = y_th


def update_f(t, y):
    """ update f by updating n(neutronics) and th(thermal-hydraulics) arrays
    """
    #print 'update_f at %f' %t
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


def f_th(t, y_th):
    """Returns the thermal hydraulics solution at time t

    :param t: the time [s] at which the update is occuring.
    :type t: float.
    :param y: TODO
    :type y: np.ndarray
    """
    #print 'time in f_th %f' %t
    t_idx = si.timer.t_idx(t*units.seconds)
    f = units.Quantity(np.zeros(shape=(n_components,), dtype=float),
                       'kelvin / second')
    power = _y[t_idx][0]
    #print 'power %f' %power
    o_i = 1+si.n_pg
    o_f = 1+si.n_pg+si.n_dg
    omegas = _y[t_idx][o_i:o_f]
    for idx, comp in enumerate(si.components):
        f[idx] = si.th.dtempdt(component=comp,
                               power=power,
                               omegas=omegas,
                               t_idx=t_idx)
    return f


def f(t, y):
    print 't in f %f' %t
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
    _y[0] = f
    return f


def y0_n():
    """The initial conditions for y_n, the neutronics sub-block of y"""
    idx = si.n_pg+si.n_dg + 1
    y = y0()[:idx]
    return y


def y0_th():
    """The initial conditions for y_th, the thermal hydraulics sub-block of
    y"""
    tidx = si.n_pg+si.n_dg + 1
    y = y0()[tidx:]
    return y


def solve():
    """Conducts the solution step, based on the dopri5 integrator in scipy"""
<<<<<<< HEAD
    #eqn = ode(f).set_integrator('vode', method='bdf', nsteps=infile.nsteps, max_step=1.0)
    
    eqn = ode(f).set_integrator('dopri5', nsteps=infile.nsteps)
    #eqn = ode(f)
    #eqn._integrator= my_vode(method='bdf', order=2, nsteps=infile.nsteps, max_step=1.0)
=======
    eqn = ode(f).set_integrator('vode', method='bdf', nsteps=infile.nsteps, max_step=1.0, order=2)
    #eqn = ode(f).set_integrator('dopri5', nsteps=infile.nsteps)
>>>>>>> f8c62ad56461d52f7449f41ce3dde677dde15739
    eqn.set_initial_value(y0(), si.timer.t0.magnitude)
    tf1=40*units.seconds
    while (eqn.successful() and eqn.t < tf1.magnitude): #si.timer.tf1.magnitude):
      #TODO: change eqn.t limit to input
<<<<<<< HEAD
        #print 'before'
        #print _y
        #si.timer.advance_one_timestep()
        #print 'mid'
        #print _y
        #eqn.integrate(si.timer.current_time().magnitude)
        eqn.integrate(tf1.magnitude)
        #assert eqn.t+0.01>si.timer.current_time().magnitude, '%f and %f' %(eqn.t, 
        #    si.timer.current_time().magnitude)
        #eqn.integrate(si.timer.tf.magnitude, step=True)
        print 'timer time %f' %si.timer.current_time().magnitude
        print 'eqn time %f' %eqn.t
        #update_f(eqn.t, eqn.y)
        #print 'after'
        #print _y
    #eqn_trans = ode(f)
    #eqn_trans._integrator= my_vode(method='bdf', nsteps=infile.nsteps*10, max_step=1.0)
    eqn_trans = ode(f).set_integrator('dopri5', nsteps=infile.nsteps)
    #eqn_trans = ode(f).set_integrator('vode', method='bdf', nsteps=infile.nsteps, max_step=1.0)
    eqn_trans.set_initial_value(eqn.y, si.timer.t0.magnitude)
=======
        si.timer.advance_one_timestep()
        eqn.integrate(si.timer.current_time().magnitude)
        #eqn.integrate(si.timer.tf.magnitude, step=True)
        print si.timer.current_time().magnitude
        update_f(eqn.t, eqn.y)
    eqn_trans = ode(f).set_integrator('vode', method='bdf', nsteps=infile.nsteps*10, max_step=1.0, order=2)
    eqn_trans.set_initial_value(eqn.y, eqn.t)
>>>>>>> f8c62ad56461d52f7449f41ce3dde677dde15739
    while (eqn_trans.successful() and eqn_trans.t < si.timer.tf.magnitude):
        #print 'before'
        #print _y
        si.timer.advance_one_timestep()
        #print 'mid'
        #print _y
        eqn_trans.integrate(si.timer.current_time().magnitude)
        #eqn_trans.integrate(si.timer.current_time().magnitude, step=True)
        #eqn_trans.integrate(si.timer.tf.magnitude, step=True)
        print 'timer time %f' %si.timer.current_time().magnitude
        update_f(eqn_trans.t, eqn_trans.y)
        #print 'after'
        #print _y
    return _y


def log_results():
    logger.info("\nReactivity : \n"+str(si.ne._rho))
    logger.info("\nFinal Result : \n"+np.array_str(_y))
    for comp in si.components:
        logger.info("\n" + comp.name + ":\n" + np.array_str(comp.T.magnitude))
    logger.info("\nPrecursor lambdas: \n"+str(si.ne._pd.lambdas()))
    logger.info("\nDelayed neutron frac: \n"+str(si.ne._pd.beta()))
    logger.info("\nPrecursor betas: \n"+str(si.ne._pd.betas()))
    logger.info("\nDecay kappas: \n"+str(si.ne._dd.kappas()))
    logger.info("\nDecay lambdas: \n"+str(si.ne._dd.lambdas()))


"""Run it as a script"""
if __name__ == "__main__":
    with open('logo.txt', 'r') as logo:
        logger.critical("\nWelcome to PyRK.\n" +
                        "(c) Kathryn D. Huff\n" +
                        "Your simulation is starting.\n" +
                        "Perhaps it's time for a coffee.\n" +
                        logo.read())
    sol = solve()
    log_results()
    plotter.plot(sol, si, si.plot_dir)
    logger.critical("\nSimulation succeeded.\n")
