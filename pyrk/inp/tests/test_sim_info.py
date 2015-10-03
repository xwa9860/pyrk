from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance, with_setup

from pyrk.inp import sim_info as si

from utilities.ur import units
import th_system
import th_component
from timer import Timer

def test_init_reasonable_sim():
    t0 = 0*units.seconds
    tf = 10*units.seconds
    dt = 0.1*units.seconds
    ti = Timer(t0=t0,tf=tf,dt=dt)
    iso = "u235"
    spectrum = "thermal"
    npg = 6
    ndg = 11
    kappa = 0.06
    tester = th_component.THComponent()
    kappa = 0.0
    info = si.SimInfo(timer=ti, components={}, iso=iso, e=spectrum,
                      n_precursors=npg, n_decay=ndg, kappa=kappa)
    assert_equal(t0, info.timer.t0)
    assert_equal(tf, info.timer.tf)
    assert_equal(dt, info.timer.dt)
    assert_equal(info.timer.timesteps(), 101)
