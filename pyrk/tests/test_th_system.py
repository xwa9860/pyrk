from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance, with_setup

import numpy as np
import th_system
import th_component
from utilities.ur import units
from materials.material import Material

def test_dtempfueldt_returns_numbers():
    components = [th_component.THComponent(),
                  th_component.THComponent()]
    T = 750.0
    th = th_system.THSystem(0, components)
    p = 1.0000002
    omegas = np.array([0, 0, 0])

    for c in components:
        obs = th.dtempdt(c, p, omegas, 0)
        assert(obs + T*units.kelvin/units.second > 0*units.kelvin/units.second)

def test_conduction_slab():
    mat = Material(k=1*units.watt/units.meter/units.kelvin)
    components = [th_component.THComponent(mat=mat, T0=700*units.kelvin),
                  th_component.THComponent(mat=mat, T0=700*units.kelvin)]
    th = th_system.THSystem(0, components)
    assert_equal(th.conduction_slab(components[0], components[1], 0,
                              1*units.meter, 1*units.meter**2), 0)
    components = [th_component.THComponent(mat=mat, T0=800*units.kelvin),
                  th_component.THComponent(mat=mat, T0=700*units.kelvin)]
    th = th_system.THSystem(0, components)
    assert(th.conduction_slab(components[0], components[1], 0,
                              1*units.meter, 1*units.meter**2) > 0)
