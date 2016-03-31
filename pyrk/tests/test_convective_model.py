from nose.tools import assert_equal, assert_almost_equal, assert_true, \
    assert_false, assert_raises, assert_is_instance, with_setup
from utilities.ur import units
from convective_model import ConvectiveModel
from materials.material import Material

def test_constant_model():
    h_constant = ConvectiveModel(
        20*units.W/units.meter**2/units.kelvin)
    h1_constant = ConvectiveModel(
        20*units.W/units.centimeter**2/units.kelvin)
    assert_equal(
        h_constant.h0,
        20*units.W/units.meter**2/units.kelvin)
    assert_equal(
        h1_constant.h0,
        200000*units.W/units.meter**2/units.kelvin)

def test_wakao_model():
    mat=Material(k=1*units.watt/units.meter/units.kelvin,
                 cp=1*units.joule/units.kg/units.kelvin,
                 mu=2*units.pascal*units.second
                 )
    h_wakao = ConvectiveModel(mat=mat,
                              m_flow=1*units.kg/units.g,
                              a_flow=1*units.meter**2,
                              length_scale=1*units.meter,
                              model='wakao')
    assert_equal(h_wakao.mu,
                 2*units.pascal*units.second)
    rho = 100*units.kg/units.meter**3
    assert_equal(h_wakao.h(rho, 0*units.pascal*units.second),
                 h_wakao.h(rho, 2*units.pascal*units.second))
