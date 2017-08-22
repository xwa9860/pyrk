'''
updated on Aug 9, 2017 by Xin Wang
changed from the initial input file
to represent a simplified model for the FHR core

The simulation has 3 stages:
- initial temperature without feedback
- turn on feedback
- turn on external reactivity
'''
from utilities.ur import units
import th_component as th
import math
from materials.material import Material
from materials.liquid_material import LiquidMaterial
from convective_model import ConvectiveModel
from density_model import DensityModel
from timer import Timer
import numpy as np
#############################################
#
# User Workspace
#
#############################################

# Simulation parameters
# Initial time
t0 = 0.00*units.seconds
# Timestep
dt = 0.02*units.seconds
# Final Time
tf = 100.0*units.seconds
# Time to turn on feedback
t_feedback = 40.0*units.seconds

# Thermal hydraulic params
# Temperature feedbacks of reactivity
alpha_fuel = -2.68*units.pcm/units.kelvin
alpha_cool = -0.138*units.pcm/units.kelvin

# initial temperature
t_fuel = (691+273.15)*units.kelvin  # from comsol steady state
t_cool = (678+273.15)*units.kelvin  # from comsol steady state

kappa = 0.0


def area_sphere(r):
    assert(r >= 0*units.meter)
    return (4.0)*math.pi*pow(r.to('meter'), 2)


def vol_sphere(r):
    assert(r >= 0*units.meter)
    return (4./3.)*math.pi*pow(r.to('meter'), 3)

# volumes
n_pebbles = 11000
r_fuel = 3.0/100.0*units.meter

vol_fuel = vol_sphere(r_fuel)
vol_cool = (vol_fuel)*0.4/0.6
a_pb = area_sphere(r_fuel)


#############################################
#
# Required Input
#
#############################################

# Total power, Watts, thermal
power_tot = 10000000.0*units.watt

# Timer instance, based on t0, tf, dt
ti = Timer(t0=t0, tf=tf, dt=dt, t_feedback=t_feedback)

# Number of precursor groups
n_pg = 6

# Number of decay heat groups
n_dg = 0

# Fissioning Isotope
fission_iso = "tmsr"
# Spectrum
spectrum = "multipt"

#two-point model
n_ref = 1
Lambda_ref = 0.000226807
ref_lambda = [1045.2433587850091]
ref_rho = [0.31650504848545435]

# Feedbacks, False to turn reactivity feedback off. True otherwise.
feedback = True

# External Reactivity
from reactivity_insertion import StepReactivityInsertion
rho_ext = StepReactivityInsertion(timer=ti,
                                  t_step=t_feedback + 10.0*units.seconds,
                                  rho_init=0.0*units.delta_k,
                                  rho_final=650.0*units.pcm)

# maximum number of internal steps that the ode solver will take
nsteps = 5000

k_fuel = 15*units.watt/(units.meter*units.kelvin)
cp_fuel = 1818.0*units.joule/units.kg/units.kelvin
rho_fuel = DensityModel(a=3220.0*units.kg/(units.meter**3), model="constant")
Fuel = Material('fuel', k_fuel, cp_fuel, rho_fuel)


mu0 = 0*units.pascal*units.second
k_cool = 1.1*units.watt/(units.meter*units.kelvin)
cp_cool = 2386*units.joule/(units.kg*units.kelvin)
rho_cool = DensityModel(a=2413 *
                        units.kg /
                        (units.meter**3), b=-0.488 *
                        units.kg /
                        (units.meter**3) /
                        units.kelvin, model="linear")
cool = LiquidMaterial('cool', k_cool, cp_cool, rho_cool, mu0)

# Coolant flow properties
# 4700TODO implement h(T) model
h_cool = ConvectiveModel(h0=6000.0*units.watt/units.kelvin/units.meter**2,
                         mat=cool,
                         model='constant')
m_flow = 120.0*units.kg/units.seconds
t_inlet = units.Quantity(672.0, units.degC)

fuel = th.THComponent(name="fuel",
                      mat=Fuel,
                      vol=vol_fuel,
                      T0=t_fuel,
                      alpha_temp=alpha_fuel,
                      timer=ti,
                      heatgen=True,
                      power_tot=power_tot/n_pebbles,
                      sph=True,
                      ri=0.0*units.meter,
                      ro=r_fuel)
  
# mesh size for the fuel pebble FVM calculation
l = 0.0005*units.meter
comp_list = fuel.mesh(l)
pebble = th.THSuperComponent('pebble', t_fuel, comp_list, timer=ti)
# Add convective boundary condition to the pebble
pebble.add_conv_bc('cool', h=h_cool)

cool = th.THComponent(name="cool",
                      mat=cool,
                      vol=vol_cool,
                      T0=t_cool,
                      alpha_temp=alpha_cool,
                      timer=ti)
# The coolant convects to the shell
cool.add_convection('pebble', h=h_cool, area=a_pb)
cool.add_advection('cool', m_flow/n_pebbles, t_inlet, cp=cool.cp)

components = []
for i in range(0, len(pebble.sub_comp)):
    components.append(pebble.sub_comp[i])
components.extend([pebble, cool])
