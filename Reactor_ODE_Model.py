from pyomo import dae as pod
from pyomo import environ as po
import math
import numpy as np
from Reactor_Design_Function_File import *


def create_mdoel(modelparameters):
    model = po.ConcreteModel()




    def _get_global_effectiveness_factor():



    """"defining fixed parameters"""
    model.Bi = modelparameters.Bi


    """defining model variables"""
    model.z = pod.ContinuousSet()
    model.temperature = po.Var()
    model.C_nitric = po.Var()
    model.C_toluene = po.Var()
    model.intrinsic_rate_constant = po.Var()
    model.dN_dz = pod.DerivativeVar(model.C_nitric, wrt=model.z)
    model.dT_dz = pod.DerivativeVar(model.C_toluene, wrt=model.z)
    model.effectiveness_factor = #call global effectiveness function

    def _calculate_intrinsic_rate_constant(m,i):
        return m.intrinsic_rate_constant == m.frequency_factor * math.exp(-m.activation_energy * m.temperature)
    model.calculate_intrinsic_rate_constant = po.Constraint(model.temperature, rule = _calculate_intrinsic_rate_constant)

    def _calcualte_global_effectiveness_factpr(m,i):
        m.catalyst_rate_constant = m.intrisic_rate_constant * m.C_nitric
        m.thielemodulus = Thiele_Modulus(Diameter_Particle = m.d_catalyst)
        m.




    def _material_balance_nitric():


    model.material_balance_nitric = po.Constraint(model.z, rule = _material_balance_nitric)


    def _material_balance_toluene():


    model.material_balance_Toluene = po.Constraint(model.z, rule = *_material_balance_toluene())


    return model


def simulate_model():
    """simulating"""
    model = create_mdoel()
    simulator = pod.Simulator(model, package= 'casadi')
    simulator.simulate(integrator = 'idas')
    simulator.initialize_model()
    discretizer = po.TransformationFactory("dae.collocation")
    discretizer.apply_to(model, nfe = 10, ncp = 3, scheme= "LAGRANGE-RADAU")
