from pyomo import dae as pod
from pyomo import environ as po
import math
import numpy as np
from Reactor_Design_Function_File import *
from matplotlib import pyplot as plt

class modelparameters:
    def __init__(self, d_catalyst , Mass_Toluene_Feed, Mass_Aqueous_Feed, diameter_ratio, Number_tubes, Initial_Voidage, length, Density_Toluene_Feed,Viscosity_Toluene ):

        # =================================#
        # Organic Feed Composition:
        self.Mass_Toluene_Feed = Mass_Toluene_Feed  # kg/h
        self.Mass_Fraction_Toluene_Feed = 1  # kg/h
        self.Density_Toluene_Feed = Density_Toluene_Feed  # kg/m3
        self.Viscosity_Toluene = Viscosity_Toluene

        # --------------------------------#
        # Aqueous Feed Composition:
        self.Mass_Aqueous_Feed = Mass_Aqueous_Feed  # kg/h
        self.Mass_Fraction_NA_Aqueous_Feed = 0.7  # kg/h - Remainder = Water
        self.Density_Aqueous_Feed = 1000  # kg/m3
        self.Viscosity_Aq = 8.9 * 10 ** (-4)

        # ===============================#
        self.MW_Toluene = 92.14  # g/mol or kg/kmol
        self.MW_NA = 63.01  # g/mol or kg/kmol

        # ---------------------------------#
        self.d_catalyst = d_catalyst  # m
        self.D_tube_minimum = 8 * self.d_catalyst
        self.D_tube = self.D_tube_minimum * diameter_ratio
        self.length = length
        self.Number_tubes = Number_tubes  # number of tubes in bundle
        self.Initial_volume = 0.007885  # m3
        self.Initial_Voidage = Initial_Voidage  # 70% free space
        self.Voidage = self.Initial_Voidage
        self.Diffusivity_Toluene_in_Water = 8.6 * 10 ** (-12)
        self.Saturated_Toluene_Conc = 515  # mg/L or g/m3
        self.d_pore = 20e-10
        self.tortuosity_particle = 1.3  # https://www.sciencedirect.com/science/article/pii/S0304389405007594
        self.intraparticle_void_fraction = 0.39
        self.frequency_factor = po.exp(62.63)  # M-1s-1
        self.activation_energy = 22830
        self.x_A = 0.7  # conversion of toluene
        self.Sb = 1.0055
        self.gas_constant = 8.314  # J/k/mol
        self.Temp = 333  # K

        # ---------------------------------#
        # Intial calculations of model
        self.initial_toluene_conc = Concentration_from_MassComp(self.MW_Toluene, self.Mass_Fraction_Toluene_Feed, self.Density_Toluene_Feed)
        self.initial_nitric_conc = Concentration_from_MassComp(self.MW_NA, self.Mass_Fraction_NA_Aqueous_Feed, self.Density_Aqueous_Feed)
        self.Area_tube = Area_Circle(self.D_tube); self.Combined_Area = self.Area_tube * self.Number_tubes

        self.v_total_hour, self.v_frac_organic, self.v_frac_aq, self.flow_density, self.flow_viscosity = Vol_Flow_proportions(self.Mass_Toluene_Feed,
                                                                                                     self.Density_Toluene_Feed,
                                                                                                     self.Mass_Aqueous_Feed,
                                                                                                     self.Density_Aqueous_Feed,
                                                                                                     self.Viscosity_Toluene,
                                                                                                     self.Viscosity_Aq)  # m3/hr

        self.v_total_second = per_hour_to_per_second(self.v_total_hour)

        self.u_super = Superficial_velocity(self.v_total_second, self.Combined_Area)

        self.Re = Reynolds_J(self.d_catalyst, self.u_super, self.flow_density, self.flow_viscosity)

        self.Re_Check = Reynolds_Assumption_Check(self.Re)
        self.j_d = J_factor_Re_Function(self.Re, self.Voidage);
        self.Sc = Schmidt_Number(self.flow_viscosity, self.flow_density, self.Diffusivity_Toluene_in_Water)
        self.Sh = Sh_number_from_j_factor(self.j_d, self.Re, self.Sc)
        self.k_toluene = MT_coeff_Surface_film_from_Sh(self.Sh, self.Diffusivity_Toluene_in_Water, self.d_catalyst)
        self.Area_Cat, self.Vol_Cat = Area_Volume_Catalyst_Spherical(self.d_catalyst)
        self.a_s = specific_area(self.Voidage, self.Area_Cat, self.Vol_Cat) * self.v_frac_aq

        self.D_p_corrected = Diffusion_coeff_pore(self.d_pore, self.gas_constant, self.Temp, self.MW_Toluene)
        self.D_ea = Get_Effective_Diffusion_constant(self.intraparticle_void_fraction, self.D_p_corrected,
                                                     self.tortuosity_particle)
        self.Bi = Biot_number(self.k_toluene, self.d_catalyst / 2, self.D_ea)
        self.product_selectivity = 1.05

        self.n_A0 = Get_Molar_Flowrate(per_hour_to_per_second(self.Mass_Toluene_Feed) * self.Mass_Fraction_Toluene_Feed,
                                  self.MW_Toluene / 1000)
        self.n_B0 = Get_Molar_Flowrate(per_hour_to_per_second(self.Mass_Aqueous_Feed) * self.Mass_Fraction_NA_Aqueous_Feed,
                                  self.MW_NA / 1000)

        self.C_A0 = initial_conc_overall(self.n_A0, self.v_total_second)
        self.C_B0 = initial_conc_overall(self.n_B0, self.v_total_second)



def create_mdoel(modelparameters):
    model = po.ConcreteModel()

    """"defining fixed parameters"""
    model.biot_number = po.Param(initialize = modelparameters.Bi)
    model.frequency_factor = po.Param(initialize = modelparameters.frequency_factor)
    model.activation_energy = po.Param(initialize = modelparameters.activation_energy)
    model.d_catalyst = po.Param(initialize = modelparameters.d_catalyst)
    model.product_selectivity = po.Param(initialize = modelparameters.product_selectivity)
    model.u_super = po.Param(initialize = modelparameters.u_super)
    model.D_ea = po.Param(initialize = modelparameters.D_ea)
    model.temperature = po.Param(initialize= modelparameters.Temp)
    model.intrinsic_rate_constant = po.Param(initialize = modelparameters.frequency_factor * po.exp(- modelparameters.activation_energy/modelparameters.Temp))
    #model.a_s = po.Param(initialize = modelparameters.a_s)
    model.v_frac_aq = po.Param(initialize = modelparameters.v_frac_aq)


    """defining model variables"""
    model.z = pod.ContinuousSet(bounds = (0,modelparameters.length))
    model.global_effectiveness_factor = po.Var(model.z)
    model.i = po.Set(initialize = ["Nitric", "Toluene"])
    model.reaction_rate = po.Var(model.z)
    model.C = po.Var(model.z, model.i)
    model.P = po.Var(model.z)
    model.dC_dz = pod.DerivativeVar(model.C , wrt = model.z)
    model.dP_dz = pod.DerivativeVar(model.P, wrt = model.z)


    """intial conditions"""
    model.C[0, "Nitric"].fix(modelparameters.C_B0)
    model.C[0, "Toluene"].fix(modelparameters.C_A0)
    model.P[0].fix(0)


    """def _calculate_intrinsic_rate_constant(m,z):
        return m.intrinsic_rate_constant[z] == m.frequency_factor * po.exp(-m.activation_energy / m.temperature)
    model.calculate_intrinsic_rate_constant = po.Constraint(model.z, rule = _calculate_intrinsic_rate_constant)"""

    def _calcualte_global_effectiveness_factor(m,z):
        catalyst_rate_constant = m.intrinsic_rate_constant * m.C[z, "Nitric"]
        thielemodulus = Thiele_Modulus(Diameter_Particle = m.d_catalyst, Intrinsic_Rate = catalyst_rate_constant, D_ea = m.D_ea)
        effectiveness_factor = po.tanh(thielemodulus)/thielemodulus
        return m.global_effectiveness_factor[z] == ((1/effectiveness_factor)+((thielemodulus**2)/m.biot_number))**(-1)
    model.calculate_global_effectiveness_factor = po.Constraint(model.z, rule = _calcualte_global_effectiveness_factor)

    def _calculate_reaction_rate(m,z):
        return m.reaction_rate[z] ==  m.v_frac_aq * m.global_effectiveness_factor[z] * m.intrinsic_rate_constant * m.C[z,"Nitric"] * m.C[z,"Toluene"]
    model.calculate_reaction_rate = po.Constraint(model.z, rule = _calculate_reaction_rate)

    def _material_balance_reactants(m,z,i):
        return m.dC_dz[z,i] == -(1/m.u_super)*(m.reaction_rate[z])*(m.product_selectivity)
    model.material_balance_nitric = po.Constraint(model.z, model.i , rule = _material_balance_reactants)

    def _material_balance_products(m,z):
        return m.dP_dz[z] == (1/m.u_super)*(m.reaction_rate[z])
    model.material_balance_products = po.Constraint(model.z, rule = _material_balance_products)

    return model


def simulate_model(modelparameters):
    """simulating"""
    model = create_mdoel(modelparameters = modelparameters )
    simulator = pod.Simulator(model, package = 'casadi')
    results = simulator.simulate(integrator = 'idas')
    discretizer = po.TransformationFactory("dae.collocation")
    discretizer.apply_to(model, nfe = 10, ncp = 3, scheme= "LAGRANGE-RADAU")
    simulator.initialize_model()
    result_list = [model.C[z,"Nitric"].value for z in model.z]
    return results, result_list

"""simulating and plotting results"""
def plot_simulation_results(results):
    resultmatrix = results[1]
    length_index = results[0]
    nitric_conc = resultmatrix[:,0]
    toluene_conc = resultmatrix[:,1]
    product_conc = resultmatrix[:,2]
    global_effectiveness_factor = resultmatrix[:,3]
    reaction_rate = resultmatrix[:,4]

    plt.figure(1)
    plt.plot(length_index,nitric_conc)
    plt.title("Nitric Acid Concentration WRT Length")

    plt.figure(2)
    plt.plot(length_index,toluene_conc)
    plt.title("Toluene Concentration WRT Length")

    plt.figure(3)
    plt.plot(length_index,product_conc)
    plt.title("Product Concentration WRT Length")

    plt.figure(4)
    plt.plot(length_index ,global_effectiveness_factor)
    plt.title("Global Effectiveness Factor WRT Length")

    plt.figure(5)
    plt.plot(length_index , reaction_rate)
    plt.title("Reaction Rate WRT Length")

    conversion = 1 - (toluene_conc/toluene_conc[0])
    plt.figure(6)
    plt.plot(length_index, conversion)
    plt.title("Toluene Conversion WRT Length")
    plt.show()
