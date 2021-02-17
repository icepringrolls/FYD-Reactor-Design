from Reactor_ODE_Model import *
from Reactor_Design_Function_File import *


#Feed Info: NEEDS CHECKING BEFORE FINAL DESIGN
Mass_Toluene_Feed= 270; Mass_Aqueous_Feed= 530 #kg/hr
Mass_Fraction_NA_Aqueous_Feed = 0.7
Density_Toluene_Feed = 867; Density_Aqueous_Feed = 1000 #kg/hr
Viscosity_Toluene = 0.56 * 10 ** (-3); Viscosity_Aq = 8.9 * 10 ** (-4) 


#Reactor Info:
diameter_ratio = 2 #min_catalyst diameter * diameter ratio
Number_tubes = 20
Initial_Voidage = 0.7
length = 10
solid_density = 2100
#Catalyst Info:
d_catalyst = 0.004 #m
D_tube_minimum = 8*d_catalyst
D_tube = 0.05; 
if D_tube>=D_tube_minimum:
    print("Tube minimum is SATISFIED")
else: print("Tube minimum is NOT satisfied")

#Chemical Info:
MW_Toluene = 92.14  # g/mol or kg/kmol
MW_NA = 63.01  # g/mol or kg/kmol

#============================================================#
Area_tube = Area_Circle(D_tube);Combined_Area = Area_tube*Number_tubes; equiv_diameter = math.sqrt(4*Combined_Area/math.pi)
v_total_hour,v_frac_organic,v_frac_aq,flow_density, flow_viscosity = Vol_Flow_proportions(Mass_Toluene_Feed,Density_Toluene_Feed,Mass_Aqueous_Feed,Density_Aqueous_Feed,Viscosity_Toluene,Viscosity_Aq)
v_total_second =  per_hour_to_per_second(v_total_hour); print("Volumetric Flowrate (m3/s) = ",v_total_second)
u_super = Superficial_velocity(v_total_second,Combined_Area); print("Superficial Velocity (m/s) = ",u_super)
U_mf = Get_U_mf(Initial_Voidage,solid_density,d_catalyst,v_total_second,v_frac_organic,v_frac_aq,flow_density, flow_viscosity)

print("==========================================")
"""Defining essential model parameters"""
modelparameters = modelparameters( d_catalyst, 
                                  Mass_Toluene_Feed, 
                                  Mass_Aqueous_Feed, 
                                  diameter_ratio, 
                                  Number_tubes, 
                                  Initial_Voidage, 
                                  length, 
                                  Density_Toluene_Feed,
                                  Viscosity_Toluene,
                                  Mass_Fraction_NA_Aqueous_Feed,
                                  Density_Aqueous_Feed,
                                  Viscosity_Aq,D_tube_minimum,D_tube)

"""Simulatiing the model """
results, result_list = simulate_model(modelparameters = modelparameters)

# """Plotting the model results """
# plot_simulation_results(results)