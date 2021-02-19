from Reactor_ODE_Model import *
from Reactor_ODE_Model_FLUIDIZED import *
from Reactor_Design_Function_File import *
from Fluidization_Reactor_Functions import *

#Feed Info: NEEDS CHECKING BEFORE FINAL DESIGN
Mass_Toluene_Feed= 270; Mass_Aqueous_Feed= 530 #kg/hr
Mass_Fraction_NA_Aqueous_Feed = 0.7
Density_Toluene_Feed = 867; Density_Aqueous_Feed = 1000 #kg/hr
Viscosity_Toluene = 0.56 * 10 ** (-3); Viscosity_Aq = 8.9 * 10 ** (-4) 

#Conversion:
x_A = 0.7
#Reactor Info:
Temp = 333
diameter_ratio = 2 #min_catalyst diameter * diameter ratio
Number_tubes = 20
Initial_Voidage = 0.7
length = 10
solid_density = 2100
#Catalyst Info:
d_catalyst = 0.004 #m
D_tube_minimum = 8*d_catalyst
D_tube = 0.01; 
d_pore = 20e-10
intraparticle_void_fraction = 0.39
if D_tube>=D_tube_minimum:
    print("Tube minimum is SATISFIED")
else: print("Tube minimum is NOT satisfied")

#Chemical Info:
MW_Toluene = 92.14  # g/mol or kg/kmol
MW_NA = 63.01  # g/mol or kg/kmol
Diffusivity_Toluene_in_Water = 8.6 * 10 ** (-12)
Saturated_Toluene_Conc = 515
tortuosity_particle = 1.3
Sb = 1.0055
Bubble_diameter = 0.01

#============================================================#




#============================================================#



Area_tube = Area_Circle(D_tube);Combined_Area = Area_tube*Number_tubes; equiv_diameter = math.sqrt(4*Combined_Area/math.pi)
v_total_hour,v_frac_organic,v_frac_aq,flow_density, flow_viscosity = Vol_Flow_proportions(Mass_Toluene_Feed,Density_Toluene_Feed,Mass_Aqueous_Feed,Density_Aqueous_Feed,Viscosity_Toluene,Viscosity_Aq)
v_total_second =  per_hour_to_per_second(v_total_hour); print("Volumetric Flowrate (m3/s) = ",v_total_second)
print("------------------------------------------")
u_super = Superficial_velocity(v_total_second,Combined_Area); print("Superficial Velocity (m/s) = ",u_super)
U_mf,Re_fluid_min = Get_U_mf(Initial_Voidage,solid_density,d_catalyst,v_total_second,v_frac_organic,v_frac_aq,flow_density, flow_viscosity)
Re = Reynolds_J(d_catalyst,u_super,flow_density,flow_viscosity); print("Reynolds number (Actual) = ",Re)
C_d = 0
if Re<1:
    C_d = 24/Re
elif Re>=1 and Re<1000:
    C_d = math.exp(-5.5+69.43/(math.log(Re)+7.99))
elif Re>=1000:
    C_d = 0.4
print("Drag Coefficient =",C_d)
U_t = terminal_velocity(C_d,d_catalyst,flow_density,solid_density); print("Terminal Velocity (m/s) = ",U_t)

if u_super>=U_mf:
    Reactor_type = "Fluidized Bed";
elif u_super>=U_t:
    Reactor_type = "Transport Reactor"
else:  Reactor_type = "Fixed Bed";
print("Reactor Type: ",Reactor_type)
print("------------------------------------------")
U_b = bubble_velocity(u_super,U_mf,D_tube,Bubble_diameter); print("Bubble Velocity (m/s) = ",U_b)



print("==========================================")
"""FIXED BED CALCULATIONS"""
if Reactor_type == "Fixed Bed":
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
                                  Viscosity_Aq,D_tube_minimum,D_tube,MW_Toluene,MW_NA,Diffusivity_Toluene_in_Water,
                                  Saturated_Toluene_Conc,d_pore,tortuosity_particle,intraparticle_void_fraction,
                                  x_A,Sb,Temp)
    """Simulatiing the model """
    results, result_list = simulate_model(modelparameters = modelparameters)

    """Plotting the model results """
    plot_simulation_results(results)

if Reactor_type == "Fluidized Bed":
    modelparameters = fluidized_modelparameters( d_catalyst, 
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
                                  Viscosity_Aq,D_tube_minimum,D_tube,MW_Toluene,MW_NA,Diffusivity_Toluene_in_Water,
                                  Saturated_Toluene_Conc,d_pore,tortuosity_particle,intraparticle_void_fraction,
                                  x_A,Sb,Temp,U_mf,Bubble_diameter,solid_density)
    """Simulatiing the model """
    results, result_list = fluidized_simulate_model(modelparameters = modelparameters)

    """Plotting the model results """
    fluidized_plot_simulation_results(results)