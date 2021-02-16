from Reactor_ODE_Model import *
from Reactor_Design_Function_File import *

#Catalyst Info:
d_catalyst = 0.004 #m

#Feed Info:
Mass_Toluene_Feed= 270; Mass_Aqueous_Feed= 530 #kg/hr
Density_Toluene_Feed = 867; 
Viscosity_Toluene = 0.56 * 10 ** (-3);


#Reactor Info:
diameter_ratio = 2 #min_catalyst diameter * diameter ratio
Number_tubes = 20
Initial_Voidage = 0.7
length = 10

"""Defining essential model parameters"""
modelparameters = modelparameters( d_catalyst, Mass_Toluene_Feed, Mass_Aqueous_Feed, diameter_ratio
                                   , Number_tubes, Initial_Voidage, length, Density_Toluene_Feed,Viscosity_Toluene)

"""Simulatiing the model """
results, result_list = simulate_model(modelparameters = modelparameters)

"""Plotting the model results """
plot_simulation_results(results)