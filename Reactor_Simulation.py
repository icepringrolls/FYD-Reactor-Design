from Reactor_ODE_Model import *

"""Defining essential model parameters"""
modelparameters = modelparameters( d_catalyst = 0.004, Mass_Toluene_Feed= 270, Mass_Aqueous_Feed= 530, diameter_ratio = 2
                                   , Number_tubes = 20, Initial_Voidage = 0.7, length = 5, feed_temperature = 333 )

"""Simulatiing the model """
results, result_list = simulate_model(modelparameters = modelparameters)

"""Plotting the model results """
plot_simulation_results(results)