#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:26:12 2021

@author: henryorlebar
"""

import numpy as np
import pandas as pd
import math
from Reactor_Design_Function_File import *
from matplotlib import pyplot as plt

#=================================#
#Organic Feed Composition:
Mass_Toluene_Feed = 270 #kg/h
Mass_Fraction_Toluene_Feed = 1 #kg/h
Density_Toluene_Feed = 867 #kg/m3
Viscosity_Toluene = 0.56*10**(-3)
#--------------------------------#
#Aqueous Feed Composition:
Mass_Aqueous_Feed = 700 #kg/h
Mass_Fraction_NA_Aqueous_Feed = 0.7 #kg/h - Remainder = Water
Density_Aqueous_Feed = 1000 #kg/m3
Viscosity_Aq = 8.9*10**(-4)
#===============================#
MW_Toluene = 92.14 #g/mol or kg/kmol
MW_NA = 63.01 #g/mol or kg/kmol
#---------------------------------#
d_catalyst = 0.003 #m
D_tube_minimum = 8*d_catalyst
D_tube = D_tube_minimum
Number_tubes = 3 #number of tubes in bundle
Initial_volume = 0.007885 #m3
Initial_Voidage = 0.7 #70% free space
Voidage = Initial_Voidage
Diffusivity_Toluene_in_Water = 8.6*10**(-12)
Saturated_Toluene_Conc = 515 #mg/L or g/m3
d_pore = 20e-10
tortuosity_particle = 1.3  #https://www.sciencedirect.com/science/article/pii/S0304389405007594
intraparticle_void_fraction = 0.39
intrinsic_rate_coeff = 1.1e-5 #M-1s-1
x_A = 0.7 #conversion of toluene
Sb = 1.0055

#================================================≠#
Area_tube = Area_Circle(D_tube);Combined_Area = Area_tube*Number_tubes

v_total_hour,v_frac_organic,v_frac_aq,flow_density, flow_viscosity = Vol_Flow_proportions(Mass_Toluene_Feed,Density_Toluene_Feed,Mass_Aqueous_Feed,Density_Aqueous_Feed,Viscosity_Toluene,Viscosity_Aq) #m3/hr
v_total_second =  per_hour_to_per_second(v_total_hour)

u_super = Superficial_velocity(v_total_second,Combined_Area)

Re = Reynolds_J(d_catalyst,u_super,flow_density,flow_viscosity); print("Re = ",Re); Re_Check = Reynolds_Assumption_Check(Re)
j_d = J_factor_Re_Function(Re,Voidage); 
Sc = Schmidt_Number(flow_viscosity,flow_density,Diffusivity_Toluene_in_Water)
Sh = Sh_number_from_j_factor(j_d,Re,Sc)
k_toluene = MT_coeff_Surface_film_from_Sh(Sh,Diffusivity_Toluene_in_Water,d_catalyst)
#Add new diffusion coeff eqn
D_ea = Get_Effective_Diffusion_constant(intraparticle_void_fraction,Diffusivity_Toluene_in_Water,tortuosity_particle)
Bi = Biot_number(k_toluene,d_catalyst/2,D_ea)
TM = Thiele_Modulus(d_catalyst,intrinsic_rate_coeff,D_ea)

effectiveness_factor,global_effectiveness_factor = Global_effectiveness_factor(TM,Bi)

n_A0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Toluene_Feed)*Mass_Fraction_Toluene_Feed,MW_Toluene/1000)
n_B0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Aqueous_Feed)*Mass_Fraction_NA_Aqueous_Feed,MW_NA/1000)

Volume = Volume_Calc(v_total_second,global_effectiveness_factor,intrinsic_rate_coeff,Sb,n_A0,n_B0,x_A)
Length_tube = Volume/Combined_Area
print("Length of tube = ",Length_tube)

X = np.linspace(0,x_A,1000)
Vol_x = []; L_x = []

for i in range(0,1000):
    Vol_x.append(Volume_Calc(v_total_second,global_effectiveness_factor,intrinsic_rate_coeff,Sb,n_A0,n_B0,X[i]))
    L_x.append(Vol_x[i]/Combined_Area)


