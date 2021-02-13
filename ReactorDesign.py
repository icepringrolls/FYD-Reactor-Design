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
Mass_Aqueous_Feed = 500 #kg/h
Mass_Fraction_NA_Aqueous_Feed = 0.7 #kg/h - Remainder = Water
Density_Aqueous_Feed = 1420 #kg/m3
Viscosity_Aq = 8.9*10**(-4)
#===============================#
MW_Toluene = 92.14 #g/mol or kg/kmol
MW_NA = 63.01 #g/mol or kg/kmol
#---------------------------------#
d_catalyst = 0.003 #m
D_tube_minimum = 8*d_catalyst
D_tube = D_tube_minimum*1.25
Number_tubes = 50 #number of tubes in bundle
Initial_volume = 0.007885 #m3
Initial_Voidage = 0.3 #30% free space
Voidage = Initial_Voidage
Diffusivity_Toluene_in_Water = 8.6*10**(-12)
Saturated_Toluene_Conc = 515 #mg/L or g/m3
d_pore = 20e-10
tortuosity_particle = 1.3  #https://www.sciencedirect.com/science/article/pii/S0304389405007594
intraparticle_void_fraction = 0.39
intrinsic_rate_coeff = 1.1e-5 #M-1s-1
x_A = 0.7 #conversion of toluene
Sb = 1.0055
gas_constant = 8.314 #J/k/mol
Temp = 333 #K
Heat_Reaction_1 = -129000 #J/mol
Conductivity = 3.63 #W/mK
#================================================≠#
#Calculations Below:
#================================================#
#Calculate Concentrations in the feeds --> Know what to buy:
if Mass_Fraction_Toluene_Feed==1: print("Toluene Feed is PURE") 
else: print("Toluene Feed is NOT pure");
print("Concentration of Toluene Feed (mol/m3) = ",Concentration_from_MassComp(MW_Toluene,Mass_Fraction_Toluene_Feed,Density_Toluene_Feed)) 
Conc_NA = Concentration_from_MassComp(MW_NA,Mass_Fraction_NA_Aqueous_Feed,Density_Aqueous_Feed)
print("Concentration of HNO3 Feed (mol/m3) = ",Conc_NA)
conc_Tol_Aq = Saturated_Toluene_Conc/MW_Toluene #Calculate equivalent molar concentration of aqueous saturated Toluene
print("Concentration of Toluene in Aqueous phase (mol/m3) = ",conc_Tol_Aq)
print("Relative driving diffusion driving force across organic-aqueous phase (mol/m3) = ",Concentration_from_MassComp(MW_Toluene,Mass_Fraction_Toluene_Feed,Density_Toluene_Feed)-conc_Tol_Aq)
print("------------------------------------------------------------")

#Assess individual tube diameter and area + Combined:
Area_tube = Area_Circle(D_tube);Combined_Area = Area_tube*Number_tubes 
print("Diameter of a tube (m) = ", D_tube); print("Minimum diameter of a tube (m) = ", D_tube_minimum)
if D_tube >= D_tube_minimum: print("Minimum tube diameter SATISFIED")
else: print("Minimum tube diameter INVALID");
print("Tube Inner Cross-sectional area (m2) = ",Area_tube);print("Total Reactor Cross-sectional area (m2) = ",Combined_Area)
print("------------------------------------------------------------")

#Calculate Area & Volume of Catalyst for further use: Assuming spherical
Area_Cat,Vol_Cat = Area_Volume_Catalyst_Spherical(d_catalyst)
a_s = specific_area(Voidage,Area_Cat,Vol_Cat)
print("Volume of catalyst Bead = ",Vol_Cat);print("Surface Area of catalyst Bead = ",Area_Cat);print("Specific Surface area (m2 Cat/m3 Reactor)= ",a_s)
print("------------------------------------------------------------")

#Calculate flow profile, composition nad features
v_total_hour,v_frac_organic,v_frac_aq,flow_density, flow_viscosity = Vol_Flow_proportions(Mass_Toluene_Feed,Density_Toluene_Feed,Mass_Aqueous_Feed,Density_Aqueous_Feed,Viscosity_Toluene,Viscosity_Aq) #m3/hr
v_total_second =  per_hour_to_per_second(v_total_hour); print("Volumetric Flowrate (m3/s) = ",v_total_second)
print("------------------------------------------------------------")

#Calculate superficial velocity:
u_super = Superficial_velocity(v_total_second,Combined_Area)
print("Superficial velocity (m/s) = ",u_super)
print("------------------------------------------------------------")

#J-factor correlation for determining film MT coeff
Re = Reynolds_J(d_catalyst,u_super,flow_density,flow_viscosity); print("Reynolds Number (across sphere) = ",Re); Re_Check = Reynolds_Assumption_Check(Re)
j_d = J_factor_Re_Function(Re,Voidage); 
Sc = Schmidt_Number(flow_viscosity,flow_density,Diffusivity_Toluene_in_Water); print("Schmidt Number = ",Sc)
Sh = Sh_number_from_j_factor(j_d,Re,Sc); print("Sherwood Number = ",Sh)
k_toluene = MT_coeff_Surface_film_from_Sh(Sh,Diffusivity_Toluene_in_Water,d_catalyst); print("Toluene Mass Transfer Coeff (m/s) = ",k_toluene)
print("------------------------------------------------------------")



D_ea = Get_Effective_Diffusion_constant(intraparticle_void_fraction,Diffusivity_Toluene_in_Water,tortuosity_particle)
Bi = Biot_number(k_toluene,d_catalyst/2,D_ea)
D_p_corrected = Diffusion_coeff_pore(d_pore,gas_constant,Temp,MW_Toluene)
TM = Thiele_Modulus(d_catalyst,(intrinsic_rate_coeff/a_s)*Conc_NA,D_p_corrected); print("Thiele Modulus = ",TM)

effectiveness_factor,global_effectiveness_factor = Global_effectiveness_factor(TM,Bi)

n_A0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Toluene_Feed)*Mass_Fraction_Toluene_Feed,MW_Toluene/1000)
n_B0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Aqueous_Feed)*Mass_Fraction_NA_Aqueous_Feed,MW_NA/1000)

Volume = Volume_Calc(v_total_second,global_effectiveness_factor,intrinsic_rate_coeff,Sb,n_A0,n_B0,x_A)
Length_tube = Volume/Combined_Area
Axial_Dispersion_Check(Length_tube,D_tube)
print("Length of tube = ",Length_tube)


print("------------------------------------------------------------")
C_A_surface  = conc_Tol_Aq*(1-(global_effectiveness_factor*intrinsic_rate_coeff*Conc_NA)/(k_toluene*a_s))
#ALTERNATE REACTOR SIZE -> Assuming Constant concentrations of Toluene & HNO3:
V_R_Const = n_A0*(1-x_A)/(global_effectiveness_factor*intrinsic_rate_coeff*Conc_NA*conc_Tol_Aq)
print("Assuming Constant concentrations of Toluene & HNO3: Volume = ",V_R_Const)
Length_tube_constant = V_R_Const/Combined_Area
print("Length of tube = ",Length_tube_constant);Axial_Dispersion_Check(Length_tube_constant,D_tube)
print("------------------------------------------------------------")
Beta = Prater_number(Heat_Reaction_1,D_ea,C_A_surface,Conductivity,Temp)
print("------------------------------------------------------------")
k_reaction = (global_effectiveness_factor*intrinsic_rate_coeff*Conc_NA)/a_s
k_MT = k_toluene
print("Coeff of Reaction: ",k_reaction);print("Coeff of MT Toluene: ",k_MT);
if k_MT*10<k_reaction:
    print("Mass Transfer Limited")
else: print("Not Mass Transfer limited")
Weisz_Platz(effectiveness_factor,TM)
print("------------------------------------------------------------")