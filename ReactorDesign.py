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

plot_TH_EF_GEF_distribution = 0 #YES=1, NO =0;
#=================================#
#Organic Feed Composition:
Mass_Toluene_Feed = 270 #kg/h
Mass_Fraction_Toluene_Feed = 1 #kg/h
Density_Toluene_Feed = 867 #kg/m3
Viscosity_Toluene = 0.56*10**(-3)
#--------------------------------#
#Aqueous Feed Composition:
Mass_Aqueous_Feed = 212.5 #kg/h
range_set,Aq_feed_test_set,Push_value = TEST_AQUEOUS_MASS_FLOWRATES(1,Mass_Aqueous_Feed,200,1000,100) #YES=1, NO = 0;
Mass_Fraction_NA_Aqueous_Feed = 0.7 #kg/h - Remainder = Water
Density_Aqueous_Feed = 1420 #kg/m3
Viscosity_Aq = 8.9*10**(-4)
#===============================#
MW_Toluene = 92.14 #g/mol or kg/kmol
MW_NA = 63.01 #g/mol or kg/kmol
MW_H20 = 18.02 #g/mol
#---------------------------------#
d_catalyst = 0.003 #m
D_tube_minimum = 8*d_catalyst
D_tube = D_tube_minimum*2
Number_tubes = 50 #number of tubes in bundle
Initial_Voidage = 0.3 #30% free space
Voidage = Initial_Voidage
Diffusivity_Toluene_in_Water = 8.6*10**(-12)
Saturated_Toluene_Conc = 515 #mg/L or g/m3
d_pore = 20e-10
tortuosity_particle = 1.3  #https://www.sciencedirect.com/science/article/pii/S0304389405007594
intraparticle_void_fraction = 0.39
intrinsic_rate_coeff = 1e-3 #M-1s-1
x_A = 0.7 #conversion of toluene
Sb = 1.0055 #Stoich Coeff Nitric Acid
S_c = 0.6535 #Stoich Coeff 2-Nitrotoluene
S_d = 0.0454 #Stoich Coeff 3-Nitrotoluene
S_e = 0.2956 #Stoich Coeff 4-Nitrotoluene
S_f = 0.0055 #Stoich Coeff Dinitrotoluene
S_g = 1.0055 #Stoich Coeff Water
gas_constant = 8.314 #J/k/mol
Temp = 333 #K
Heat_Reaction_1 = -129000 #J/mol
Conductivity = 3.63 #W/mK
#================================================≠#
#Setting up test space:
VOL_DIST_ORG = []; VOL_DIST_AQ = []; TOTAL_VOL = []; SUPER_VEL = []; RE_NO = []; RE_CHECK = []; K_TOL = []; BI_NO = []; TM_NO = []; VOL_STORE = []; LENGTH_TUBE_STORE = []; AXIAL_DIST = []; WP_CRIT = [];
EF_RANGE_STORE = [];GEF_RANGE_STORE = []; TM_RANGE_STORE = []; EF_MAX = []; GEF_MAX = []; TM_MAX = []
#===================================================#
for i in range_set:
    #================================================#
    #Calculate Concentrations in the feeds --> Know what to buy:
    print("============================================================")
    print("Mass Flowate Aqueous Feed (kg/hr) = ",Aq_feed_test_set[i])
    print("------------------------------------------------------------")
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
    Area_tube = Area_Circle(D_tube);Combined_Area = Area_tube*Number_tubes; equiv_diameter = math.sqrt(4*Combined_Area/math.pi)
    print("Diameter of a tube (m) = ", D_tube); print("Minimum diameter of a tube (m) = ", D_tube_minimum)
    if D_tube >= D_tube_minimum: print("Minimum tube diameter SATISFIED")
    else: print("***** Minimum tube diameter INVALID *****");
    print("Tube Inner Cross-sectional area (m2) = ",Area_tube);print("Total Reactor Cross-sectional area (m2) = ",Combined_Area); print("Equivalent Internal Diameter (m) = ",equiv_diameter)
    print("------------------------------------------------------------")
    
    #Calculate flow profile, composition nad features
    v_total_hour,v_frac_organic,v_frac_aq,flow_density, flow_viscosity = Vol_Flow_proportions(Mass_Toluene_Feed,Density_Toluene_Feed,Aq_feed_test_set[i],Density_Aqueous_Feed,Viscosity_Toluene,Viscosity_Aq) #m3/hr
    v_total_second =  per_hour_to_per_second(v_total_hour); print("Volumetric Flowrate (m3/s) = ",v_total_second)
    print("------------------------------------------------------------")
    
        #Calculate Area & Volume of Catalyst for further use: Assuming spherical
    Area_Cat,Vol_Cat = Area_Volume_Catalyst_Spherical(d_catalyst)
    a_s = specific_area(Voidage,Area_Cat,Vol_Cat)*v_frac_aq
    print("Volume of catalyst Bead = ",Vol_Cat);print("Surface Area of catalyst Bead = ",Area_Cat);print("Specific Surface area (m2 Cat/m3 Reactor)= ",a_s)
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
    k_toluene = MT_coeff_Surface_film_from_Sh(Sh,Diffusivity_Toluene_in_Water,d_catalyst); print("Toluene Mass Transfer Coeff, k' (m/s) = ",k_toluene)
    print("------------------------------------------------------------")
    
    #Find effective diffusion coefficient through catalyst:
    D_p_corrected = Diffusion_coeff_pore(d_pore,gas_constant,Temp,MW_Toluene); print("Knudsen Corrected Diffusion Coeff (m/s) = ",D_p_corrected)
    D_ea = Get_Effective_Diffusion_constant(intraparticle_void_fraction,D_p_corrected,tortuosity_particle); print("Effective Diffusion Coeff (m2/s) = ",D_ea)
    Bi = Biot_number(k_toluene,d_catalyst/2,D_ea); print("Biot Number = ",Bi)
    TM = Thiele_Modulus(d_catalyst,(intrinsic_rate_coeff/a_s)*Conc_NA,D_ea); print("Thiele Modulus = ",TM)
    if TM>1: print("Rate of Diffusion > Rate of Reaction") 
    elif TM==1: print("Rate of Diffusion = Rate of Reaction") 
    else: print("Rate of Diffusion < Rate of Reaction") 
    print("------------------------------------------------------------")
    effectiveness_factor,global_effectiveness_factor = Global_effectiveness_factor(TM,Bi)
    print("Catalyst Effectivness Factor = ",effectiveness_factor);print("Global Effectivness Factor = ",global_effectiveness_factor)
    print("------------------------------------------------------------")
    n_A0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Toluene_Feed)*Mass_Fraction_Toluene_Feed,MW_Toluene/1000); print("Input Molar flowrate: Toluene (mol/s) = ",n_A0)
    n_B0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Aqueous_Feed)*Mass_Fraction_NA_Aqueous_Feed,MW_NA/1000); print("Input Molar flowrate: Nitric Acid (mol/s) = ",n_B0)
    n_G0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Aqueous_Feed)*(1-Mass_Fraction_NA_Aqueous_Feed),MW_H20/1000); print("Input Molar flowrate: Water (mol/s) = ",n_G0)
    print("------------------------------------------------------------")
    #Mass Balances:
    n_A = n_A0*(1-x_A); print("Output Flowrate: Toluene (mol/s) = ",n_A)
    n_B = n_B0 - Sb*n_A0*x_A; print("Output Flowrate: Nitric Acid (mol/s) = ",n_B)
    n_C = S_c*n_A0*x_A; print("Output Flowrate: 2-Nitrotoluene (mol/s) = ",n_C)
    n_D = S_d*n_A0*x_A; print("Output Flowrate: 3-Nitrotoluene (mol/s) = ",n_D)
    n_E = S_e*n_A0*x_A; print("Output Flowrate: 4-Nitrotoluene (mol/s) = ",n_E)
    n_F = S_f*n_A0*x_A; print("Output Flowrate: Dinitrotoluene (mol/s) = ",n_F)
    n_G = n_G0 + S_g*n_A0*x_A; print("Output Flowrate: Water (mol/s) = ",n_G)
    
    print("------------------------------------------------------------")
    C_A_surface  = conc_Tol_Aq*(1-(global_effectiveness_factor*intrinsic_rate_coeff*Conc_NA)/(k_toluene*a_s))
    #ALTERNATE REACTOR SIZE -> Assuming Constant concentrations of Toluene & HNO3:
    #V_R_Const = n_A0*(1-x_A)/(global_effectiveness_factor*intrinsic_rate_coeff*Conc_NA*conc_Tol_Aq)
    V_R_Const = (-v_total_second/(global_effectiveness_factor*intrinsic_rate_coeff*conc_Tol_Aq*Sb))*math.log(1-(Sb*n_A0/n_B0)*x_A)
    print("Assuming Constant concentrations of Toluene & HNO3: Volume = ",V_R_Const)
    Length_tube_constant = V_R_Const/Combined_Area
    print("Length of tube = ",Length_tube_constant);AD_check = Axial_Dispersion_Check(Length_tube_constant,D_tube)
    print("------------------------------------------------------------")
    Beta = Prater_number(Heat_Reaction_1,D_ea,C_A_surface,Conductivity,Temp)
    print("------------------------------------------------------------")
    k_reaction = (global_effectiveness_factor*intrinsic_rate_coeff*Conc_NA)/a_s
    k_MT = k_toluene
    print("Coeff of Reaction: ",k_reaction);print("Coeff of MT Toluene: ",k_MT);
    if k_MT*10<k_reaction:
        print("Mass Transfer Limited")
    else: print("Not Mass Transfer limited")
    Weisz_platz = Weisz_Platz(effectiveness_factor,TM)
    print("------------------------------------------------------------")
    
    
    #======================================================#
    Conc_NA_end = n_B/v_total_second; print("Concentration Nitric Acid @ end (mol/m3) = ",Conc_NA_end)
    Conc_Distribution = np.linspace(Conc_NA,Conc_NA_end,100); TM_distribution = []; effectiveness_factor_distribution=[];global_effectiveness_factor_distribution=[]
    for i in range(0,len(Conc_Distribution)):
        TM_distribution.append(Thiele_Modulus(d_catalyst,(intrinsic_rate_coeff/a_s)*Conc_Distribution[i],D_ea));
        effectiveness_factor_iteration,global_effectiveness_factor_iteration = Global_effectiveness_factor(TM_distribution[i],Bi)
        effectiveness_factor_distribution.append(effectiveness_factor_iteration);global_effectiveness_factor_distribution.append(global_effectiveness_factor_iteration)
        
    ef_max = max(effectiveness_factor_distribution);ef_min = min(effectiveness_factor_distribution); ef_range = ef_max-ef_min; print(''.join("Effectiveness factor reactor range = "+str(ef_range)));
    G_ef_max = max(global_effectiveness_factor_distribution);G_ef_min = min(global_effectiveness_factor_distribution); G_ef_range = G_ef_max-G_ef_min;print("Global Effectiveness factor reactor range = ",G_ef_range)
    Tm_dist_max = max(TM_distribution);Tm_dist_min = min(TM_distribution); Tm_dist_range = Tm_dist_max-Tm_dist_min; print("Thiele Modulus reactor range = ",Tm_dist_range)
    
    if plot_TH_EF_GEF_distribution ==1:
        plt.figure()
        plt.plot(Conc_Distribution,TM_distribution);
        plt.plot(Conc_Distribution,effectiveness_factor_distribution)
        plt.plot(Conc_Distribution,global_effectiveness_factor_distribution)
        
    
    #CALLING TEST VARIABLES
    VOL_DIST_ORG.append(v_frac_organic); VOL_DIST_AQ.append(v_frac_aq); TOTAL_VOL.append(v_total_second); SUPER_VEL.append(u_super); RE_NO.append(Re); RE_CHECK.append(Re_Check); K_TOL.append(k_toluene);BI_NO.append(Bi); TM_NO.append(TM); VOL_STORE.append(V_R_Const); 
    LENGTH_TUBE_STORE.append(Length_tube_constant); AXIAL_DIST.append(AD_check); WP_CRIT.append(Weisz_platz); EF_RANGE_STORE.append(ef_range); GEF_RANGE_STORE.append(G_ef_range); TM_RANGE_STORE.append(Tm_dist_range); EF_MAX.append(ef_max); GEF_MAX.append(G_ef_max); TM_MAX.append(Tm_dist_max) 
#================================================
df = pd.DataFrame(); 
df["Mass Flowrate Aqueous"] = Aq_feed_test_set; df["Aqueous Vol Distribution"] = VOL_DIST_AQ; df["Organic Vol Distribution"] = VOL_DIST_ORG; df["Total Volumetric Flowrate"] = TOTAL_VOL; df["Reynolds Number Check"] = RE_CHECK
df["Superficial Velocity"] = SUPER_VEL; df["Reynolds Number over Sphere"] = RE_NO; df["Thiele Modulus"] = TM_NO; df["Volume of Reactor"] = VOL_STORE; df["Length of reactor"] = LENGTH_TUBE_STORE;
df["Axial Distribution Check"] = AXIAL_DIST; df["Weisz-Platz Criterion"] =WP_CRIT; df["Effectiveness Factor Range"] = EF_RANGE_STORE; df["Effectiveness Factor Max"] = EF_MAX;
df["Global Effectiveness Factor Range"] =  GEF_RANGE_STORE; df["Global Effectiveness Factor Range"] =  GEF_RANGE_STORE;
df["Thiele Modulus Range"] = TM_RANGE_STORE; df["Thiele Modulus Max"] = TM_MAX
#================================================
df = df[df["Reynolds Number Check"]==1]; df = df[df["Axial Distribution Check"]==1]
min_value_vol = df["Volume of Reactor"].min(); print("Determined Minimum Volume = ",min_value_vol); minimised_vol = df[df["Volume of Reactor"]==min_value_vol]



Re = Reynolds_J(d_catalyst,u_super,flow_density,flow_viscosity); print("Re = ",Re); Re_Check = Reynolds_Assumption_Check(Re)
j_d = J_factor_Re_Function(Re,Voidage); 
Sc = Schmidt_Number(flow_viscosity,flow_density,Diffusivity_Toluene_in_Water)
Sh = Sh_number_from_j_factor(j_d,Re,Sc)
k_toluene = MT_coeff_Surface_film_from_Sh(Sh,Diffusivity_Toluene_in_Water,d_catalyst)
#Add new diffusion coeff eqn
D_p_corrected = Diffusion_coeff_pore(d_pore,gas_constant,Temp,MW_Toluene)
D_ea = Get_Effective_Diffusion_constant(intraparticle_void_fraction,D_p_corrected,tortuosity_particle)
Bi = Biot_number(k_toluene,d_catalyst/2,D_ea)
TM = Thiele_Modulus(d_catalyst,intrinsic_rate_coeff,D_ea)

n_A0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Toluene_Feed)*Mass_Fraction_Toluene_Feed,MW_Toluene/1000)
n_B0 = Get_Molar_Flowrate(per_hour_to_per_second(Mass_Aqueous_Feed)*Mass_Fraction_NA_Aqueous_Feed,MW_NA/1000)

C_A0 = initial_conc_overall(n_A0 , v_total_second)
C_B0 = initial_conc_overall(n_B0, v_total_second)

Volume = Volume_Calc(v_total_second,global_effectiveness_factor,intrinsic_rate_coeff,Sb,n_A0,n_B0,x_A)
Length_tube = Volume/Combined_Area
print("Length of tube = ",Length_tube)

X = np.linspace(0,x_A,1000)
Vol_x = []; L_x = []

for i in range(0,1000):
    Vol_x.append(Volume_Calc(v_total_second,global_effectiveness_factor,intrinsic_rate_coeff,Sb,n_A0,n_B0,X[i]))
    L_x.append(Vol_x[i]/Combined_Area)

