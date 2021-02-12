#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 11:32:46 2021

@author: henryorlebar
"""
def Concentration_from_MassComp(Molar_Mass,Mass_Fraction,Mass_Density):
    Concentration = Mass_Fraction*(Mass_Density/Molar_Mass)
    print(Concentration)
    return Concentration

def Get_Molar_Flowrate(Mass_Flowrate,molar_weight):
    Molar_Flowrate = Mass_Flowrate/molar_weight
    return Molar_Flowrate


def per_hour_to_per_second(per_hour_variable):
    per_second = per_hour_variable/3600
    return per_second

def Area_Circle(diameter):
    import math
    A = math.pi*0.25*(diameter**2)
    return A

def Superficial_velocity(volumetric_flowrate,Wetted_Area):
    u_super = volumetric_flowrate/Wetted_Area
    return u_super

def Reynolds_J(dp,u_superficial,rho,mu):
    Re = (dp*u_superficial*rho)/mu
    return Re

def Reynolds_Assumption_Check(Re):
    if Re>=10 and Re<=2000:
        print("Reynolds Limits Satisified")
        Re_Check =1 
    else:
        print("Reynolds Limits NOT Satisified")
        Re_Check = 0
    return Re_Check
  
def Particle_Voidage(Volume_Reactor, Density_Catalyst, Mass_Catalyst):
    Voidage = (1-(Mass_Catalyst/Density_Catalyst))/Volume_Reactor
    return Voidage

def J_factor_Re_Function(Re,voidage):
    j_d = 0.4548/(voidage*(Re**(0.4069)))
    return j_d

def Schmidt_Number(mu,rho,diffusivity):
    Sc = mu/(rho*diffusivity)
    return Sc

def Vol_Flow_proportions(Mass_flowrate_Organic,Density_Organic,Mass_flowrate_Aq,Density_Aq,Viscosity_Organic,Viscosity_Aq):
    v_organic = Mass_flowrate_Organic/Density_Organic; v_aq = Mass_flowrate_Aq/Density_Aq
    v_total = v_organic+v_aq
    v_frac_organic = v_organic/v_total; v_frac_aq = v_aq/v_total
    density_weighted = v_frac_organic*Density_Organic+v_frac_aq*Density_Aq
    viscosity_weighted = v_frac_organic*Viscosity_Organic+v_frac_aq*Viscosity_Aq
    print("Vol% Organic = ",v_frac_organic);print("Vol% Aqueous = ",v_frac_aq); print("Total Volumetric Flowrate (m3/hr) = ", v_total)
    return v_total,v_frac_organic,v_frac_aq,density_weighted,viscosity_weighted

def Sh_number_from_j_factor(j_d,Re,Sc):
    Sh = j_d*Re*(Sc**(1/3))
    return Sh

def MT_coeff_Surface_film_from_Sh(Sh,Diffusivity,dp):
    k = Sh*Diffusivity/dp
    return k

def Get_Effective_Diffusion_constant(voidage_particle,diffusion_coeff,tortuosity):
    Dea = (voidage_particle*diffusion_coeff)/tortuosity
    return Dea

def Biot_number(kmc,particle_radius,D_ea):
    Bi = kmc*particle_radius/D_ea
    return Bi

def Thiele_Modulus(Diameter_Particle,Intrinsic_Rate,D_ea):
    import math
    TM = (Diameter_Particle/6)*math.sqrt(Intrinsic_Rate/D_ea)
    return TM

def Global_effectiveness_factor(Thiele_modulus,Biot_Number):
    import numpy as np
    effectiveness_factor = np.tanh(Thiele_modulus)/Thiele_modulus
    global_effectiveness_factor = ((1/effectiveness_factor)+((Thiele_modulus**2)/Biot_Number))**(-1)
    return effectiveness_factor,global_effectiveness_factor

def Volume_Integral_Component(n_A0,Sb,n_B0,x):
    import math
    A = Sb*n_A0; B = n_B0
    Comp = (math.log(abs(A*x-B))-math.log(abs(x-1)))/(B-A)
    return Comp

def Volume_Calc(v_T,E_Factor,rate_coef,Sb,n_A0,n_B0,x_A):
    Coef = ((v_T**2)/(E_Factor*rate_coef))
    Comp1 = Volume_Integral_Component(n_A0,Sb,n_B0,x_A)
    Comp2 = Volume_Integral_Component(n_A0,Sb,n_B0,0)
    Volume = Coef*(Comp1-Comp2)
    return Volume
    

def Get_effectiveness_factor(intrinsic_rate_constant, conc_nitronium,p_d, gas_constant,temperature, Toluene_MolarMass, d_pore, tortuosity,intraparticle_void_fraction):
    import math
    import numpy as np
    kv = intrinsic_rate_constant * conc_nitronium
    D_p = (d_pore/3) * math.sqrt(8*gas_constant * temperature/(math.pi * Toluene_MolarMass))
    D_p_corrected = D_p / 20           #knudesen overestimated diffusion by a factor of 20 https://www.osti.gov/servlets/purl/1424576                                                              #https://www.sciencedirect.com/science/article/pii/S0304389405007594
    D_ea = (D_p_corrected*intraparticle_void_fraction)/tortuosity
    thiele_modulus = 8/3*p_d*math.sqrt((kv/(D_ea)))
    effectiveness_factor = np.tanh(thiele_modulus)/thiele_modulus
    return thiele_modulus, effectiveness_factor, D_ea


def get_global_effectiveness_factor(effectiveness_factor, thiele_modulus, kmc, particle_radius, D_ea):
    biot_number = kmc*particle_radius/D_ea
    global_effectiveness_factor = 1/(1/effectiveness_factor + thiele_modulus**2/biot_number)
    return global_effectiveness_factor