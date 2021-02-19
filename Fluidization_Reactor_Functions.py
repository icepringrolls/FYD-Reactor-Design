#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:13:34 2021

@author: henryorlebar
"""
def bubble_velocity(u_super,u_min_fluid,tube_diameter,bubble_diameter):
    if tube_diameter<0.1: alpha = 0.64
    elif tube_diameter>=0.1 and tube_diameter<1: alpha = 1.6*tube_diameter**0.4
    elif tube_diameter>=1: alpha = 1.6
    
    U_br = alpha*(bubble_diameter*9.81)**(1/2)
    U_b = (u_super-u_min_fluid)+U_br
    return U_b

def get_parameters(bubble_velocity,u_super,u_min_fluid,bed_voidage):
    f_b = (u_super-u_min_fluid)/bubble_velocity
    voidage_fluid = 1-(1-f_b)*(1-bed_voidage)
    f_e = voidage_fluid -f_b
    return f_b,voidage_fluid, f_e