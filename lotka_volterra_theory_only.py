# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 11:16:51 2023

@author: joshh
Semester 2 extension of MPhys project where general correlations are introduced


"""
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.legend_handler import HandlerPathCollection
from matplotlib.text import Text
import matplotlib.cm as cm

cmap = cm.get_cmap('viridis')

data = pd.read_csv('row_bg_-0.8_lg_0.csv', sep=",", skiprows= 4)
ROWS_DOT = np.sort(data['ROW'].unique())

# Round them to 1 decimal place and convert to string
ROWS_rounded = [f'{ROW:.1f}' for ROW in ROWS_DOT]

def W_calc(delta):
    output = np.zeros((delta.size,3), dtype = float)
    for i in range(0,delta.size):
        output[i,0] = 0.5*(1 + math.erf(delta[i]/np.sqrt(2)))
        output[i,1] = 0.5*((np.sqrt(2/np.pi)*np.exp(-(delta[i]**2)/2))+delta[i]*(1 + math.erf(delta[i]/np.sqrt(2))))
        output[i,2] = output[i,0] + delta[i]*output[i,1]
    return output

def Variance(w,big_gamma,delta,r): #This function calculates the variance of the set of the data in terms of the w variables
    var_array = np.zeros((delta.size), dtype = float)
    for i in range(0,delta.size):
        var_array[i] = w[i,2]/((w[i,2] + ((big_gamma * w[i,0])/(1 + r*(w[i,1]**2/w[i,2]))))**2 * (1 + r*(w[i,1]**2/w[i,2])))
        #var_array[i] = 1/(w[i,2] - r*(w[i,1]**2)) #REDUCED VARIANCE

    return var_array

def Biomass(w,big_gamma,delta,r,gamma, MU):
    bio_array = np.zeros((delta.size), dtype = float)
    for i in range(0,delta.size):
        b = (w[i,2] * delta[i])*(1 + r*(w[i,1]**2/w[i,2])) - gamma*w[i,0]*w[i,1]
        b/= (w[i,1] * (w[i,2]*(1 + r*(w[i,1]**2/w[i,2])) + (big_gamma * w[i,0])))
        b-= MU
        bio_array[i] = 1/b
    return bio_array

def Chi(w,gamma,delta,r):
    chi_array = np.zeros((delta.size), dtype = float)
    for i in range(0,delta.size):
        chi_array[i] = w[i,0] + ((gamma/(1 - r*(w[i,1]**2/w[i,2]))) * ((w[i,0]**2)/w[i,2]))
    return chi_array

def Q(w, bio,delta):
    q_array = np.zeros((delta.size), dtype = float)
    for i in range(0,delta.size):
        q_array[i] = bio**2 * ((w[i,1])**2/w[i,2])
    return q_array 


MU = -1
gam_arr = np.arange(0,1.0,0.2) # ONE INDEX
BIG_GAMMA = -0.8
GAMMA = 0
row_arr = np.arange(0,3,0.5)
delta0 = np.linspace(0,500, num = 100000)
w_array0 = W_calc(delta0)
norm = plt.Normalize(min(row_arr), max(row_arr))

plt.title("Average Abundance of Species for Varied $r$",fontsize = 16)
plt.xlabel(r"$\sigma^2$",fontsize = 14)
plt.ylabel(r"M",fontsize = 14)

for ROWS in row_arr:
    plt.plot(Variance(w_array0,BIG_GAMMA,delta0,ROWS), Biomass(w_array0, BIG_GAMMA,delta0, ROWS, GAMMA, MU), color=cmap(norm(ROWS)))
scatter = sns.scatterplot(x='VARIANCE', y='BIOMASS', hue='ROW', palette='viridis', data=data)
handles, labels = scatter.get_legend_handles_labels()
scatter.legend(handles=handles[1:], labels=ROWS_rounded, title=r'$r$', handler_map={type(handles[0]): HandlerPathCollection()})

plt.xlim(0,10)
# plt.ylim(0,1.005)
plt.show()
plt.legend()
plt.show()


