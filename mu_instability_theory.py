import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
from math import erf, exp

def W_calc(delta):
    output = np.zeros((delta.size,3), dtype = float)
    for i in range(0,delta.size):
        output[i,0] = 0.5*(1 + erf(delta[i]/np.sqrt(2)))
        output[i,1] = 0.5*((np.sqrt(2/np.pi)*np.exp(-(delta[i]**2)/2))+delta[i]*(1 + erf(delta[i]/np.sqrt(2))))
        output[i,2] = output[i,0] + delta[i]*output[i,1]
    return output

def Variance(w,big_gamma,delta,r): #This function calculates the variance of the set of the data in terms of the w variables
    var_array = np.zeros((delta.size), dtype = float)
    for i in range(delta.size):
        temp = 1+r*(w[i,1]**2)/w[i,2]
        var_array[i] = w[i,2]*temp/((w[i,2]*temp+w[i,0]*big_gamma)**2)
    return var_array

def Mu(w,big_gamma,delta,r,gamma):
    Mu_array = np.zeros((delta.size), dtype = float)
    for i in range(delta.size):
        temp = 1+r*(w[i,1]**2)/w[i,2]
        mu = delta[i]*w[i,2]*temp - gamma*w[i,0]*w[i,1]
        mu /= w[i,1]*(w[i,0]*big_gamma+w[i,2]*temp)
        Mu_array[i] = mu
    return Mu_array

def delta_function(delta, r):
    omega_0 = 0.5*(1+erf(delta/(np.sqrt(2))))
    omega_1 = np.exp(-(delta**2)/2)/(np.sqrt(np.pi*2)) + delta*omega_0
    omega_2 = omega_0 + delta*omega_1
    # omega_2 + r*(omega_1**2)-omega_0
    return delta+r*omega_1

def delta_solver (row):
    delta = scipy.optimize.bisect(delta_function, -5, 0, args = row)
    return delta

def variance_opper (delta, big_gamma,r):
    phi = 0.5*(1+erf(delta/np.sqrt(2)))
    variance = 1/(phi*((1+big_gamma)**2))
    return variance 



def vary_gamma(number_lines, w_array0, delta0):
    gamma_array = np.array([-.8,0,0.8])
    ROWS = 0.0
    SMALL_GAMMA = 0.0
    plt.title(r'Instability Plot in Parameter Space for Varied $\Gamma$', fontsize = 16)
    cmap = plt.get_cmap('jet_r')

    for i, BIG_GAMMA in enumerate(gamma_array):
        delta_opper = delta_solver(ROWS)
        opper_point = variance_opper(delta_opper, BIG_GAMMA, ROWS)

        mu_array = Mu(w_array0,BIG_GAMMA,delta0,ROWS,SMALL_GAMMA)
        variance_array = Variance(w_array0,BIG_GAMMA,delta0,ROWS)
        plot_1 = np.column_stack((mu_array, variance_array))
        indx1 = np.where(plot_1[:,1] > opper_point)[0]
        plot_1 = np.delete(plot_1, indx1, axis = 0)
        
        indx2 = np.argmax(plot_1[:,1], axis = 0)
        max_mu = plot_1[indx2,0]
        num_points = 1000
        x_plot = np.linspace(-10,max_mu, num_points)
        y_plot = np.repeat(opper_point,num_points)
        plt.plot(plot_1[:,0], plot_1[:,1], label = r'$\Gamma$= %.2f' %BIG_GAMMA, color = cmap(float(i)/len(gamma_array)))
        plt.plot(x_plot, y_plot,'--', color = cmap(float(i)/len(gamma_array)))
        plt.yscale("log")
        plt.ylim(10**(-2),10**(2))


def vary_l_gamma(number_lines, w_array0, delta0):
    MAX_lgamma = 2
    lgamma_array = np.arange(1,2.5,0.5)
    ROWS = 1
    BIG_GAMMA = 0
    plt.title(r'Instability Plot in Parameter Space for Varied $\gamma$', fontsize=16)
    cmap = plt.get_cmap('jet_r')

    delta_opper = delta_solver(ROWS)
    opper_point = variance_opper(delta_opper, BIG_GAMMA, ROWS)
    variance_array = Variance(w_array0,BIG_GAMMA,delta0,ROWS)
    # indx2 = np.argmax(plot_1[:,1], axis = 0)
    # max_mu = plot_1[indx2,0]
    num_points = 1000
    x_plot = np.linspace(-3.5,1, num_points)
    y_plot = np.repeat(opper_point,num_points)

    for i, SMALL_GAMMA in enumerate(lgamma_array):
        mu_array = Mu(w_array0,BIG_GAMMA,delta0,ROWS,SMALL_GAMMA)
        plot_1 = np.column_stack((mu_array, variance_array))
        indx1 = np.where(plot_1[:,1] > opper_point)[0]
        plot_1 = np.delete(plot_1, indx1, axis = 0)
        # print("GAMMA:", SMALL_GAMMA, "/", opper_point)
        # idx = np.where((mu_array >-0.01)&(mu_array <0.01))[0]
        # Print the corresponding value of y1
        # if len(idx) > 0:
        #     for i in idx:
        #         print(f'mu={mu_array[i]:.2f}, var={variance_array[i]:.2f}')
        
        
        plt.plot(plot_1[:,0], plot_1[:,1], label = r'$\gamma$= %.1f' %SMALL_GAMMA, color = cmap(float(i)/len(lgamma_array)))
        plt.plot(x_plot, y_plot,'--', color = cmap(float(i)/len(lgamma_array)))

def vary_r(number_lines, w_array0, delta0):
    MAX_R = 5
    SMALL_GAMMA = 0.0
    BIG_GAMMA = 0.0
    rows_array = np.arange(0,MAX_R + MAX_R/number_lines, MAX_R/number_lines)
    plt.title(r'Instability Plot in Parameter Space for Varied $r$', fontsize=16)
    cmap = plt.get_cmap('jet_r')

    for i, ROWS in enumerate(rows_array):
        delta_opper = delta_solver(ROWS)
        opper_point = variance_opper(delta_opper, BIG_GAMMA, ROWS)

        mu_array = Mu(w_array0,BIG_GAMMA,delta0,ROWS,SMALL_GAMMA)
        variance_array = Variance(w_array0,BIG_GAMMA,delta0,ROWS)
        plot_1 = np.column_stack((mu_array, variance_array))
        indx1 = np.where(plot_1[:,1] > opper_point)[0]
        plot_1 = np.delete(plot_1, indx1, axis = 0)
        # print("r:", ROWS, "/", opper_point)
        # idx = np.where((mu_array > -1.01) & (mu_array < -0.99))[0]
        # Print the corresponding value of y1
        # if len(idx) > 0:
        #     for i in idx:
        #         print(f'mu={mu_array[i]:.2f}, var={variance_array[i]:.1f}')
        
        indx2 = np.argmax(plot_1[:,1], axis = 0)
        max_mu = plot_1[indx2,0]
        num_points = 1000
        x_plot = np.linspace(-30,max_mu, num_points)
        y_plot = np.repeat(opper_point,num_points)

        plt.plot(plot_1[:,0], plot_1[:,1], label = 'r= %.d' %ROWS, color = cmap(float(i)/len(rows_array)))
        plt.plot(x_plot, y_plot,'--', color = cmap(float(i)/len(rows_array)))
        plt.fill_between(x, y, 0, color='blue', alpha=.1)
        # plt.yscale("log")
        # plt.ylim(10**(-2),10**(2))
    

number_lines = 5
delta0 = np.linspace(0,50, num = 100000)
w_array0 = W_calc(delta0)
fig = plt.figure()
plt.xlabel(r"$\mu$",fontsize=16)
plt.ylabel(r"$\sigma^2$",fontsize=16)

#vary_r(number_lines, w_array0, delta0)
#vary_l_gamma(number_lines, w_array0, delta0)
vary_gamma(number_lines, w_array0, delta0)

# plt.ylim(0.001,100)
plt.xlim(-1,2)
plt.legend(fontsize = 14)
plt.show()