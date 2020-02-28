# -*- coding: utf-8 -*-
"""Util functions"""
import numpy as np
import pandas as pd
import time

# STANDARD CURVE
def W_map_func(u):
    a = [-0.57722, 0.99999, -0.24991, 0.05519, -0.00976, 0.00108]
    b = [0.26777, 8.63476, 18.05902, 8.57333]
    c = [3.95850, 21.09965, 25.63296, 9.57332]
    if u <= 1:
        W = -np.log(u)+ a[0] + a[1]*u + a[2]*u**2 + a[3]*u**3 + a[4]*u**4 + a[5]*u**5
    else:
        W = 1/(u*np.exp(u))*(b[0]+b[1]*u+b[2]*u**2+b[3]*u**3+u**4)/(c[0]+c[1]*u+c[2]*u**2+c[3]*u**3+u**4)
    return W
    # return -0.577216-np.log(u)

def W(ua):
    """
    returns the value of well function when given
    the value of u or a series of ua by using the
    polynomial approximation equation of well
    function.
    """
    #return -0.577216-np.log(ua)

    try:
        len(ua)
    except:
        ua = np.array([ua])
    Ws = list(map(W_map_func, ua))
    return np.array(Ws)


def W_u_std(u_from, u_to, u_step):
    """
    returns the series of data of the standard line of w-(1/u)
    """
    u_std = np.arange(u_from, u_to, u_step)
    rec_u_std = 1/u_std
    W_std = W(u_std)
    return np.log10(rec_u_std), np.log10(W_std)


# READ MEASSURED DATA
def ex_data():
    df = pd.read_excel("s-t.xlsx")

    t = df["t(min)"].values
    r = df.columns[2:]
    s = np.reshape(np.transpose(df[r].values), (1, -1))[0]
    Q = df["Q(m3/h)"].values

    t_r2 = np.reshape(np.array([t/r_i**2 for r_i in r]), (1, -1))[0]

    # data clean (drop the NAN values and 0)
    cl_df = pd.DataFrame(dict(a=t_r2, b=s)).dropna()
    cl_df = cl_df.drop(cl_df[cl_df['b'].isin([0])].index)
    t_r2 = cl_df["a"].values
    s = cl_df["b"].values
    return np.log10(t_r2), np.log10(s), Q

# AUTOFIT ALGORITHM
def AutoFit(t_r2, s, algorithm):
    """
    fitting algorithm
        the loss function is RMSE (Root-Mean-Squared Error) between
        the experimental data and standard curve.
    returns the optimal x_bias, y_bias and the optimized RMSE
    """
    print("======= START AUTOMATICALLY FITTING =======")
    t1 = time.time()

    if algorithm == 0:
        """ Using Traverse algorithm to find the optimized solution """
        rmse = 100
        opt = [0, 0, rmse]
        sols = []
        for x_bias in np.arange(1e-6 , 6, 0.1):
            #print("x_bias = ", x_bias)
            for y_bias in np.arange(0, 5, 0.1):
                W_bias = (W(1/(10**(t_r2+float('%2f' %x_bias)))))
                rmse = np.sqrt(np.mean((W_bias - 10**(s+float('%.2f' %y_bias))) ** 2))
                sols.append([x_bias, y_bias, rmse])
                if rmse < opt[2]:
                    opt = [x_bias, y_bias, rmse]
            #print(x_bias, y_bias, rmse, opt)

    elif algorithm == 1:
        """ Using Mass Center algorithm to finthe optimized solution """
        rmse = 100
        opt = [0, 0, rmse]
        sols = [] ########
        strts = [] #########

        #default settigns
        total_scale = 8+2.5
        accuracy = 0.01

        series_scale = np.max(t_r2)-np.min(t_r2)
        total_step = (total_scale - series_scale)*accuracy
        strt_series = np.arange(-2.5, (total_scale - series_scale)-2.5, accuracy)
        for strt in strt_series:
            strts.append(strt) #######
            test_scale = np.arange(strt, strt + series_scale, accuracy)

            # calculate the mass center of the meassured data
            # and the theis standard line within the test scale
            mc_mesr_x, mc_mesr_y = MassCenter(t_r2, s)
            mc_std_x, mc_std_y = MassCenter(test_scale, np.log10(W(1/10**test_scale)))

            # calculate the biases of x and y
            x_bias = mc_std_x - mc_mesr_x
            y_bias = mc_std_y - mc_mesr_y

            # calculate the RMSE between meassured data and the
            # standard line data within the test scale
            W_bias = (W(1/(10**(t_r2+float('%2f' %x_bias)))))
            rmse = np.sqrt(np.mean((W_bias - 10**(s+float('%.2f' %y_bias))) ** 2))

            # find the optimal solutions
            sols.append([x_bias, y_bias, rmse])
            if rmse < opt[2]:
                opt = [x_bias, y_bias, rmse]
                #print(x_bias, y_bias, rmse, opt)

    elif algorithm == 2:
        """ Using Slope algorithm to find the optimal solution """
        # calculate the slope of the meassured data
        slp_mesr, itrcpt_mesr = LinFit(t_r2, s)

        #default settigns
        total_scale = 8+2.5
        accuracy = 0.01

        series_scale = np.max(t_r2)-np.min(t_r2)
        strt_series = np.arange(-2.5, (total_scale - series_scale)-2.5, accuracy)
        for strt in strt_series:
            test_scale = np.arange(strt, strt + series_scale, accuracy)

            # calculate the slope of the theis standard line within the test scale
            slp_std, itrcpt_std = LinFit(test_scale, np.log10(W(1/10**test_scale)))

            if abs(slp_mesr - slp_std) < 1e-2:
                break
        # calculate the biases of x and y
        x_bias = strt - np.min(t_r2)
        y_bias = slp_mesr * x_bias + (itrcpt_std - itrcpt_mesr)

        # calculate the RMSE between meassured data and the
        # standard line data within the test scale
        W_bias = (W(1/(10**(t_r2+float('%2f' %x_bias)))))
        rmse = np.sqrt(np.mean((W_bias - 10**(s+float('%.2f' %y_bias))) ** 2))

        opt = [x_bias, y_bias, rmse]
    
    print("Time costs: %.4f secs" %(time.time() - t1))
    return opt


def MassCenter(xs, ys):
    """
    returns the coordinates of the mass center of a two-dimensional data
    """
    return np.mean(xs), np.mean(ys)


def LinFit(xs, ys):
    """
    returns the slope of the linear fitting function of xs and ys
    """
    return np.polyfit(xs, ys, 1)


# RUN AUTOFIT
def calResult(t_r2, s, Q, x_bias, y_bias, num=2):
    """
    returns the W 1/u s t_r2 T S
    """
    W_p = 1
    rac_u_p = 10
    s_p = 10**(np.log10(W_p)-y_bias)
    t_r2_p = 10**(np.log10(rac_u_p)-x_bias)

    Q_p = Q[num]

    T = Q_p*W_p*24/(4*np.pi*s_p)
    S = 4*T*(t_r2_p)/(rac_u_p*60*24)
    return W_p, rac_u_p, s_p, t_r2_p, T, S


def run_AutoFit(algorithm):
    """
    run the AutoFit algorithm
    """
    t_r2, s, Q = ex_data()
    x_bias, y_bias, RMSE = AutoFit(t_r2, s, algorithm)
    W_p, rac_u_p, s_p, t_r2_p, T, S = calResult(t_r2, s, Q, x_bias, y_bias, num=2)

    # print the results
    print("======== Auto Fit Results ========")
    print('    W = %.2f  rac_u = %.2f  s = %.2f  t_r2 = %.4f' %(W_p,rac_u_p,s_p,t_r2_p))
    print('    T = %.2f m^3/d  S = %.6f  x_bias = %.2f  y_bias = %.2f' %(T, S, x_bias, y_bias))
    print("======== Statistics of Fitting ========")
    print("    RMSE = %.4F" %(RMSE))
    return t_r2, s, Q, x_bias, y_bias, RMSE, W_p, rac_u_p, s_p, t_r2_p, T, S
