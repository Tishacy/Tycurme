# -*- coding: utf-8 -*-
"""
=================
    Tycurme
Type-Curve Method
=================

A Graphic Method for Solving Hydrogeological Parameters
by Matching the Measured Curve with the Theoretical Curve
in the Pumping Test.

Author: TIANShuo
    from Ocean University of China
    Institute of Environmental Science and Engineering
Date:   2018/6/12
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.widgets import Slider, Button, TextBox, RadioButtons


# STANDARD CURVE
def W(ua):
    """
    returns the value of well function when given
    the value of u or a series of ua by using the
    polynomial approximation equation of well
    function.
    """
    #return -0.577216-np.log(ua)
    a = [-0.57722, 0.99999, -0.24991, 0.05519, -0.00976, 0.00108]
    b = [0.26777, 8.63476, 18.05902, 8.57333]
    c = [3.95850, 21.09965, 25.63296, 9.57332]
    try:
        len(ua)
    except:
        ua = np.array([ua])
    Ws = []
    for u in ua:
        if u <= 1:
            W = -np.log(u)+ a[0] + a[1]*u + a[2]*u**2 + a[3]*u**3 + a[4]*u**4 + a[5]*u**5
        else:
            W = 1/(u*np.exp(u))*(b[0]+b[1]*u+b[2]*u**2+b[3]*u**3+u**4)/(c[0]+c[1]*u+c[2]*u**2+c[3]*u**3+u**4)        
        Ws.append(W)
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
    
    t_r2 = []
    for r_i in r:
        t_r2.append(t/r_i**2)
    t_r2 = np.reshape(np.array(t_r2), (1, -1))[0]
    
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
        the loss function is RMSE(Root-Mean-Squared Error) between
        the experimental data and standard curve.
    returns the optimal x_bias, y_bias and the optimized RMSE
    """
    print("======= START AUTOMATICALLY FITTING =======")

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
    
##    # print the results
##    print("======== Auto Fit Results ========")
##    print('    W = %.2f  rac_u = %.2f  s = %.2f  t_r2 = %.4f' %(W_p,rac_u_p,s_p,t_r2_p))
##    print('    T = %.2f m^3/d  S = %.6f  x_bias = %.2f  y_bias = %.2f' %(T, S, x_bias, y_bias))
##    print("======== Statistics of Fitting ========")
##    print("    RMSE = %.4F" %(RMSE))
    return t_r2, s, Q, x_bias, y_bias, RMSE, W_p, rac_u_p, s_p, t_r2_p, T, S



# GUI
def GUI(algorithm = 0):
    """
    User interface of type curve method.
    """
    ######## Create a GUI and get data ########
    axcolor = 'lightgrey'
    fig = plt.figure(figsize=(8, 8))
    t_r2, s, Q, x_bias, y_bias, RMSE, W_p, rac_u_p, s_p, t_r2_p, T, S = run_AutoFit(algorithm)


    ######## Funtions in GUI ########
    def update(val):
        """ updates the figure using slider"""
        x_bias = s_xbias.val
        y_bias = s_ybias.val
        t_xbias.set_val(float('%.2f' %x_bias))
        t_ybias.set_val(float('%.2f' %y_bias))
        l.set_xdata(t_r2+x_bias)
        l.set_ydata(s+y_bias)

        chg_result(t_r2, s, x_bias, y_bias)
        fig.canvas.draw_idle()

    def reset(event):
        """ reset the figure """
        radio.set_active(0)
        s_xbias.reset()
        s_ybias.reset()

    def dec_x_bias(event):
        """ decrease the x_bias using '-' button """
        x_bias = s_xbias.val
        t_xbias.set_val(float('%.2f' %(x_bias-0.01)))

    def inc_x_bias(event):
        """ increase the x_bias using '+' button """
        x_bias = s_xbias.val
        t_xbias.set_val(float('%.2f' %(x_bias+0.01)))

    def dec_y_bias(event):
        """ decrease the y_bias using '-' button """
        y_bias = s_ybias.val
        t_ybias.set_val(float('%.2f' %(y_bias-0.01)))

    def inc_y_bias(event):
        """ increase the y_bias using '+' button """
        y_bias = s_ybias.val
        t_ybias.set_val(float('%.2f' %(y_bias+0.01)))
    
    def submit_x(text):
        """ updates the figure using textbox (x_bias)"""
        x_bias = float('%.2f' %float(text))
        s_xbias.set_val(x_bias)

    def submit_y(text):
        """ updates the figure using textbox (y_bias)"""
        y_bias = float('%.2f' %float(text))
        s_ybias.set_val(y_bias)

    def chg_result(t_r2, s, x_bias, y_bias):
        """ change the results in the results zone """
        W_p, rac_u_p, s_p, t_r2_p, T, S = calResult(t_r2, s, Q, x_bias, y_bias, num=2)
        W_bias = (W(1/(10**(t_r2+float('%2f' %x_bias)))))
        RMSE = np.sqrt(np.mean((W_bias - 10**(s+float('%.2f' %y_bias))) ** 2))
        
        W_l.set_text(r'$W: %.2f $' %(W_p))
        u_l.set_text(r'$1/u: %.2f $' %(rac_u_p))
        s_l.set_text(r'$s: %.2f $' %(s_p))
        tr2_l.set_text(r'$t/r^2: %.4f $' %(t_r2_p))
        T_l.set_text(r'$T: %.2f m^3/d$' %(T))
        S_l.set_text(r'$S: %.6f $' %(S))
        RMSE_l.set_text(r'$RMSE: %.4f$' %(RMSE))

    def chg_algorithm(label):
        algorithm = np.argwhere(algos==label)[0][0]
        t_r2, s, Q, x_bias, y_bias, RMSE, W_p, rac_u_p, s_p, t_r2_p, T, S = run_AutoFit(algorithm)
        t_xbias.set_val(float('%.2f' %x_bias))
        t_ybias.set_val(float('%.2f' %y_bias))
        l.set_xdata(t_r2+x_bias)
        l.set_ydata(s+y_bias)
        chg_result(t_r2, s, x_bias, y_bias)
        fig.canvas.draw_idle()

    ######## Plot Zone ########
    ax = fig.add_subplot(1,1,1)
    plt.subplots_adjust(bottom=0.42, left=0.18, right=0.65)

    rec_u, W_std = W_u_std(1e-8, 10, 1e-4)
    l, = plt.plot(t_r2+x_bias, s+y_bias, 's', color='maroon',
                  markersize=5, label='measurement data')
    plt.plot(rec_u, W_std, color='k',
             linewidth=0.85, label='theis type-curve')
    plt.xlabel(r'lg$(\frac{1}{u})$ ', fontsize=10)
    plt.ylabel(r'lg$(W)$', fontsize=10)
    plt.xlim(np.min(t_r2+x_bias)-1, np.max(t_r2+x_bias)+1)
    plt.ylim(np.min(s+y_bias)-1, np.max(s+y_bias)+0.5)
    plt.grid(True, ls='--', lw=.5, c='k', alpha=0.6)
    plt.legend(fontsize=8)
    plt.suptitle("Theis Type-curve Fitting")

    ######## Control Zone ########
    contrl_ax = plt.axes([0.05, 0.05, 0.90, 0.27])
    plt.text(0.035, 0.85, 'Control Zone')
    contrl_ax.set_xticks([])
    contrl_ax.set_yticks([])
    
    # slider
    ax_xbias = plt.axes([0.14, 0.21, 0.35, 0.03], facecolor=axcolor)
    ax_ybias = plt.axes([0.14, 0.16, 0.35, 0.03], facecolor=axcolor)

    s_xbias = Slider(ax_xbias, 'x_bias', 0.0, 5, valinit=x_bias)
    s_ybias = Slider(ax_ybias, 'y_bias', 0.0, 5, valinit=y_bias)

    s_xbias.on_changed(update)
    s_ybias.on_changed(update)

    # button
    aftax = plt.axes([0.8, 0.06, 0.1, 0.04])
    button = Button(aftax, 'RESET', color=axcolor, hovercolor='mistyrose')
    button.on_clicked(reset)
    
    x_bias_dec_tax = plt.axes([0.57, 0.21, 0.03, 0.03])
    x_bias_dec_button = Button(x_bias_dec_tax, '-', color=axcolor, hovercolor='mistyrose')
    x_bias_dec_button.on_clicked(dec_x_bias)

    x_bias_inc_tax = plt.axes([0.72, 0.21, 0.03, 0.03])
    x_bias_inc_button = Button(x_bias_inc_tax, '+', color=axcolor, hovercolor='mistyrose')
    x_bias_inc_button.on_clicked(inc_x_bias)

    y_bias_dec_tax = plt.axes([0.57, 0.16, 0.03, 0.03])
    y_bias_dec_button = Button(y_bias_dec_tax, '-', color=axcolor, hovercolor='mistyrose')
    y_bias_dec_button.on_clicked(dec_y_bias)

    y_bias_inc_tax = plt.axes([0.72, 0.16, 0.03, 0.03])
    y_bias_inc_button = Button(y_bias_inc_tax, '+', color=axcolor, hovercolor='mistyrose')
    y_bias_inc_button.on_clicked(inc_y_bias)

    # textbox
    xtextax = plt.axes([0.61, 0.21, 0.1, 0.03])
    ytextax = plt.axes([0.61, 0.16,  0.1, 0.03])
    t_xbias = TextBox(xtextax, ' ', str(float('%.2f' %x_bias)))
    t_ybias = TextBox(ytextax, ' ', str(float('%.2f' %y_bias)))
    t_xbias.on_submit(submit_x)
    t_ybias.on_submit(submit_y)

    # radio button
    rax = plt.axes([0.78, 0.12, 0.15, 0.17], facecolor=axcolor)
    plt.title('AutoFit Algorithms', fontsize=8)
    algos = np.array(['Traverse', 'Mass Center', 'Slope'])
    radio = RadioButtons(rax, algos,
                         active=algorithm)
    radio.on_clicked(chg_algorithm)
    
    ######## Result Zone ########
    re_ax = plt.axes([0.70, 0.37, 0.20, 0.50])
    plt.text(0.035, 0.95, 'Results Zone')
    re_ax.set_xticks([])
    re_ax.set_yticks([])
    
    W_l = plt.text(0.2, 0.85, r'$W: %.2f $' %(W_p))
    u_l = plt.text(0.2, 0.75, r'$1/u: %.2f $' %(rac_u_p))
    s_l = plt.text(0.2, 0.65, r'$s: %.2f $' %(s_p))
    tr2_l = plt.text(0.2, 0.55, r'$t/r^2: %.4f $' %(t_r2_p))
    T_l = plt.text(0.2, 0.45, r'$T: %.2f m^3/d$' %(T))
    S_l = plt.text(0.2, 0.35, r'$S: %.6f $' %(S))
    plt.text(0.035, 0.2, 'Statistics of Fitting')
    RMSE_l = plt.text(0.2, 0.1, r'$RMSE: %.4f$' %(RMSE))

    plt.show()
    
    
if __name__=="__main__":
    GUI()
    
    
