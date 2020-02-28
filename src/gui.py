# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox, RadioButtons

from utils import *

class GUI():
    def __init__(self):
        self.axcolor = 'lightgrey'
        self.fig = plt.figure(figsize=(8, 6))
        self.updating_status = False
        self.main()
    
    def main(self):
        self.calc_infos()
        self.init_figure()
        self.init_control_zone()
        self.init_result_zone()
        self.add_events()
        plt.show()

    def calc_infos(self, algorithm=0):
        self.t_r2, self.s, self.Q, self.x_bias, self.y_bias, self.RMSE, self.W_p, self.rac_u_p, self.s_p, self.t_r2_p, self.T, self.S = run_AutoFit(algorithm)
    
    def init_figure(self):
        ######## Plot Zone ########
        ax = self.fig.add_subplot(1,1,1)
        plt.subplots_adjust(bottom=0.42, left=0.18, right=0.65)

        rec_u, W_std = W_u_std(1e-8, 10, 0.01)
        self.l, = plt.plot(self.t_r2+self.x_bias, self.s+self.y_bias, 's', color='maroon',
                    markersize=5, label='measurement data')
        plt.plot(rec_u, W_std, color='k',
                linewidth=0.85, label='theis type-curve')
        plt.xlabel(r'lg$(\frac{1}{u})$ ', fontsize=10)
        plt.ylabel(r'lg$(W)$', fontsize=10)
        plt.xlim(np.min(self.t_r2+self.x_bias)-1, np.max(self.t_r2+self.x_bias)+1)
        plt.ylim(np.min(self.s+self.y_bias)-1, np.max(self.s+self.y_bias)+0.5)
        plt.grid(True, ls='--', lw=.5, c='k', alpha=0.6)
        plt.legend(fontsize=8)
        plt.suptitle("Theis Type-curve Fitting")

    def init_control_zone(self):
        ######## Control Zone ########
        contrl_ax = plt.axes([0.05, 0.05, 0.90, 0.27])
        plt.text(0.035, 0.85, 'Control Zone')
        contrl_ax.set_xticks([])
        contrl_ax.set_yticks([])

        # slider
        ax_xbias = plt.axes([0.14, 0.21, 0.35, 0.03], facecolor=self.axcolor)
        ax_ybias = plt.axes([0.14, 0.16, 0.35, 0.03], facecolor=self.axcolor)

        self.s_xbias = Slider(ax_xbias, 'x_bias', 0.0, 5, valinit=self.x_bias)
        self.s_ybias = Slider(ax_ybias, 'y_bias', 0.0, 5, valinit=self.y_bias)

        # buttons
        aftax = plt.axes([0.8, 0.06, 0.1, 0.04])
        self.reset_button = Button(aftax, 'RESET', color=self.axcolor, hovercolor='mistyrose')

        x_bias_dec_tax = plt.axes([0.57, 0.21, 0.03, 0.03])
        self.x_bias_dec_button = Button(x_bias_dec_tax, '-', color=self.axcolor, hovercolor='mistyrose')

        x_bias_inc_tax = plt.axes([0.72, 0.21, 0.03, 0.03])
        self.x_bias_inc_button = Button(x_bias_inc_tax, '+', color=self.axcolor, hovercolor='mistyrose')

        y_bias_dec_tax = plt.axes([0.57, 0.16, 0.03, 0.03])
        self.y_bias_dec_button = Button(y_bias_dec_tax, '-', color=self.axcolor, hovercolor='mistyrose')

        y_bias_inc_tax = plt.axes([0.72, 0.16, 0.03, 0.03])
        self.y_bias_inc_button = Button(y_bias_inc_tax, '+', color=self.axcolor, hovercolor='mistyrose')

        # textbox
        xtextax = plt.axes([0.61, 0.21, 0.1, 0.03])
        ytextax = plt.axes([0.61, 0.16,  0.1, 0.03])
        self.t_xbias = TextBox(xtextax, ' ', str(float('%.2f' %self.x_bias)))
        self.t_ybias = TextBox(ytextax, ' ', str(float('%.2f' %self.y_bias)))

        # radio button
        rax = plt.axes([0.78, 0.12, 0.15, 0.17], facecolor=self.axcolor)
        plt.title('AutoFit Algorithms', fontsize=8)
        self.algos = np.array(['Traverse', 'Mass Center', 'Slope'])
        self.radio = RadioButtons(rax, self.algos, active=0)

    def init_result_zone(self, ):
        ######## Result Zone ########
        re_ax = plt.axes([0.70, 0.37, 0.20, 0.50])
        plt.text(0.035, 0.95, 'Results Zone')
        re_ax.set_xticks([])
        re_ax.set_yticks([])

        self.W_l = plt.text(0.2, 0.85, r'$W: %.2f $' %(self.W_p))
        self.u_l = plt.text(0.2, 0.75, r'$1/u: %.2f $' %(self.rac_u_p))
        self.s_l = plt.text(0.2, 0.65, r'$s: %.2f $' %(self.s_p))
        self.tr2_l = plt.text(0.2, 0.55, r'$t/r^2: %.4f $' %(self.t_r2_p))
        self.T_l = plt.text(0.2, 0.45, r'$T: %.2f m^3/d$' %(self.T))
        self.S_l = plt.text(0.2, 0.35, r'$S: %.6f $' %(self.S))
        plt.text(0.035, 0.2, 'Statistics of Fitting')
        self.RMSE_l = plt.text(0.2, 0.1, r'$RMSE: %.4f$' %(self.RMSE))

    def add_events(self):
        self.s_xbias.on_changed(self.update_x)
        self.s_ybias.on_changed(self.update_y)
        self.reset_button.on_clicked(self.reset)
        self.x_bias_dec_button.on_clicked(self.dec_x_bias)
        self.x_bias_inc_button.on_clicked(self.inc_x_bias)
        self.y_bias_dec_button.on_clicked(self.dec_y_bias)
        self.y_bias_inc_button.on_clicked(self.inc_y_bias)
        self.t_xbias.on_submit(self.submit_x)
        self.t_ybias.on_submit(self.submit_y)
        self.radio.on_clicked(self.chg_algorithm)

    def chg_result(self, t_r2, s, x_bias, y_bias):
        """ change the results in the results zone """
        self.W_p, self.rac_u_p, self.s_p, self.t_r2_p, self.T, self.S = calResult(self.t_r2, self.s, self.Q, self.x_bias, self.y_bias, num=2)
        self.W_bias = (W(1/(10**(t_r2+float('%2f' %x_bias)))))
        self.RMSE = np.sqrt(np.mean((self.W_bias - 10**(s+float('%.2f' %y_bias))) ** 2))

        self.W_l.set_text(r'$W: %.2f $' %(self.W_p))
        self.u_l.set_text(r'$1/u: %.2f $' %(self.rac_u_p))
        self.s_l.set_text(r'$s: %.2f $' %(self.s_p))
        self.tr2_l.set_text(r'$t/r^2: %.4f $' %(self.t_r2_p))
        self.T_l.set_text(r'$T: %.2f m^3/d$' %(self.T))
        self.S_l.set_text(r'$S: %.6f $' %(self.S))
        self.RMSE_l.set_text(r'$RMSE: %.4f$' %(self.RMSE))



    def update_x(self, x_bias):
        """update all widgets with new x_bias"""
        # 防抖
        if self.updating_status:
            return 

        x_bias = float('%.2f' %x_bias)
        if x_bias == self.x_bias:
            return

        self.updating_status = True
        self.x_bias = x_bias
        self.s_xbias.set_val(float('%.2f' %self.x_bias))
        self.t_xbias.set_val(float('%.2f' %self.x_bias))
        self.l.set_xdata(self.t_r2+self.x_bias)

        self.fig.canvas.draw_idle()
        self.chg_result(self.t_r2, self.s, self.x_bias, self.y_bias)
        self.updating_status = False

    def update_y(self, y_bias):
        """update all widgets with new y_bias"""
        # 防抖
        if self.updating_status:
            return 

        y_bias = float('%.2f' %y_bias)
        if y_bias == self.y_bias:
            return

        self.updating_status = True
        self.y_bias = y_bias
        self.s_ybias.set_val(float('%.2f' %self.y_bias))
        self.t_ybias.set_val(float('%.2f' %self.y_bias))
        self.l.set_ydata(self.s+self.y_bias)

        self.fig.canvas.draw_idle()
        self.chg_result(self.t_r2, self.s, self.x_bias, self.y_bias)
        self.updating_status = False

    def reset(self, event):
        """ reset the figure """
        self.radio.set_active(0)

    def dec_x_bias(self, event):
        """ decrease the x_bias using '-' button """
        self.update_x(self.x_bias-0.01)

    def inc_x_bias(self, event):
        """ increase the x_bias using '+' button """
        self.update_x(self.x_bias+0.01)

    def dec_y_bias(self, event):
        """ decrease the y_bias using '-' button """
        self.update_y(self.y_bias-0.01)

    def inc_y_bias(self, event):
        """ increase the y_bias using '+' button """
        self.update_y(self.y_bias+0.01)

    def submit_x(self, text):
        """ updates the figure using textbox (x_bias)"""
        self.update_x(float(text))

    def submit_y(self, text):
        """ updates the figure using textbox (y_bias)"""
        self.update_y(float(text))

    def chg_algorithm(self, label):
        algorithm = np.argwhere(self.algos==label)[0][0]
        # t_r2, s, Q, x_bias, y_bias, RMSE, W_p, rac_u_p, s_p, t_r2_p, T, S = run_AutoFit(algorithm)
        self.calc_infos(algorithm)
        self.update_x(self.x_bias)
        self.update_y(self.y_bias)




